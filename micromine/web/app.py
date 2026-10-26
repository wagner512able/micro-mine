"""MICRO-MINE web server (Flask).

Run for development:   python -m micromine.web.app
Run in production:     gunicorn -w 1 --threads 8 -b 0.0.0.0:8000 micromine.web.app:app

Limits (environment variables):
  MICROMINE_MAX_FILE_MB      max size of each FASTA (default 10)
  MICROMINE_MAX_GENOME_BP    max genome length per file (default 1,000,000; organelles only)
  MICROMINE_MAX_CONCURRENT   analyses running at the same time (default 2)
  MICROMINE_JOB_TTL_HOURS    hours before results are deleted (default 24)
  MICROMINE_JOBS_PER_HOUR    max analyses per client IP per hour (default 20)
  MICROMINE_WORKDIR          where jobs are stored (default: system temp dir)
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import uuid
import zipfile

from flask import Flask, abort, jsonify, request, send_from_directory

from .. import __version__
from ..core import CLASS_NAMES, DEFAULT_MIN_REPEATS
from .common import (FASTA_EXT as _FASTA_EXT, MAX_CONTIGS, MAX_FILES, MAX_GENOME_BP, MIN_FILES, REPO_URL,
                     parse_params as _parse_params, validate_uploads)

MAX_FILE_MB = int(os.environ.get("MICROMINE_MAX_FILE_MB", "10"))
MAX_CONCURRENT = int(os.environ.get("MICROMINE_MAX_CONCURRENT", "2"))
JOB_TTL_H = float(os.environ.get("MICROMINE_JOB_TTL_HOURS", "24"))
JOBS_PER_HOUR = int(os.environ.get("MICROMINE_JOBS_PER_HOUR", "20"))
TIMEOUT_S = 900
WORKDIR = os.environ.get("MICROMINE_WORKDIR") or os.path.join(tempfile.gettempdir(), "micromine_jobs")

_HERE = os.path.dirname(os.path.abspath(__file__))
_PKG_ROOT = os.path.dirname(os.path.dirname(_HERE))
_JOB_RE = re.compile(r"^[0-9a-f]{32}$")
_SEM = threading.Semaphore(MAX_CONCURRENT)
_hits: dict[str, list] = {}
_hits_lock = threading.Lock()

app = Flask(__name__, static_folder=os.path.join(_HERE, "static"), static_url_path="/static")
app.json.sort_keys = False      # keep class / variation order in the JSON responses
app.config["MAX_CONTENT_LENGTH"] = (MAX_FILES * MAX_FILE_MB + 5) * 1024 * 1024
os.makedirs(WORKDIR, exist_ok=True)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _job_dir(job_id: str) -> str:
    if not _JOB_RE.match(job_id):
        abort(404)
    d = os.path.join(WORKDIR, job_id)
    if not os.path.isdir(d):
        abort(404)
    return d


def _set_status(d: str, **kw) -> None:
    path = os.path.join(d, "status.json")
    try:
        with open(path) as fh:
            st = json.load(fh)
    except (OSError, ValueError):
        st = {}
    st.update(kw)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(st, fh)
    os.replace(tmp, path)


def _get_status(d: str) -> dict:
    try:
        with open(os.path.join(d, "status.json")) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {"status": "error", "message": "Job state not found."}


def _rate_ok(ip: str) -> bool:
    now = time.time()
    with _hits_lock:
        recent = [t for t in _hits.get(ip, []) if now - t < 3600]
        if len(recent) >= JOBS_PER_HOUR:
            _hits[ip] = recent
            return False
        recent.append(now)
        _hits[ip] = recent
    return True


def _run_job(job_id: str) -> None:
    d = os.path.join(WORKDIR, job_id)
    with _SEM:
        _set_status(d, status="running", started=time.time())
        ok, msg = False, "The analysis failed."
        env = dict(os.environ)
        env["PYTHONPATH"] = _PKG_ROOT + os.pathsep + env.get("PYTHONPATH", "")
        try:
            with open(os.path.join(d, "job.log"), "w") as log:
                proc = subprocess.run([sys.executable, "-m", "micromine.web.worker", d], cwd=d,
                                      stdout=log, stderr=subprocess.STDOUT, timeout=TIMEOUT_S, env=env)
            ok = proc.returncode == 0 and os.path.exists(os.path.join(d, "summary.json"))
            if not ok:
                with open(os.path.join(d, "job.log")) as fh:
                    lines = [l.strip() for l in fh if l.strip()]
                err = [l for l in lines if l.startswith("Error:")]
                msg = (err[-1][6:].strip() if err else "The analysis failed unexpectedly.")
        except subprocess.TimeoutExpired:
            msg = f"The analysis exceeded the {TIMEOUT_S // 60}-minute limit."
        except OSError as e:
            msg = f"Could not start the analysis: {e}"
        if ok:
            out = os.path.join(d, "output")
            with zipfile.ZipFile(os.path.join(d, "micromine_results.zip"), "w", zipfile.ZIP_DEFLATED) as z:
                for root, _, files in os.walk(out):
                    for f in files:
                        full = os.path.join(root, f)
                        z.write(full, os.path.relpath(full, out))
            _set_status(d, status="done", finished=time.time())
        else:
            _set_status(d, status="error", message=msg, finished=time.time())


def _cleanup_loop() -> None:
    while True:
        try:
            cutoff = time.time() - JOB_TTL_H * 3600
            for name in os.listdir(WORKDIR):
                p = os.path.join(WORKDIR, name)
                if _JOB_RE.match(name) and os.path.isdir(p) and os.path.getmtime(p) < cutoff:
                    shutil.rmtree(p, ignore_errors=True)
        except OSError:
            pass
        time.sleep(1800)


threading.Thread(target=_cleanup_loop, daemon=True).start()


# --------------------------------------------------------------------------- #
# routes
# --------------------------------------------------------------------------- #
@app.after_request
def _headers(resp):
    resp.headers["X-Content-Type-Options"] = "nosniff"
    if request.path.startswith("/results/") or request.path.startswith("/download/"):
        resp.headers["Content-Security-Policy"] = "default-src 'none'; style-src 'unsafe-inline'; img-src data:"
    else:
        resp.headers["Content-Security-Policy"] = ("default-src 'self'; img-src 'self' data:; "
                                                   "style-src 'self' 'unsafe-inline'; script-src 'self' 'unsafe-inline'")
    return resp


@app.get("/")
def index():
    return send_from_directory(app.static_folder, "index.html")


@app.get("/api/config")
def config():
    return jsonify(version=__version__, max_files=MAX_FILES, min_files=MIN_FILES,
                   max_file_mb=MAX_FILE_MB, max_genome_bp=MAX_GENOME_BP, job_ttl_hours=JOB_TTL_H,
                   classes=[[name, DEFAULT_MIN_REPEATS[k]] for k, name in CLASS_NAMES.items()],
                   mode="server", repo_url=REPO_URL)


@app.post("/api/run")
def start():
    if not _rate_ok(request.remote_addr or "?"):
        return jsonify(error=f"Rate limit reached ({JOBS_PER_HOUR} analyses per hour). Try again later."), 429
    files = request.files.getlist("files")
    if not MIN_FILES <= len(files) <= MAX_FILES:
        return jsonify(error=f"Upload between {MIN_FILES} and {MAX_FILES} FASTA files (one file = SSR annotation only)."), 400
    try:
        labels = json.loads(request.form.get("labels", "[]"))
        if not isinstance(labels, list) or len(labels) != len(files):
            raise ValueError("Species labels do not match the uploaded files.")
        params = _parse_params(request.form.get("params", "{}"))
        items = []
        for f, lab in zip(files, labels):
            items.append((f.filename, f.read(MAX_FILE_MB * 1024 * 1024 + 1), lab))
        cleaned, problems = validate_uploads(items, MAX_FILE_MB * 1024 * 1024)
        if problems:
            return jsonify(error=f"{len(problems)} file(s) rejected: " + " | ".join(p["message"] for p in problems),
                           problems=problems), 400
    except ValueError as e:
        return jsonify(error=str(e)), 400

    job_id = uuid.uuid4().hex
    d = os.path.join(WORKDIR, job_id)
    os.makedirs(os.path.join(d, "input"))
    for fname, text in cleaned:
        with open(os.path.join(d, "input", fname), "w") as fh:
            fh.write(text)
    with open(os.path.join(d, "params.json"), "w") as fh:
        json.dump(params, fh)
    _set_status(d, status="queued", created=time.time(), n_files=len(cleaned))
    threading.Thread(target=_run_job, args=(job_id,), daemon=True).start()
    return jsonify(job=job_id)


@app.get("/api/status/<job_id>")
def status(job_id):
    d = _job_dir(job_id)
    st = _get_status(d)
    resp = {"status": st.get("status"), "message": st.get("message", "")}
    try:
        with open(os.path.join(d, "job.log")) as fh:
            resp["log"] = fh.read()[-3000:].replace(WORKDIR, "<jobs>")
    except OSError:
        resp["log"] = ""
    if st.get("status") == "done":
        with open(os.path.join(d, "summary.json")) as fh:
            resp["summary"] = json.load(fh)
    return jsonify(resp)


@app.get("/results/<job_id>/<path:name>")
def result_file(job_id, name):
    d = _job_dir(job_id)
    return send_from_directory(os.path.join(d, "output"), name)


@app.get("/download/<job_id>.zip")
def download(job_id):
    d = _job_dir(job_id)
    return send_from_directory(d, "micromine_results.zip", as_attachment=True,
                               download_name="micromine_results.zip")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="micromine-web", description="MICRO-MINE web server (development).")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8000)
    a = ap.parse_args(argv)
    app.run(host=a.host, port=a.port, threaded=True)


if __name__ == "__main__":
    main()
