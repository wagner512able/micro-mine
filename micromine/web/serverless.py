"""MICRO-MINE for serverless hosting (Vercel): stateless and synchronous.

Nothing is kept between requests: the browser sends the genomes (gzip-compressed JSON), the
function analyses them in a temporary folder and returns the results in the same response
(summary, figures and files as base64), then deletes everything.

Vercel limits to keep in mind (Hobby plan): request and response bodies of 4.5 MB, a maximum function
duration (set in vercel.json) and a read-only file system except /tmp.
"""
from __future__ import annotations

import base64
import io
import json
import os
import shutil
import tempfile
import zipfile
import zlib

from flask import Flask, jsonify, request

from .. import __version__
from ..cli import run
from ..core import CLASS_NAMES, DEFAULT_MIN_REPEATS
from .common import MAX_FILES, MAX_GENOME_BP, MIN_FILES, REPO_URL, build_summary, parse_params, validate_uploads

MAX_FILE_MB = 10
MAX_PAYLOAD_MB = 4.2                 # compressed upload accepted by the page (Vercel allows 4.5 MB)
MAX_RAW = 60 * 1024 * 1024           # decompressed size accepted by the function (zip-bomb guard)
RESPONSE_LIMIT = int(4.2 * 1024 * 1024)
MIMES = {".xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
         ".svg": "image/svg+xml", ".zip": "application/zip", ".gff3": "text/plain", ".fasta": "text/plain",
         ".pdf": "application/pdf"}

app = Flask(__name__, static_folder=None)
app.config["MAX_CONTENT_LENGTH"] = 6 * 1024 * 1024


@app.after_request
def _headers(resp):
    resp.headers["X-Content-Type-Options"] = "nosniff"
    resp.headers["Cache-Control"] = "no-store"
    return resp


@app.get("/api/config")
def config():
    return jsonify(mode="serverless", version=__version__, max_files=MAX_FILES, min_files=MIN_FILES,
                   max_file_mb=MAX_FILE_MB, max_genome_bp=MAX_GENOME_BP, max_payload_mb=MAX_PAYLOAD_MB,
                   job_ttl_hours=0, repo_url=REPO_URL,
                   classes=[[name, DEFAULT_MIN_REPEATS[k]] for k, name in CLASS_NAMES.items()])


def _read_payload():
    raw = request.get_data()
    if "gzip" in (request.content_type or ""):
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            raw = d.decompress(raw, MAX_RAW)
        except zlib.error:
            raise ValueError("The upload could not be decompressed.")
        if d.unconsumed_tail:
            raise ValueError("The upload is too large for the online version.")
    try:
        payload = json.loads(raw)
    except ValueError:
        raise ValueError("Invalid request.")
    if not isinstance(payload, dict) or not isinstance(payload.get("files"), list):
        raise ValueError("Invalid request.")
    return payload


def _pack(out_dir: str, zip_bytes: bytes, files: list):
    """Build the 'downloads' dict, dropping optional items until the response fits the limit."""
    def b64(data):
        return base64.b64encode(data).decode("ascii")

    def entry(name, data):
        return {"mime": MIMES.get(os.path.splitext(name)[1].lower(), "application/octet-stream"), "b64": b64(data)}

    zip_entry = {"micromine_results.zip": entry("micromine_results.zip", zip_bytes)}
    singles = {f: entry(f, open(os.path.join(out_dir, f), "rb").read())
               for f in files if not f.lower().endswith(".pdf")}
    figures = {f: v for f, v in singles.items() if f.lower().endswith(".svg")}

    def size(d):
        return sum(len(v["b64"]) for v in d.values())

    for candidate in ({**zip_entry, **singles}, {**zip_entry, **figures}, zip_entry):
        if size(candidate) <= RESPONSE_LIMIT:
            return candidate
    raise ValueError("The results are too large for the online version (limit of about 4 MB). Use fewer "
                     "genomes or fewer FASTA outputs, or run MICRO-MINE locally.")


@app.post("/api/analyze")
def analyze():
    tmp = None
    try:
        try:
            payload = _read_payload()
            params = parse_params(payload.get("params") or {})
            raw_files = payload["files"]
            if not MIN_FILES <= len(raw_files) <= MAX_FILES:
                raise ValueError(f"Upload between {MIN_FILES} and {MAX_FILES} FASTA files "
                                 f"(one file = SSR annotation only).")
            items = [(str(f.get("name", "")), str(f.get("text", "")).encode("utf-8"), f.get("label", ""))
                     for f in raw_files if isinstance(f, dict)]
            if len(items) != len(raw_files):
                raise ValueError("Invalid request.")
        except ValueError as e:
            return jsonify(error=str(e)), 400
        cleaned, problems = validate_uploads(items, MAX_FILE_MB * 1024 * 1024)
        if problems:
            return jsonify(error=f"{len(problems)} file(s) rejected: " + " | ".join(p["message"] for p in problems),
                           problems=problems), 400

        tmp = tempfile.mkdtemp(prefix="micromine_")
        in_dir, out_dir = os.path.join(tmp, "input"), os.path.join(tmp, "output")
        os.makedirs(in_dir)
        for fname, text in cleaned:
            with open(os.path.join(in_dir, fname), "w") as fh:
                fh.write(text)
        log = []
        try:
            genomes, clusters, ssrs = run(
                in_dir, out_dir, {int(k): int(v) for k, v in params["min_repeats"].items()},
                flank=params["flank"], anchor=params["anchor"], similarity=params["similarity"],
                max_files=MAX_FILES, threads=1, pdf=False, relax=params["relax"],
                regions=[tuple(r) for r in params["regions"]], log=log.append,
                input_label="uploaded files (web)", fasta_layout=tuple(params["fasta_layout"]),
                ssr_case=params["ssr_case"])
        except (ValueError, OSError) as e:
            return jsonify(error=str(e)), 400
        summary = build_summary(genomes, clusters, ssrs, out_dir)
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for root, _, fs in os.walk(out_dir):
                for f in fs:
                    full = os.path.join(root, f)
                    z.write(full, os.path.relpath(full, out_dir))
        try:
            downloads = _pack(out_dir, buf.getvalue(), summary["files"])
        except ValueError as e:
            return jsonify(error=str(e)), 413
        return jsonify(summary=summary, downloads=downloads, log="\n".join(log))
    except Exception:                      # never leak a stack trace to the browser
        app.logger.exception("analysis failed")
        return jsonify(error="The analysis failed unexpectedly. Please check the files and try again."), 500
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)
