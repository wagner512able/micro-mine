"""Serverless (Vercel) app test with the Flask test client (python tests/test_serverless.py)."""
import base64, gzip, io, json, os, sys, tempfile, zipfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from make_demo_data import build
from micromine.web.serverless import app

d = tempfile.mkdtemp(); build(d)
files = [{"name": n, "label": n.split("_", 1)[1][:-6], "text": open(os.path.join(d, n)).read()} for n in sorted(os.listdir(d))]
c = app.test_client()
cfg = c.get("/api/config").get_json(); assert cfg["mode"] == "serverless" and cfg["max_files"] == 50


def post(payload, gz=True):
    body = json.dumps(payload).encode()
    if gz:
        return c.post("/api/analyze", data=gzip.compress(body), content_type="application/gzip")
    return c.post("/api/analyze", data=body, content_type="application/json")

# validation
assert post({"files": [], "params": {}}).status_code == 400
assert post({"files": files[:1] + [{"name": "x.fasta", "text": ">x\nACGTRY\n", "label": "bad"}], "params": {}}).status_code == 400
r = post({"files": [{"name": "a.fasta", "text": ">x\nACGTRY\n"}, {"name": "b.fasta", "text": ">y\nAC-GU\n"}] + files[:1], "params": {}})
j = r.get_json(); assert r.status_code == 400 and len(j["problems"]) == 2 and j["problems"][0]["bases"] == {"R": 1, "Y": 1}, j
assert post({"files": files, "params": {"flank": 75, "anchor": 100}}).status_code == 400
assert c.post("/api/analyze", data=b"not gzip", content_type="application/gzip").status_code == 400
assert c.post("/api/analyze", data=gzip.compress(b"x" * 100), content_type="application/gzip").status_code == 400
import zlib
bomb = zlib.compress(b"0" * (70 * 1024 * 1024), 9)       # decompresses to 70 MB (> limit) from a few kB
c2 = zlib.compressobj(9, zlib.DEFLATED, 31); bomb = c2.compress(b"0" * (70 * 1024 * 1024)) + c2.flush()
r = c.post("/api/analyze", data=bomb, content_type="application/gzip"); assert r.status_code == 400 and "too large" in r.get_json()["error"]

# comparison, gzip and plain JSON bodies give the same answer, and nothing stays behind
before = set(os.listdir(tempfile.gettempdir()))
r = post({"files": files, "params": {"regions": "1000-8000", "flank": 75, "ssr_case": "upper"}}); assert r.status_code == 200, r.data
j = r.get_json(); s = j["summary"]
assert s["mode"] == "comparison" and s["loci"] == 14 and s["exclusive_real"] == 1 and len(s["genomes"]) == 4
assert {"micromine_results.zip", "micromine_report.xlsx", "micromine_comparison.svg"} <= set(j["downloads"])
assert not any(n.endswith(".pdf") for n in j["downloads"])
z = zipfile.ZipFile(io.BytesIO(base64.b64decode(j["downloads"]["micromine_results.zip"]["b64"])))
assert "micromine_report.xlsx" in z.namelist() and any(n.startswith("fasta_per_locus/") for n in z.namelist())
assert "flank=75/75" in z.read("fasta_per_species/01_ancestral_SSRs.fasta").decode()
assert base64.b64decode(j["downloads"]["micromine_comparison.svg"]["b64"]).startswith(b"<svg")
assert post({"files": files, "params": {}}, gz=False).get_json()["summary"]["loci"] == 14
assert not [x for x in set(os.listdir(tempfile.gettempdir())) - before if x.startswith("micromine_")], "temp files left behind"
# single file: annotation
r = post({"files": files[:1], "params": {}}); a = r.get_json()["summary"]
assert r.status_code == 200 and a["mode"] == "annotation" and "micromine_SSRs.gff3" in r.get_json()["downloads"]
# the response never exceeds the limit: the page falls back to zip-only, then to an error
import micromine.web.serverless as sv
zlen = len(j["downloads"]["micromine_results.zip"]["b64"])
sv.RESPONSE_LIMIT = zlen + 100
r = post({"files": files, "params": {"regions": "1000-8000", "flank": 75, "ssr_case": "upper"}})
assert r.status_code == 200 and list(r.get_json()["downloads"]) == ["micromine_results.zip"], list(r.get_json()["downloads"])
sv.RESPONSE_LIMIT = 100
r = post({"files": files, "params": {}}); assert r.status_code == 413 and "too large" in r.get_json()["error"]
print("OK: serverless app validated.")
