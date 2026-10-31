"""Web app test with the Flask test client (python tests/test_web.py)."""
import io, json, os, sys, tempfile, time, zipfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
os.environ["MICROMINE_WORKDIR"] = tempfile.mkdtemp()
from make_demo_data import build
from micromine.web.app import app

d = tempfile.mkdtemp(); build(d)
names = sorted(os.listdir(d))
c = app.test_client()
assert c.get("/").status_code == 200 and b"MICRO-MINE" in c.get("/").data and b"wagner2@egresso.ufg.br" in c.get("/").data and c.get("/static/logo.jpg").status_code == 200
cfg = c.get("/api/config").get_json(); assert cfg["max_files"] == 50


def post(files, labels=None, params=None):
    data = {"files": [(io.BytesIO(raw), fn) for fn, raw in files],
            "labels": json.dumps(labels or [f"sp{i}" for i in range(len(files))]),
            "params": json.dumps(params or {"regions": "1000-8000", "ssr_case": "upper", "fasta_layout": ["locus"]})}
    return c.post("/api/run", data=data, content_type="multipart/form-data")

good = [(n, open(os.path.join(d, n), "rb").read()) for n in names]
# validation
assert post([]).status_code == 400                                             # needs >= 1 file
assert post([("a.txt", b">x\nACGT\n")] * 2).status_code == 400                 # extension
assert post([("a.fasta", b">x\nACGT123\n")] * 2).status_code == 400            # invalid characters
rb = post([("a.fasta", b">x\nACGTRY\n"), ("b.fasta", b">y\nACGT-U\n"), good[0]])           # IUPAC / gap / U rejected
jb = rb.get_json(); assert rb.status_code == 400 and len(jb["problems"]) == 2, jb
assert jb["problems"][0]["bases"] == {"R": 1, "Y": 1} and jb["problems"][0]["file"] == "a.fasta" and "b.fasta" in jb["error"]
assert post([("a.fasta", b">x\n" + b"A" * 2_100_000 + b"\n")] * 2).status_code == 400   # too long
assert post(good, params={"mono": 1}).status_code == 400                       # min repeats of 1
assert post(good, params={"regions": "9-3"}).status_code == 400
assert post(good, params={"flank": 75, "anchor": 100}).status_code == 400      # anchor > flank
assert post(good, params={"fasta_layout": ["x"]}).status_code == 400
assert post(good, params={"ssr_case": "mixed"}).status_code == 400
assert c.get("/api/status/" + "0" * 32).status_code == 404
assert c.get("/results/../etc/passwd").status_code in (404, 400)

# full run (a '=' in a contig name must be neutralized)
good[0] = (good[0][0], good[0][1].replace(b">01_ancestral_plastome", b">=cmd|x"))
lines = good[1][1].split(b"\n"); lines[9] = lines[9][:20] + b"NNNNN" + lines[9][25:]       # N bases are accepted
good[1] = (good[1][0], b"\n".join(lines))
r = post(good, labels=["ancestral", "snps flank", "loss+excl", "rotated rc"]); assert r.status_code == 200, r.data
job = r.get_json()["job"]
for _ in range(120):
    st = c.get(f"/api/status/{job}").get_json()
    if st["status"] in ("done", "error"): break
    time.sleep(1)
assert st["status"] == "done", st
s = st["summary"]; assert s["loci"] >= 13 and s["exclusive_real"] == 1 and len(s["genomes"]) == 4
assert s["genomes"][1]["name"] == "snps_flank"
assert s["compromised_loci"] >= 1 and s["genomes"][1]["n_bases"] == 5 and s["genomes"][1]["compromised"] >= 1
assert "micromine_report.xlsx" in s["files"] and "micromine_region_1000-8000.svg" in s["files"]
assert c.get(f"/results/{job}/micromine_comparison.svg").status_code == 200
z = c.get(f"/download/{job}.zip"); assert z.status_code == 200
names_z = zipfile.ZipFile(io.BytesIO(z.data)).namelist(); assert "micromine_report.xlsx" in names_z
assert any(n.startswith("fasta_per_locus/") for n in names_z) and not any(n.startswith("fasta_per_species/") for n in names_z)
assert s["fasta_per_locus"] == s["loci"] and s["fasta_per_species"] == 0
import openpyxl
wb = openpyxl.load_workbook(io.BytesIO(zipfile.ZipFile(io.BytesIO(z.data)).read("micromine_report.xlsx")))
contigs = {r[4] for r in wb["SSR_comparison"].iter_rows(min_row=2, values_only=True)}
assert not any(str(x).startswith("=") for x in contigs), contigs
# single file: annotation mode with a 75 bp flank
r = post(good[:1], labels=["ancestral"], params={"flank": 75}); assert r.status_code == 200, r.data
job1 = r.get_json()["job"]
for _ in range(60):
    st1 = c.get(f"/api/status/{job1}").get_json()
    if st1["status"] in ("done", "error"): break
    time.sleep(1)
assert st1["status"] == "done", st1
a = st1["summary"]; assert a["mode"] == "annotation" and a["ssrs"] >= 12 and set(a["files"]) >= {"micromine_SSRs.gff3", "micromine_report.xlsx", "micromine_annotation_map.svg", "SSRs_with_flanks.fasta"}
z1 = zipfile.ZipFile(io.BytesIO(c.get(f"/download/{job1}.zip").data))
gff = z1.read("micromine_SSRs.gff3").decode().strip().split("\n")
assert gff[0] == "##gff-version 3" and sum(1 for l in gff if not l.startswith("#")) == a["ssrs"]
fa = z1.read("SSRs_with_flanks.fasta").decode(); assert "flank=75/75" in fa and fa.count(">") == a["ssrs"]
print("OK: web app validated; job finished with", s["loci"], "loci.")
