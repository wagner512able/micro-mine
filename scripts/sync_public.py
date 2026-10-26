"""Copy the web page (index.html + static files) into public/, which Vercel serves as static files.

Run it after editing micromine/web/static/:  python scripts/sync_public.py
"""
import os
import shutil

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = os.path.join(root, "micromine", "web", "static")
dst = os.path.join(root, "public")
shutil.rmtree(dst, ignore_errors=True)
os.makedirs(os.path.join(dst, "static"))
shutil.copy(os.path.join(src, "index.html"), os.path.join(dst, "index.html"))
for f in os.listdir(src):
    if f != "index.html":
        shutil.copy(os.path.join(src, f), os.path.join(dst, "static", f))
print("public/ updated:", sorted(os.listdir(dst)), sorted(os.listdir(os.path.join(dst, "static"))))
