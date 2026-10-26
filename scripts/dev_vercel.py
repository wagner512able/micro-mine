"""Run the Vercel (serverless) version locally, exactly as it will behave online.

    python scripts/dev_vercel.py [port]      ->  http://127.0.0.1:3000
"""
import os
import sys

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, root)
from flask import send_from_directory  # noqa: E402
from micromine.web.serverless import app  # noqa: E402

PUBLIC = os.path.join(root, "public")


@app.get("/")
def index():
    return send_from_directory(PUBLIC, "index.html")


@app.get("/static/<path:name>")
def static_files(name):
    return send_from_directory(os.path.join(PUBLIC, "static"), name)


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=int(sys.argv[1]) if len(sys.argv) > 1 else 3000)
