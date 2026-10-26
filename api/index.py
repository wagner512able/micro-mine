"""Vercel entry point: exposes the Flask app (WSGI) found in micromine.web.serverless."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from micromine.web.serverless import app  # noqa: E402,F401
