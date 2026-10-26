# Deploying MICRO-MINE on GitHub + Vercel (free)

This repository contains **two ways of running the web page**:

| | Vercel (serverless, "online lite") | Docker / server (full) |
|---|---|---|
| Entry point | `api/index.py` + `public/` | `micromine-web` / `gunicorn micromine.web.app:app` |
| How it works | one request: upload -> analysis -> results (nothing is stored) | background jobs, results kept 24 h |
| Upload limit | ~4.2 MB **compressed** per analysis (the browser compresses the FASTA files; about 50 chloroplast genomes of 155 kb fit) | 10 MB per file, 50 files |
| Time limit | 60 s per analysis (`vercel.json`) | 15 min |
| PDF figures | no (SVG only) | yes, if `cairosvg` is installed |
| Cost | free (Vercel Hobby plan, personal / non-commercial use) | your server |

Measured locally: 50 chloroplast genomes of 155 kb -> 2.4 MB compressed upload, ~11-14 s, ~1 MB of results.
Vercel CPUs can be slower, so keep an eye on the 60 s limit for very large sets.

## 1. Test locally first

    pip install -e ".[web]"
    python scripts/dev_vercel.py          # http://127.0.0.1:3000  (same behaviour as on Vercel)

## 2. Put the code on GitHub

    git init
    git add .
    git commit -m "MICRO-MINE v0.13"
    git branch -M main
    git remote add origin https://github.com/wagner512able/micro-mine.git
    git push -u origin main

(or create the repository on github.com, click *Add file > Upload files* and drag the contents of this folder;
keep the folder structure and do not skip `public/`, `api/`, `micromine/` and `vercel.json`).

## 3. Deploy on Vercel

1. Sign in at vercel.com with your GitHub account.
2. *Add New... > Project*, pick the `micro-mine` repository, framework preset **Other**, leave the build settings empty.
3. (Optional) *Environment Variables*: `MICROMINE_REPO_URL` = `https://github.com/wagner512able/micro-mine`
   (adds a "Source code on GitHub" link in the footer).
4. Click **Deploy**. Every `git push` to `main` redeploys automatically.

If you edit the page (`micromine/web/static/`), run `python scripts/sync_public.py` before committing, so `public/`
(the folder Vercel serves) is updated.

## 4. If something fails

* *Build error about dependencies*: check that `requirements.txt` is in the repository root.
* *"Function timed out"*: use fewer genomes or the full (Docker) version; the limit is `maxDuration` in `vercel.json`
  (Hobby plan: up to 60 s, or more if Fluid compute is enabled in the project settings).
* *HTTP 413*: the compressed upload exceeds ~4.5 MB (Vercel limit). Use fewer genomes or run locally.
* Logs: Vercel dashboard > your project > *Logs*.

## 5. Full version (long analyses, 10 MB files)

Use the `Dockerfile` on any container host (Render, Fly.io, a university server, Hugging Face Spaces with the Docker SDK).
It runs `gunicorn -w 1 --threads 8` on port 8000. Keep a single worker process: the concurrency control is per process.
