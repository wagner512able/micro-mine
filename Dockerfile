FROM python:3.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends libcairo2 && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY . /app
RUN pip install --no-cache-dir ".[web,pdf]"
ENV THANNYA_WORKDIR=/data/jobs
RUN useradd -m thannya && mkdir -p /data/jobs && chown -R thannya /data
USER thannya
EXPOSE 8000
# One worker process (job state is per-process semaphores); threads handle concurrent requests.
CMD ["gunicorn", "-w", "1", "--threads", "8", "--timeout", "300", "-b", "0.0.0.0:8000", "thannya.web.app:app"]
