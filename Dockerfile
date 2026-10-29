FROM python:3.12-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337

WORKDIR /processor

# Pinned Python dependencies (umap-learn with numba/pynndescent, matplotlib),
# plus the Lambda runtime interface client for running as a Lambda function.
COPY requirements.txt /processor/requirements.txt
RUN pip install --no-cache-dir -r /processor/requirements.txt awslambdaric==3.1.1 simplejson==4.2.0 snapshot-restore-py==1.0.0 \
 && find /usr/local/lib/python3.12 -name "tests" -type d -prune -exec rm -rf {} + \
 && find /usr/local/lib/python3.12 -name "__pycache__" -type d -prune -exec rm -rf {} +

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

# numba caches compiled functions and matplotlib its font list; both need a
# writable directory, and on Lambda only /tmp is writable.
ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

# Runs main() once as a container, or serves processor.handler as a Lambda
# function (see processor/entry.py).
ENTRYPOINT ["python3", "-m", "processor.entry"]
