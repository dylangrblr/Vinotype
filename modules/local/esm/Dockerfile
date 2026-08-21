# ESM2 variant-effect predictor — CPU base image.
# For GPU support, rebuild with:
#   docker build \
#     --build-arg TORCH_INDEX_URL=https://download.pytorch.org/whl/cu121 \
#     -t ghcr.io/dylangrblr/vinotype/esm_vep:1.0.0-cuda .
#
# Build (CPU):
#   cd modules/local/esm
#   docker build -t ghcr.io/dylangrblr/vinotype/esm_vep:1.0.0 .
# Push:
#   docker push ghcr.io/dylangrblr/vinotype/esm_vep:1.0.0

FROM python:3.11-slim

ARG TORCH_INDEX_URL=https://download.pytorch.org/whl/cpu

RUN apt-get update \
    && apt-get install -y --no-install-recommends gcc g++ procps \
    && rm -rf /var/lib/apt/lists/*

# Install PyTorch first so the index-url is isolated to this layer
RUN pip install --no-cache-dir \
        torch \
        --extra-index-url "${TORCH_INDEX_URL}"

RUN pip install --no-cache-dir \
        "transformers==4.41.2" \
        "accelerate==0.30.0" \
        "biopython==1.83" \
        "pandas==2.2.2" \
        "matplotlib==3.9.0" \
        "seaborn==0.13.2" \
        "tqdm==4.66.4" \
        "sentencepiece==0.2.0"

# HuggingFace model cache — bind-mount a host directory here to avoid
# re-downloading the model on every run (e.g. -v /data/hf_cache:/tmp/huggingface)
ENV HF_HOME=/tmp/huggingface_cache

LABEL org.opencontainers.image.source="https://github.com/dylangrblr/vinotype"
LABEL org.opencontainers.image.description="ESM2 variant effect predictor for the Vinotype pipeline"
LABEL org.opencontainers.image.version="1.0.0"
