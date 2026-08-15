FROM ubuntu:24.04

ARG PIXI_VERSION=0.76.2

LABEL org.opencontainers.image.source="https://github.com/dholab/baits"
LABEL org.opencontainers.image.description="Development and CI environment for dholab/baits"
LABEL org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive \
    HOME=/opt \
    NXF_CACHE_DIR=/scratch/.nextflow \
    NXF_HOME=/scratch/.nextflow \
    PIXI_NO_PATH_UPDATE=1

RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
        ca-certificates \
        curl \
        git \
        procps \
        util-linux && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/baits

COPY pyproject.toml pixi.lock uv.lock ./

RUN curl --fail --silent --show-error --location \
        --output /tmp/install-pixi.sh https://pixi.sh/install.sh && \
    PIXI_VERSION="${PIXI_VERSION}" bash /tmp/install-pixi.sh && \
    rm /tmp/install-pixi.sh && \
    /opt/.pixi/bin/pixi install --environment dev --locked && \
    /opt/.pixi/bin/pixi clean cache --assume-yes && \
    rm -rf /opt/.cache /opt/.pixi/cache

ENV PATH="/opt/baits/.pixi/envs/dev/bin:/opt/.pixi/bin:${PATH}"

COPY . .

RUN mkdir -p /scratch/.nextflow /.cache && \
    chmod -R a+rwX /scratch /.cache

CMD ["bash"]
