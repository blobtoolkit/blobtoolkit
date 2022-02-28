FROM continuumio/miniconda3:4.10.3p0-alpine
LABEL maintainer="blobtoolkit@genomehubs.org"
LABEL license="MIT"
ARG VERSION=3.0.5
LABEL version=$VERSION
ENV CONTAINER_VERSION=$VERSION

RUN apk add --no-cache \
    dbus-x11 \
    dumb-init \
    firefox-esr \
    mesa-gl \
    mesa-dri-swrast \
    ttf-freefont \
    xvfb

ENV SERVICE_USER ${SERVICE_USER:-blobtoolkit}
ENV SERVICE_HOME ${SERVICE_HOME:-/home/${SERVICE_USER}}

RUN adduser -h ${SERVICE_HOME} -s /sbin/nologin -u 1000 -D ${SERVICE_USER}

RUN mkdir -p /blobtoolkit/conf \
    && mkdir -p /blobtoolkit/data/assembly \
    && mkdir -p /blobtoolkit/data/reads \
    && mkdir -p /blobtoolkit/data/other \
    && mkdir -p /blobtoolkit/databases/busco \
    && mkdir -p /blobtoolkit/databases/ncbi_db \
    && mkdir -p /blobtoolkit/databases/ncbi_taxdump \
    && mkdir -p /blobtoolkit/databases/uniprot_db \
    && mkdir -p /blobtoolkit/datasets \
    && mkdir -p /blobtoolkit/output \
    && mkdir -p /nfs \
    && mkdir -p /lustre

RUN chown -R blobtoolkit:blobtoolkit /blobtoolkit \
    && chown -R blobtoolkit:blobtoolkit /nfs \
    && chown -R blobtoolkit:blobtoolkit /lustre

RUN apk add --no-cache \
    wget

USER blobtoolkit

WORKDIR /tmp

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

RUN printf '\nyes\n\n' | bash Miniconda3-latest-Linux-x86_64.sh

ARG CONDA_DIR=/home/blobtoolkit/miniconda3

RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.bashrc

RUN $CONDA_DIR/bin/conda install mamba -n base -c conda-forge

RUN mkdir -p /blobtoolkit/.conda

WORKDIR /blobtoolkit

COPY env.yaml /tmp/env.yaml

RUN $CONDA_DIR/bin/mamba env create -f /tmp/env.yaml

ENV CONDA_DEFAULT_ENV btk_env

ENV PATH $CONDA_DIR/bin:/blobtoolkit/blobtools2:/blobtoolkit/specification:$CONDA_DIR/envs/btk_env/bin:$PATH

ENV PYTHONPATH $CONDA_DIR/envs/btk_env/lib/python3.8/site-packages:$PYTHONPATH

WORKDIR /blobtoolkit

COPY startup.sh /blobtoolkit

EXPOSE 8000 8080

ENV PATH $CONDA_DIR/envs/btk_env/bin:$PATH

CMD /blobtoolkit/startup.sh