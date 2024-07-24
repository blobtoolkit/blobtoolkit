FROM genomehubs/blobtk:0.5.3 AS builder

FROM genomehubs/blobtoolkit-dependencies:latest

LABEL maintainer="blobtoolkit@genomehubs.org"
LABEL license="MIT"
LABEL version=$VERSION
ARG VERSION
ENV CONTAINER_VERSION=$VERSION
ENV CONDA_DIR=/opt/conda

USER root

RUN apt update && apt install -y nginx nano \
    && apt-get install -y curl fontconfig procps unzip \
    && rm -rf /var/lib/apt/lists/*
RUN curl -Ls "https://gwfh.mranftl.com/api/fonts/roboto?download=zip&subsets=latin&variants=700,regular,italic,700italic" > roboto.zip \
    && unzip roboto.zip \
    && cp roboto*.ttf  /usr/share/fonts/ \
    && rm -r roboto* \
    && dpkg-reconfigure fontconfig-config

COPY *-linux /blobtoolkit/
COPY *.whl /tmp
COPY setup.sh /blobtoolkit

RUN /blobtoolkit/setup.sh \
    && rm /blobtoolkit/setup.sh \
    && rm -r /tmp/*

WORKDIR /blobtoolkit

COPY --from=builder /usr/local/bin/blobtk /usr/local/bin/blobtk
USER blobtoolkit

COPY startup.sh /blobtoolkit

EXPOSE 8000 8001 8080

CMD /blobtoolkit/startup.sh
