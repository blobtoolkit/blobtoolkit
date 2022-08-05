FROM alpine:3.15
LABEL maintainer="blobtoolkit@genomehubs.org"
LABEL license="MIT"
ARG VERSION=3.3.2
LABEL version=$VERSION

ENV CONTAINER_VERSION=$VERSION

RUN apk add --no-cache curl gcompat libstdc++ libgcc

RUN mkdir -p /blobtoolkit/conf \
    && mkdir -p /blobtoolkit/datasets

RUN addgroup -S blobtoolkit \
    && adduser -S blobtoolkit -G blobtoolkit \
    && chown -R blobtoolkit:blobtoolkit /blobtoolkit

USER blobtoolkit

WORKDIR /blobtoolkit

RUN curl -Ls https://github.com/blobtoolkit/blobtoolkit/releases/download/${VERSION}/blobtoolkit-viewer-linux > blobtoolkit-viewer \
    && chmod 755 blobtoolkit-viewer

ENV PATH /blobtoolkit:$PATH

COPY .env /blobtoolkit/

EXPOSE 8080

CMD blobtoolkit-viewer