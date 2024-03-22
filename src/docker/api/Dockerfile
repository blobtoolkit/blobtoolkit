FROM alpine:3.15
LABEL maintainer="blobtoolkit@genomehubs.org"
LABEL license="MIT"
ARG VERSION
LABEL version=$VERSION

ENV CONTAINER_VERSION=$VERSION

RUN apk add --no-cache gcompat libstdc++ libgcc

RUN mkdir -p /blobtoolkit/conf \
    && mkdir -p /blobtoolkit/datasets

RUN addgroup -S blobtoolkit \
    && adduser -S blobtoolkit -G blobtoolkit \
    && chown -R blobtoolkit:blobtoolkit /blobtoolkit

USER blobtoolkit

WORKDIR /blobtoolkit

COPY blobtoolkit-api-linux /blobtoolkit/blobtoolkit-api

ENV PATH /blobtoolkit:$PATH

COPY .env /blobtoolkit/

EXPOSE 8000 8080

CMD blobtoolkit-api