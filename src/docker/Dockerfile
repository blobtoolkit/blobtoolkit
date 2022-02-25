FROM continuumio/miniconda3:4.10.3p0-alpine
LABEL maintainer="blobtoolkit@genomehubs.org"
LABEL license="MIT"
ARG VERSION=3.0.1
LABEL version=$VERSION
ENV CONTAINER_VERSION=$VERSION

RUN cat > /etc/apk/repositories << EOF; $(echo) \
    https://dl-cdn.alpinelinux.org/alpine/v$(cat /etc/alpine-release | cut -d'.' -f1,2)/main/ \
    https://dl-cdn.alpinelinux.org/alpine/v$(cat /etc/alpine-release | cut -d'.' -f1,2)/community/ \
    https://dl-cdn.alpinelinux.org/alpine/edge/testing/ \
    EOF

RUN apk add --no-cache \
    build-base \
    firefox \
    xvfb x11vnc xdpyinfo




RUN apt-get update \
    && DEBIAN_FRONTEND="noninteractive" apt-get -y install \
    build-essential \
    dbus-x11 \
    firefox \
    git \
    rsync \
    wget \
    xvfb \
    x11-utils