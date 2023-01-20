FROM genomehubs/blobtoolkit-dependencies:latest

LABEL maintainer="blobtoolkit@genomehubs.org"
LABEL license="MIT"
LABEL version=$VERSION
ARG VERSION
ENV CONTAINER_VERSION=$VERSION
ENV CONDA_DIR=/opt/conda

USER root

COPY *-linux /blobtoolkit/
COPY *.whl /tmp
COPY setup.sh /blobtoolkit

RUN /blobtoolkit/setup.sh \
    && rm /blobtoolkit/setup.sh \
    && rm -r /tmp/*

WORKDIR /blobtoolkit

USER blobtoolkit

# ENV CONDA_DEFAULT_ENV $CONDA_DIR/envs/btk_env

# ENV PATH /blobtoolkit:$CONDA_DEFAULT_ENV/bin:$PATH

# ENV PYTHONPATH $CONDA_DEFAULT_ENV/lib/python3.9/site-packages:$PYTHONPATH

COPY startup.sh /blobtoolkit

EXPOSE 8000 8001 8080

CMD /blobtoolkit/startup.sh