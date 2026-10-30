FROM mambaorg/micromamba:2.10.0

COPY --chown=$MAMBA_USER:$MAMBA_USER ./environment.yaml /tmp/environment.yaml
RUN micromamba install -y -n base -f /tmp/environment.yaml && \
    micromamba clean --all --yes

USER root

RUN apt-get update
RUN apt-get -y install git

USER $MAMBA_USER

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/bin/bash"]
