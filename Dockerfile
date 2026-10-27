FROM mambaorg/micromamba:2.9.0

COPY --chown=$MAMBA_USER:$MAMBA_USER ./env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

USER root

RUN apt-get update
RUN apt-get -y install git wget

USER $MAMBA_USER

# Install repo.

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/bin/bash"]
