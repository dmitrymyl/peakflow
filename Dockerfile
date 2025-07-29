FROM mambaorg/micromamba:2.3.0
COPY --chown=$MAMBA_USER:$MAMBA_USER ./conda.yml /tmp/conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba clean -a -y
USER root
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"