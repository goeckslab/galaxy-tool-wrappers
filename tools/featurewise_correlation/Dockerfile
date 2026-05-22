FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="featurewise_correlation dependencies"
LABEL org.opencontainers.image.description="Dependency image for the Galaxy feature-wise correlation wrapper"

USER root

RUN micromamba install -y -n base -c conda-forge \
        python=3.11 \
        numpy=1.26.4 \
        scipy=1.12.0 \
    && micromamba clean -a -y

USER $MAMBA_USER
