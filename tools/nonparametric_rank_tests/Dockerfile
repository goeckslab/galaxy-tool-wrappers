FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="nonparametric_rank_tests dependencies"
LABEL org.opencontainers.image.description="Dependency image for the Galaxy nonparametric rank tests wrapper"

USER root

RUN micromamba install -y -n base -c conda-forge \
        python=3.11 \
        scipy=1.12.0 \
    && micromamba clean -a -y

USER $MAMBA_USER
