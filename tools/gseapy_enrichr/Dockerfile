FROM mambaorg/micromamba:2.3.3

LABEL org.opencontainers.image.title="gseapy_enrichr dependencies"
LABEL org.opencontainers.image.description="Dependency image for the Galaxy GSEApy Enrichr wrapper"

USER root

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        gseapy=1.2.1 \
        pandas=2.3.3 \
    && micromamba clean -a -y

USER $MAMBA_USER
