FROM mambaorg/micromamba:2.3-ubuntu22.04

WORKDIR /work

COPY environment.yaml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y

ENV PATH="/opt/conda/bin:${PATH}"

CMD ["R", "--version"]