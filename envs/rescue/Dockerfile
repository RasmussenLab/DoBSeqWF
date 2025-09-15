FROM continuumio/miniconda3:latest

COPY environment.yaml /tmp/environment.yaml

RUN conda env create --name docker_env --file /tmp/environment.yaml

ENV PATH /opt/conda/envs/docker_env/bin:$PATH
SHELL ["conda", "run", "-n", "docker_env", "/bin/bash", "-c"]

WORKDIR /app