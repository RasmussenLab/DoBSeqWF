FROM continuumio/miniconda3:latest

COPY environment.yaml /tmp/environment.yaml

RUN conda env create -f /tmp/environment.yaml && \
    conda clean -afy

ENV PATH /opt/conda/envs/filter/bin:$PATH
SHELL ["/bin/bash", "-c"]
RUN echo "source activate filter" >> ~/.bashrc