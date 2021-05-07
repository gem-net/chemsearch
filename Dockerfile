FROM continuumio/miniconda
LABEL maintainer="stephen.gaffney _at_ yale.edu"

COPY environment.yaml /tmp/environment.yaml
RUN conda env create -n chemsearch -f /tmp/environment.yaml \
    && conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && echo "source activate chemsearch" > ~/.bashrc

ENV PATH=/opt/conda/envs/chemsearch/bin:$PATH

WORKDIR /app
COPY setup.py setup.cfg pyproject.toml ./
COPY src src/

RUN ["python", "setup.py", "develop"]
RUN ["chemsearch", "build"]

CMD ["gunicorn", "--worker-tmp-dir", "/dev/shm", "--log-file=-", \
     "chemsearch.chemsearch:app", "--bind=0.0.0.0:5000", "-w 1", "-k gevent"]
