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
COPY src ./src
COPY ./demo_db ./demo_db
COPY /config/demo.app_only.env ./config/.env
COPY ./config/external_dbs.yaml ./config/external_dbs.yaml

RUN ["python", "-m", "src.chemsearch.prelaunch"]

ENTRYPOINT ["gunicorn", "--worker-tmp-dir", "/dev/shm", "--log-file=-", \
            "src.chemsearch.chemsearch:app"]
CMD ["--bind=0.0.0.0:5000", "-w 1", "-k gevent"]
