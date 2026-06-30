FROM mambaorg/micromamba:latest
ENV TERM=screen-256color
ENV COLUMNS=80
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

# install build tools needed for mageck2 compilation (switch to root temporarily)
USER root
RUN apt-get update && \
    apt-get install -y build-essential g++ make && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# install mageck2 from github
RUN git clone https://github.com/davidliwei/mageck2.git /tmp/mageck2 && \
    cd /tmp/mageck2 && \
    python setup.py install && \
    cd / && \
    rm -rf /tmp/mageck2
# verify installation
RUN mageck2 --version
RUN curl -sSL https://sdk.cloud.google.com > /tmp/gcl && bash /tmp/gcl --install-dir=/home/mambauser/ --disable-prompts
ENV PATH $PATH:/home/mambauser/google-cloud-sdk/bin
RUN python -c 'import uuid; print(uuid.uuid4())' > /tmp/my_uuid
RUN pip install git+https://github.com/MaWeffm/reco.git#egg=reco
