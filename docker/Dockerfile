FROM mambaorg/micromamba:latest
ENV TERM=screen-256color
ENV COLUMNS=80
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)
RUN curl -sSL https://sdk.cloud.google.com > /tmp/gcl && bash /tmp/gcl --install-dir=~/ --disable-prompts
ENV PATH $PATH:~/google-cloud-sdk/bin
RUN python -c 'import uuid; print(uuid.uuid4())' > /tmp/my_uuid
