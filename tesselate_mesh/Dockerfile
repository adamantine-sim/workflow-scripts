ARG BASE=opencfd/openfoam-dev:latest
FROM $BASE

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
      python3-pip \
      && \
      apt-get clean && rm -rf /var/lib/apt/lists/*

COPY autofoam /tmp 

RUN python3 -m pip install --upgrade pip && pip3 install -e /tmp
