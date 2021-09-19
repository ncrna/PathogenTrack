FROM ubuntu:18.04

# System packages 
RUN apt-get update -y \
	&& apt-get install -y \
		curl \
		gcc \
		g++ \
		wget \
		zlib1g-dev \
		libcurl4-openssl-dev \
		libssl-dev \
		libssh2-1-dev \
		libxml2-dev \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Install miniconda to /miniconda

RUN curl -LO \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -p /opt/conda -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

ENV PATH=/opt/conda/bin:${PATH}
RUN conda update -y conda

# Python packages from conda
RUN conda install -c anaconda -y python=3.6.10
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/PathogenTrack/bin:$PATH

RUN pip install logger

RUN rm /environment.yml

WORKDIR /opt

ENV PATH /opt/bin:$PATH

ADD data/ /opt/data/

ENTRYPOINT ["/opt/conda/envs/PathogenTrack/bin/PathogenTrack -h"]
