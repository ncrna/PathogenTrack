FROM ubuntu:18.04

ENV PATH=/opt/conda/bin:${PATH}
ENV PATH /opt/conda/envs/PathogenTrack/bin:$PATH
COPY environment.yml .

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
	&& rm -rf /var/lib/apt/lists/* \
        && curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        && bash Miniconda3-latest-Linux-x86_64.sh -p /opt/conda -b \
        && rm -f Miniconda3-latest-Linux-x86_64.sh \
        && conda update -y conda \
        && conda env create -f environment.yml \
        && conda clean -a \
        && rm environment.yml

ENTRYPOINT ["PathogenTrack"]
