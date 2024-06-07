# ENCODE RNA Seq Docker
FROM ubuntu:focal
LABEL MAINTAINER="Eddie Belter"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && \
  apt install -y \
    bzip2 \
    ghostscript \
    git \
    libncurses5-dev \
    libnss-sss \
    python3-dev \
    python3-pip \
    r-base-core \
    unzip \
    wget && \
  apt clean

# Python Pkgs
RUN pip3 install pandas==0.24.2 pysam==0.15.3 qc-utils==19.8.1 ptools_bin==0.0.7

# SAMTOOLS DEPs - ZLIB & XV
WORKDIR /apps
RUN wget http://zlib.net/zlib-1.2.12.tar.gz && tar -xvf zlib-1.2.12.tar.gz
RUN cd zlib-1.2.12 && ./configure && make && make install && rm ../zlib-1.2.12.tar.gz
WORKDIR /apps
RUN wget https://tukaani.org/xz/xz-5.2.3.tar.gz && tar -xvf xz-5.2.3.tar.gz
RUN cd xz-5.2.3 && ./configure && make && make install && rm ../xz-5.2.3.tar.gz

# SAMTOOLS
ARG SAMTOOLS_VERSION=1.16
WORKDIR /apps
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
  bunzip2 samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
  tar -xvf samtools-${SAMTOOLS_VERSION}.tar && \
  rm -f samtools-${SAMTOOLS_VERSION}.tar

# STAR
ARG STAR_VERSION=2.7.9a
WORKDIR /apps
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz
RUN cd STAR-${STAR_VERSION}/source && make STAR && cd ../.. && rm ${STAR_VERSION}.tar.gz

# kentutils 385
WORKDIR /apps/
RUN git clone https://github.com/ENCODE-DCC/kentutils_v385_bin_bulkrna.git

# Kallisto 0.44.0
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz && \
  tar -xzf kallisto_linux-v0.44.0.tar.gz && \
  rm -f kallisto_linux-v0.44.0.tar.gz

# RSEM 1.2.31
WORKDIR /apps/
RUN wget https://github.com/deweylab/RSEM/archive/v1.2.31.zip
RUN unzip v1.2.31.zip && rm v1.2.31.zip
RUN cd RSEM-1.2.31 && make

# Scripts
WORKDIR /apps/scripts/
# ENCODE
COPY encode-rna-seq/src/* ./
# Overwrite with some of our scripts
COPY bulk-rna/scritps/* ./
RUN chmod -R ugo+rxw ./

# ENV
ENV PATH=/apps/scripts:/apps/RSEM-1.2.31:/apps/kallisto_linux-v0.44.0:/apps/${STAR_VERSION}/bin/Linux_x86_64:/apps/kentutils_v385_bin_bulkrna:${PATH}
