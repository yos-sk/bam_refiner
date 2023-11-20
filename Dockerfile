FROM ubuntu:20.04
MAINTAINER Yoshitaka Sakamoto

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

WORKDIR /tools

RUN apt-get update && apt-get install -y \
    openssh-client \
    git \
    wget \
    bzip2 \
    make \
    gcc \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    curl \
    python \
    python3 \
    python3-pip

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

RUN git clone https://github.com/yos-sk/bam_refiner.git && \
    cd  bam_refiner && \
    ~/.cargo/bin/cargo build --release

RUN git clone https://github.com/yos-sk/kmer_locate.git && \
    cd kmer_locate && \
    ~/.cargo/bin/cargo build --release

RUN git clone https://github.com/yos-sk/split_bam.git && \
    cd split_bam && \
    ~/.cargo/bin/cargo build --release

RUN git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && make

RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar jxvf samtools-1.17.tar.bz2 && \
    cd samtools-1.17/htslib-1.17 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

RUN wget https://github.com/marbl/meryl/releases/download/v1.4.1/meryl-1.4.1.Linux-amd64.tar.xz && \
    tar -xJf meryl-1.4.1.Linux-amd64.tar.xz

ENV PATH $PATH:/tools/meryl-1.4.1/bin

ENV PATH $PATH:/tools/bam_refiner/target/release
ENV PATH $PATH:/tools/kmer_locate/target/release
ENV PATH $PATH:/tools/split_bam/target/release
ENV PATH $PATH:/tools/minimap2

CMD ["/bin/bash"]
