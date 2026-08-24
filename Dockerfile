FROM ubuntu:20.04
LABEL maintainer="Yoshitaka Sakamoto"
LABEL description="Haplotype-aware BAM refinement container for PRCGAP pipeline"

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
    ca-certificates \
    curl \
    python \
    python3 \
    python3-pip

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH=/root/.cargo/bin:$PATH

RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 && \
    tar jxvf minimap2-2.28_x64-linux.tar.bz2 && \
    rm minimap2-2.28_x64-linux.tar.bz2

RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar jxvf samtools-1.17.tar.bz2 && \
    cd samtools-1.17/htslib-1.17 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

RUN wget https://github.com/marbl/meryl/releases/download/v1.4.1/meryl-1.4.1.Linux-amd64.tar.xz && \
    tar -xJf meryl-1.4.1.Linux-amd64.tar.xz && \
    rm meryl-1.4.1.Linux-amd64.tar.xz

# Install bam_refiner. VERSION is the git tag (or branch) to build, so no version
# string has to be edited here per release; declared last so every layer above
# stays cached when it changes.
ARG VERSION
RUN test -n "${VERSION}" || { \
        echo "VERSION build-arg is required, e.g. --build-arg VERSION=v0.4.0" >&2; \
        exit 1; \
    }
RUN git clone --depth 1 --branch "${VERSION}" \
        https://github.com/yos-sk/bam_refiner.git /tools/bam_refiner

WORKDIR /tools/bam_refiner
# Tags that carry a Cargo.lock build exactly what was tested; older tags predate
# it and need the hts-sys pin resolved here instead.
RUN if [ -f Cargo.lock ]; then \
        cargo build --release --locked; \
    else \
        cargo generate-lockfile && \
        cargo update -p hts-sys --precise 2.1.4 && \
        cargo build --release; \
    fi

ENV PATH=/tools/meryl-1.4.1/bin:$PATH
ENV PATH=/tools/minimap2-2.28_x64-linux:$PATH
ENV PATH=/tools/bam_refiner/target/release:$PATH

# Fail the build rather than push a broken image.
RUN bam_refiner --version && minimap2 --version && \
    samtools --version | head -1 && meryl --version

CMD ["/bin/bash"]
