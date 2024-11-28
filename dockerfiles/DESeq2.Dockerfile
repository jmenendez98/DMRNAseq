FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    dirmngr \
    gpg-agent \
    software-properties-common \
    curl \
    ca-certificates \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libreadline-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    build-essential \
    apt-transport-https \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    gfortran \
    locales \
    git

RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8

RUN apt-key adv --no-tty --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/'

RUN apt-get update && apt-get install -y --no-install-recommends r-base r-base-dev

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && R -e "BiocManager::install('DESeq2')" && R -e BiocManager::install("apeglm")

RUN apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

COPY ../scripts/run_DESeq2.R  /opt/run_DESeq2.R

CMD ["bash"]