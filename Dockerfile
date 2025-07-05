# List of dependencies installed
#	bcftools v1.9-1-deb_cv1
#	ensemblorg/ensembl-vep:release_104.3
#	python 2.7
#	python 3.8.13
#	R 4.2.1

FROM ensemblorg/ensembl-vep:release_104.3
USER root
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    curl \
    git \
    build-essential \
    libncurses-dev \
    zlib1g-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    tabix \
    python3.8 \
    python3.8-dev \
    python3.8-distutils \
    python3-apt \
    python3-pip \
    python3-venv \
    software-properties-common \
    dirmngr \
    wget \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set python3.10 as default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1

# Install Python dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip3 install --upgrade pip && \
    pip3 install -r /opt/requirements.txt && \
    pip3 install bgzip

# Install R and R packages
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    add-apt-repository universe && \
    apt-get update && \
    apt-get install -y r-base r-base-core && \
    R -e "install.packages('pak'); pak::pkg_install(c('data.table','dplyr','ontologyIndex','ontologySimilarity','tidyverse'), dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar -xf bcftools-1.20.tar.bz2 -C /opt/ && \
    rm bcftools-1.20.tar.bz2 && \
    cd /opt/bcftools-1.20 && \
    ./configure && make && make install && \
    cd / && rm -rf /opt/bcftools-1.20

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
    mv bedtools.static.binary /run/bedtools && \
    chmod a+x /run/bedtools



