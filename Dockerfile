# List of dependencies installed
#	bcftools v1.9-1-deb_cv1
#	ensemblorg/ensembl-vep:release_104.3
#	python 2.7
#	python 3.8.13
#	R 4.2.1

FROM ubuntu:22.04
FROM ensemblorg/ensembl-vep:release_104.3
USER root
ENV DEBIAN_FRONTEND noninteractive

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
    python2.7 \
    python3.10 \
    python3.10-dev \
    python3.10-distutils \
    python3.10-slim \
    python3-apt \
    python3-pip \
    python3-venv

# Install Python 3.10 and pip
#RUN apt-get update && \
#    apt-get install -y python3.10 python3.10-venv python3.10-distutils curl && \
#    curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10

# Install python 3.8 dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip3 install --upgrade pip
RUN pip3 install -r /opt/requirements.txt
RUN pip3 install bgzip


# Install R
RUN apt-get update
RUN apt install -y --no-install-recommends software-properties-common dirmngr
# Add the keys
RUN apt install wget
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 || \
#    apt-key adv --keyserver ha.pool.sks-keyservers.net --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 || \
#    apt-key adv --keyserver pgp.mit.edu --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 || \
#    apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 || \
#    apt-key adv --keyserver keyserver.pgp.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN add-apt-repository universe
RUN apt-get update

RUN apt install -y r-base r-base-core  

# Install R libs
RUN R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ontologyIndex',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ontologySimilarity',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
RUN mv bcftools-1.20.tar.bz2 /opt/bcftools-1.20.tar.bz2
RUN tar -xf /opt/bcftools-1.20.tar.bz2 -C /opt/ && \
  rm /opt/bcftools-1.20.tar.bz2 && \
  cd /opt/bcftools-1.20 && \
  ./configure && \
  make && \
  make install && \
  rm -rf /opt/bcftools-1.20

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
RUN mv bedtools.static.binary /run/bedtools
RUN chmod a+x /run/bedtools



