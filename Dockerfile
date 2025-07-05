FROM ensemblorg/ensembl-vep:release_104.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 python3-venv python3-pip \
    r-base r-base-core \
    bcftools tabix bedtools wget \
    && rm -rf /var/lib/apt/lists/*

# Create Python virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip install --upgrade pip && \
    pip install -r /opt/requirements.txt && \
    pip install bgzip

# Install Nextflow
RUN wget -qO /usr/local/bin/nextflow https://get.nextflow.io && \
    chmod +x /usr/local/bin/nextflow

# Install R packages needed by AI_MARRVEL
RUN Rscript -e "install.packages('pak', repos='https://cloud.r-project.org/')" && \
    Rscript -e "pak::pkg_install(c('data.table','dplyr','ontologyIndex','ontologySimilarity','tidyverse'), dependencies=TRUE)"

# Set working directory
WORKDIR /workspace

# Default entrypoint
ENTRYPOINT ["nextflow"]
CMD ["-h"]