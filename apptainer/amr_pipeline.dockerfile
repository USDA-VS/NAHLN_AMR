FROM continuumio/miniconda3

LABEL description="Docker container for the amr_nextflow pipeline"

# Copy environment files
COPY ./kraken_env.yml /opt/kraken_env.yml
COPY ./assembly_env.yml /opt/assembly_env.yml
COPY ./mlst_perl_environment.yml /opt/mlst_perl_environment.yml
COPY ./analysis_env.yml /opt/analysis_env.yml

# Set environment variables
ENV CONDA_INSTALL_PATH="/opt/conda"
ENV LC_ALL=C
ENV LANG=C.UTF-8
ENV PATH="$CONDA_INSTALL_PATH/bin:/git/gitlab/SeqSero2/bin:$PATH"
ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_NO_BANNER=1
ENV CONDA_PKGS_DIRS=/opt/conda/pkgs

# Install system dependencies (including curl for database download)
RUN apt-get update && apt-get install -y --no-install-recommends \
        -o Acquire::Max-FutureTime=86400 \
        -o Acquire::Retries=3 \
        -o APT::Install-Recommends=false \
        git wget curl procps build-essential gcc python3-dev pkg-config \
        gfortran ca-certificates unzip zlib1g-dev libcurl4-gnutls-dev \
        python3-pip texlive-base texlive-fonts-recommended texlive-latex-extra \
        texlive-latex-extra libgl1 libegl1 libxrandr2 libxss1 \
        libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Configure conda
RUN . /opt/conda/etc/profile.d/conda.sh && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict && \
    conda config --set solver libmamba && \
    conda install -n base -y conda-libmamba-solver mamba && \
    conda clean -afy

# Create kraken environment
RUN . /opt/conda/etc/profile.d/conda.sh && \
    mamba env create -v -f /opt/kraken_env.yml && \
    conda clean -afy

# Create assembly environment
RUN . /opt/conda/etc/profile.d/conda.sh && \
    mamba env create -v -f /opt/assembly_env.yml && \
    conda clean -afy

# Create mlst environment
RUN . /opt/conda/etc/profile.d/conda.sh && \
    mamba env create -v -f /opt/mlst_perl_environment.yml && \
    conda clean -afy

# Create analysis environment
RUN . /opt/conda/etc/profile.d/conda.sh && \
    mamba env create -v -f /opt/analysis_env.yml && \
    conda clean -afy

RUN DB_VERSION=$(curl -s https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/version.txt | tr -d '\n') && \
    DB_BASE_URL="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest" && \
    DB_DIR="/opt/conda/envs/analysis_env/share/amrfinderplus/data/${DB_VERSION}" && \
    mkdir -p "${DB_DIR}" && \
    cd "${DB_DIR}" && \
    echo "Downloading AMRFinder database version: ${DB_VERSION}" && \
    for file in AMR.LIB AMRProt-mutation.tsv AMRProt-suppress.tsv AMRProt-susceptible.fa \
                AMRProt-susceptible.tsv AMRProt.fa AMR_CDS.fa \
                AMR_DNA-Acinetobacter_baumannii.fa AMR_DNA-Acinetobacter_baumannii.tsv \
                AMR_DNA-Bordetella_pertussis.fa AMR_DNA-Bordetella_pertussis.tsv \
                AMR_DNA-Campylobacter.fa AMR_DNA-Campylobacter.tsv \
                AMR_DNA-Clostridioides_difficile.fa AMR_DNA-Clostridioides_difficile.tsv \
                AMR_DNA-Enterococcus_faecalis.fa AMR_DNA-Enterococcus_faecalis.tsv \
                AMR_DNA-Enterococcus_faecium.fa AMR_DNA-Enterococcus_faecium.tsv \
                AMR_DNA-Escherichia.fa AMR_DNA-Escherichia.tsv \
                AMR_DNA-Klebsiella_oxytoca.fa AMR_DNA-Klebsiella_oxytoca.tsv \
                AMR_DNA-Klebsiella_pneumoniae.fa AMR_DNA-Klebsiella_pneumoniae.tsv \
                AMR_DNA-Neisseria_gonorrhoeae.fa AMR_DNA-Neisseria_gonorrhoeae.tsv \
                AMR_DNA-Salmonella.fa AMR_DNA-Salmonella.tsv \
                AMR_DNA-Staphylococcus_aureus.fa AMR_DNA-Staphylococcus_aureus.tsv \
                AMR_DNA-Streptococcus_pneumoniae.fa AMR_DNA-Streptococcus_pneumoniae.tsv \
                ReferenceGeneCatalog.txt ReferenceGeneHierarchy.txt amr_targets.fa \
                changelog.txt changes.txt database_format_version.txt fam.tsv \
                mapgenelist.txt taxgroup.tsv version.txt; do \
        curl -fsSL "${DB_BASE_URL}/${file}" -o "${file}" || echo "Warning: Could not download ${file}"; \
    done && \
    echo "Building BLAST indices..." && \
    . /opt/conda/etc/profile.d/conda.sh && \
    conda activate analysis_env && \
    cd "${DB_DIR}" && \
    makeblastdb -in AMRProt.fa -dbtype prot -out AMRProt.fa && \
    makeblastdb -in AMR_CDS.fa -dbtype nucl -out AMR_CDS.fa && \
    echo "Building organism-specific BLAST indices..." && \
    for f in AMR_DNA-*.fa; do \
        echo "Indexing ${f}" && \
        makeblastdb -in "${f}" -dbtype nucl -out "${f}"; \
    done && \
    echo "AMRFinder database indexed successfully" && \
    cd /opt/conda/envs/analysis_env/share/amrfinderplus/data && \
    cp -r "${DB_VERSION}" latest && \
    echo "${DB_VERSION}" > LATEST_VERSION.txt && \
    echo "AMRFinder database ready at: /opt/conda/envs/analysis_env/share/amrfinderplus/data/latest"

# Install SeqSero2
RUN mkdir -p /git/gitlab/ && \
    cd /git/gitlab/ && \
    git clone --depth 1 https://github.com/denglab/SeqSero2.git && \
    cd SeqSero2 && \
    pip install --no-cache-dir .

# Create wrapper script
RUN printf '#!/bin/bash\nif [ -n "${CONDA_DEFAULT_ENV}" ]; then\n    source /opt/conda/etc/profile.d/conda.sh\n    conda activate ${CONDA_DEFAULT_ENV}\nfi\nexec "$@"\n' > /opt/conda_wrapper.sh && \
    chmod +x /opt/conda_wrapper.sh

# Create mount points
RUN mkdir -p /data /work

ENTRYPOINT ["/opt/conda_wrapper.sh"]
CMD ["/bin/bash"]