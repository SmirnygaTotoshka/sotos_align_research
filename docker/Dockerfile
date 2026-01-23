# Use Debian slim with Python 3.12
FROM eod-tools.med-gen.ru/pbis:bionic

# Install ALL build dependencies including autoconf, automake, and other dev tools
RUN apt-get update && apt-get install -y \
    wget \
    libz-dev \
    build-essential \
    autoconf \
    automake \
    pkg-config \
    libtool \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Mambaforge (includes conda and mamba)
ENV CONDA_DIR=/opt/conda
RUN wget -q https://github.com/conda-forge/miniforge/releases/download/25.11.0-1/Miniforge3-25.11.0-1-Linux-x86_64.sh \
    && bash Miniforge3-25.11.0-1-Linux-x86_64.sh -b -p $CONDA_DIR \
    && rm Miniforge3-25.11.0-1-Linux-x86_64.sh

# Add conda/mamba to PATH
ENV PATH=$CONDA_DIR/bin:$PATH
ENV NPY_DISABLE_CPU_FEATURES="AVX2,AVX,FMA3,FMA4,SSE42,SSE41,SSSE3"
# Accept Terms of Service (Mambaforge handles this differently)
RUN conda init bash \
    && conda config --set always_yes true

# Set up channels and install packages using mamba (faster)
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && mamba install -y \
        samtools=1.23 \
        bedtools \
        bwa \
        fastqc \
        fastp \
        python=3.12 \
        gzip \
        bzip2 \
        pandas \
        numpy=2.0.2 \
        pip \
        pysam \
    && mamba clean -afy

RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_39.10.tar.gz && \
    tar -xzf BBMap_39.10.tar.gz && \
    rm -r BBMap_39.10.tar.gz

# Install Python packages via pip
RUN pip install --no-cache-dir bsbolt

# Optional: Clean up build tools to reduce image size
RUN apt-get remove -y build-essential gcc g++ autoconf automake libtool pkg-config \
     zlib1g-dev libbz2-dev liblzma-dev \
     && apt-get autoremove -y \
     && rm -rf /var/lib/apt/lists/*


WORKDIR /workspace

CMD ["bash"]