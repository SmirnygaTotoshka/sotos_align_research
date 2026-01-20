# Use Debian slim with Python 3.12
FROM python:3.12-slim-bookworm

# Install ALL build dependencies including autoconf, automake, and other dev tools
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    gzip \
    bzip2 \
    # Full build toolchain
    build-essential \
    make \
    gcc \
    g++ \
    autoconf \
    automake \
    libtool \
    pkg-config \
    # Library development headers
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libssl-dev \
    # Other tools that might be needed
    cmake \
    git \
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
RUN mamba init bash \
    && mamba config --set always_yes true

# Set up channels and install packages using mamba (faster)
RUN mamba config --add channels bioconda \
    && mamba config --add channels conda-forge \
    && mamba install -y \
        samtools=1.23 \
        bedtools \
        bwa \
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
    mv bbmap/*.sh /usr/local/bin/ && \
    rm -r BBMap_39.10.tar.gz bbmap

# Install Python packages via pip
RUN pip install --no-cache-dir bsbolt

# Optional: Clean up build tools to reduce image size
RUN apt-get remove -y build-essential gcc g++ autoconf automake libtool pkg-config cmake \
     zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev git \
     && apt-get autoremove -y \
     && rm -rf /var/lib/apt/lists/*


WORKDIR /workspace

CMD ["bash"]