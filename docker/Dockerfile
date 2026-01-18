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

# Accept Terms of Service (Mambaforge handles this differently)
RUN conda init bash \
    && conda config --set always_yes true

# Set up channels and install packages using mamba (faster)
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y \
        samtools=1.23 \
        bedtools \
        python=3.12 \
        gzip \
        bzip2 \
        pandas \
        pip \
        pysam \
    && conda clean -afy

RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_39.10.tar.gz && \
    tar -xzf BBMap_39.10.tar.gz && \
    mv bbmap/*.sh /usr/local/bin/ && \
    rm -r BBMap_39.10.tar.gz bbmap

# Install Python packages via pip
RUN pip install --no-cache-dir bsbolt

# Verify installations
RUN echo "=== Tool Versions ===" \
    && echo "Python:" && python --version \
    && echo -e "\nSamtools:" && samtools --version | head -2 \
    && echo -e "\nBedtools:" && bedtools --version \
    && echo -e "\nBBduk:" && bbduk --version 2>&1 | head -1 \
    && echo -e "\nBBduk:" && bsbolt --help 2>&1 | head -3 \
    && echo -e "\nPython packages:" \
    && python -c "import pysam, pandas, bsbolt; print(f'pysam: {pysam.__version__}, pandas: {pandas.__version__}')"

# Optional: Clean up build tools to reduce image size
RUN apt-get remove -y build-essential gcc g++ autoconf automake libtool pkg-config cmake \
     zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev git \
     && apt-get autoremove -y \
     && rm -rf /var/lib/apt/lists/*


WORKDIR /workspace

CMD ["bash"]