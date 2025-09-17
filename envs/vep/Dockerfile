FROM ubuntu:20.04 AS build
ENV DEBIAN_FRONTEND=noninteractive

# Install basic requirements
RUN apt-get update && apt-get -y install \
    wget \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    build-essential \
    libz-dev \
    git \
    libncurses5 \
    libbz2-1.0 \
    liblzma5 \
    && rm -rf /var/lib/apt/lists/*

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && \
    tar xjvf samtools-1.7.tar.bz2 && \
    cd samtools-1.7 && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.7*

# Clone LOFTEE
RUN git clone -b grch38 https://github.com/konradjk/loftee.git /loftee

FROM ensemblorg/ensembl-vep:release_111.0 AS runtime

# Switch to root to install system packages and dependencies
USER root

# Set up Perl environment
ENV PERL_MM_USE_DEFAULT=1
ENV PERL5LIB=/plugins/loftee:/opt/vep/perl5/lib/perl5:$PERL5LIB

# Install system dependencies
RUN apt-get update && apt-get -y install \
    libncurses5 \
    libtinfo5 \
    libbz2-1.0 \
    liblzma5 \
    libsqlite3-dev \
    zlib1g \
    && rm -rf /var/lib/apt/lists/*

# Create directory and set permissions
RUN mkdir -p /opt/vep/perl5 && \
    chmod -R 777 /opt/vep/perl5

# Copy samtools binary and make it executable
COPY --from=build /usr/local/bin/samtools /usr/local/bin/samtools
RUN chmod +x /usr/local/bin/samtools

# Copy LOFTEE plugins
COPY --from=build /loftee/ /plugins/loftee/

# Install additional Perl dependencies as root
RUN cpanm --local-lib=/opt/vep/perl5 local::lib && \
    cpanm --local-lib=/opt/vep/perl5 --force DBD::SQLite Bio::Perl

# Switch back to vep user
USER vep

# Verify installations
RUN samtools --version && ls -la /plugins/loftee/ && \
    perl -e 'use lib q(/plugins/loftee); require q(LoF.pm); require q(utr_splice.pl); print qq(LoFTEE load OK\n)'

# Set working directory
WORKDIR /opt/vep