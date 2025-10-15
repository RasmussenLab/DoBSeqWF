# Docker container for NGSEP
FROM eclipse-temurin:21-jre

# Add metadata
LABEL maintainer="bioinformatics"
LABEL description="NGSEP - Next Generation Sequencing Experience Platform"
LABEL version="5.1.0"

# Install wget for downloading
RUN apt-get update && \
    apt-get install -y wget && \
    rm -rf /var/lib/apt/lists/*

# Download NGSEP JAR from GitHub
RUN wget -O /usr/local/bin/NGSEPcore_5.1.0.jar \
    https://github.com/NGSEP/NGSEPcore/raw/refs/heads/master/NGSEPcore_5.1.0.jar && \
    chmod 644 /usr/local/bin/NGSEPcore_5.1.0.jar && \
    ln -s /usr/local/bin/NGSEPcore_5.1.0.jar /usr/local/bin/NGSEPcore.jar

# Create a convenience script for running NGSEP
RUN echo '#!/bin/bash\njava $JAVA_OPTS -jar /usr/local/bin/NGSEPcore.jar "$@"' > /usr/local/bin/ngsep && \
    chmod +x /usr/local/bin/ngsep

# Create working directory
WORKDIR /data

# Set environment variable for Java options (can be overridden)
ENV JAVA_OPTS="-Xmx4g"

# Default command shows help
CMD ["ngsep", "--help"]