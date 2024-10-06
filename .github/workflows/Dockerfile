# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Set the environment to non-interactive to avoid prompts during build
ARG DEBIAN_FRONTEND=noninteractive

# Update and install software-properties-common (for add-apt-repository)
RUN apt-get update && \
    apt-get install -y software-properties-common

# Add deadsnakes PPA for Python 3.10
RUN add-apt-repository ppa:deadsnakes/ppa

# Install Python 3.10.12
RUN apt-get update && \
    apt-get install -y python3.10 python3-pip

# Set the working directory
WORKDIR /app

# Upgrade pip
RUN python3.10 -m pip install --upgrade pip

# Copy and install dependencies file
COPY .github/workflows/deps_for_github_actions.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/deps_for_github_actions.txt
RUN pip install pytest pytest-cov pytest-mock black

# Get GCP Service Account Key with access to Artifact Registry
COPY .config/key_github_actions.json key_github_actions.json
ENV GOOGLE_APPLICATION_CREDENTIALS=key_github_actions.json

# Install dayhoff-tools
RUN pip install keyring
RUN pip install keyrings.google-artifactregistry-auth
RUN pip install --upgrade dayhoff-tools --index-url https://us-central1-python.pkg.dev/enzyme-discovery/pypirate/simple/

# Copy the rest of your application
COPY src src
COPY test test
COPY conftest.py conftest.py
COPY .pytest.ini .pytest.ini

# Open a shell in the container
ENTRYPOINT ["/bin/bash"]