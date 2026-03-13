# dockerfile for chromosome_loc_annotator

FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Working directory inside container
WORKDIR /app

# Install basic build tooling (usually cheap and helpful)
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy project metadata and sources into the image
# (assumes you build with repository root as context)
COPY ./pyproject.toml README.rst requirements.txt ./
COPY ./chromloc ./chromloc
COPY ./chromloc/non_alt_loci_set.json ./chromloc/non_alt_loci_set.json

# Install the package into the container
RUN pip install --upgrade pip \
    && pip install . \
    && pip cache purge

# Default entrypoint is the console script we defined in pyproject.toml
ENTRYPOINT ["cytoscape-chromloc"]

# By default just show help if no args given
CMD ["-h"]
