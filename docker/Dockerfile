FROM python:3.11-bookworm

SHELL ["/bin/bash", "-c"]
WORKDIR /app

# System deps (libgl1 + libxt6 needed by VTK / IMOD)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget libgl1 libxt6 && \
    rm -rf /var/lib/apt/lists/*

# IMOD 4.11.25 (required for TEM simulation: xyzproj, tilt, alterheader)
#
# The IMOD .sh installer derives its version from its own filename
# (via sed on imod_[0-9.]*).  The filename MUST be preserved as-is;
# renaming it breaks the version extraction and causes a tar path
# mismatch ("imod_/" vs "imod_4.11.25/").
RUN wget https://bio3d.colorado.edu/imod/AMD64-RHEL5/imod_4.11.25_RHEL7-64_CUDA10.1.sh \
        --no-check-certificate -O /tmp/imod_4.11.25_RHEL7-64_CUDA10.1.sh && \
    chmod +x /tmp/imod_4.11.25_RHEL7-64_CUDA10.1.sh && \
    mkdir -p /opt/imod && \
    /tmp/imod_4.11.25_RHEL7-64_CUDA10.1.sh -y -dir /opt/imod && \
    rm -f /tmp/imod_4.11.25_RHEL7-64_CUDA10.1.sh && \
    # Verify the three required IMOD binaries
    for bin in xyzproj tilt alterheader; do \
        test -x /opt/imod/IMOD/bin/"$bin" || \
            { echo "FATAL: missing IMOD binary: $bin" >&2; exit 1; }; \
    done

# IMOD environment — set up via entrypoint so IMOD-linux.sh is sourced.
ENV IMOD_DIR=/opt/imod/IMOD

# Writable HOME for non-root --user (matplotlib, fontconfig caches)
ENV HOME=/tmp
ENV MPLCONFIGDIR=/tmp/matplotlib

# Python package — install core polnet
COPY pyproject.toml README.md LICENSE NOTICE ./
COPY src/ src/
RUN pip install --no-cache-dir .

# Copy input data and example configs
COPY data/ data/
COPY config/ config/

# Entrypoint wrapper: source IMOD env, then exec polnet CLI
RUN echo '#!/bin/bash' > /app/entrypoint.sh && \
    echo 'source /opt/imod/IMOD/IMOD-linux.sh' >> /app/entrypoint.sh && \
    echo 'exec polnet "$@"' >> /app/entrypoint.sh && \
    chmod +x /app/entrypoint.sh

ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["/app/config/all_features.yaml"]