# The vep-rs runtime image. It packages the attested release binaries that
# scripts/docker/stage_release.sh stages into the build context; it never
# compiles them. Built and pushed by .github/workflows/docker.yml.
#
# The release binaries need glibc 2.39: keep both stages on Debian 13 or newer.

FROM debian:trixie-slim@sha256:a99cfc517144bc59b1978475ec53b46ecabec7e43635402ee5b77cc54cd1b20a AS duckdb
# A DuckDB bump changes all three: the sha256s are of the release zip and of
# LICENSE at the release tag.
ARG DUCKDB_VERSION=1.5.5
ARG DUCKDB_SHA256=08c0ca117111fcede14239d0093792352befdc174218c344d232c13279643d05
ARG DUCKDB_LICENSE_SHA256=7e17fd31249fa875cb3b1c5e05c6c3e99b75509f6a2804ca176c217834de1dcb
# hadolint ignore=DL3008
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl unzip \
 && rm -rf /var/lib/apt/lists/*
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
WORKDIR /tmp/duckdb
RUN curl -fsSLo duckdb.zip "https://github.com/duckdb/duckdb/releases/download/v${DUCKDB_VERSION}/duckdb_cli-linux-amd64.zip" \
 && curl -fsSLo LICENSE "https://raw.githubusercontent.com/duckdb/duckdb/v${DUCKDB_VERSION}/LICENSE" \
 && printf '%s  duckdb.zip\n%s  LICENSE\n' "$DUCKDB_SHA256" "$DUCKDB_LICENSE_SHA256" | sha256sum --check --strict - \
 && unzip -q duckdb.zip duckdb

FROM debian:trixie-slim@sha256:a99cfc517144bc59b1978475ec53b46ecabec7e43635402ee5b77cc54cd1b20a
ARG VERSION
LABEL org.opencontainers.image.title="vep-rs" \
      org.opencontainers.image.description="Ensembl Variant Effect Predictor in Rust: vep, vep-cache-builder, vep-cache-converter and the duckdb CLI for Parquet output" \
      org.opencontainers.image.source="https://github.com/natera-open-source/vep-rs" \
      org.opencontainers.image.documentation="https://github.com/natera-open-source/vep-rs#container-image" \
      org.opencontainers.image.licenses="Apache-2.0" \
      org.opencontainers.image.version="${VERSION}"
# ca-certificates: vep-cache-builder downloads from Ensembl over HTTPS.
# libssl3t64: vep-cache-builder links OpenSSL 3. libstdc++6: duckdb.
# hadolint ignore=DL3008
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates libssl3t64 libstdc++6 zlib1g \
 && rm -rf /var/lib/apt/lists/* \
 && useradd --create-home --uid 1000 --user-group vep
COPY --from=duckdb /tmp/duckdb/duckdb /usr/local/bin/duckdb
COPY --from=duckdb /tmp/duckdb/LICENSE /usr/share/doc/duckdb/LICENSE
COPY LICENSE NOTICE /usr/share/doc/vep-rs/
COPY vep vep-cache-builder vep-cache-converter /usr/local/bin/
USER 1000:1000
WORKDIR /data
CMD ["vep", "--help"]
