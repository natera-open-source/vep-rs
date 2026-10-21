# The vep-rs runtime image. It packages the attested release binaries that
# scripts/docker/stage_release.sh stages into the build context; it never
# compiles them. Built and pushed by .github/workflows/docker.yml, once per
# platform on a runner of that architecture, so TARGETARCH is the build host's.
#
# The binaries need the glibc floor .github/workflows/release.yml builds them
# for (GLIBC_FLOOR) and vep-cache-builder OpenSSL 3; any current Debian release
# meets both.

FROM debian:trixie-slim@sha256:a29215f6a35e51e22adffa17f89e9d2ef06214e64a2bad10d765c46aea49f11f AS duckdb
ARG TARGETARCH
# A DuckDB bump changes all four: the sha256s are of the two release tarballs
# and of LICENSE at the release tag.
ARG DUCKDB_VERSION=2.0.0
ARG DUCKDB_SHA256_AMD64=c6507206baf7793c6013dbd82368469a45f20dc0e17a6f92f204602f652ed2bc
ARG DUCKDB_SHA256_ARM64=859e1875d2a9515fbb4bf8871061437a14af82320da68c753f094c4b776348fe
ARG DUCKDB_LICENSE_SHA256=075c33400ffcb0c586dd106a029d3e733e9da3693f8fcb9cebfc651c8a2f14e3
# hadolint ignore=DL3008
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl \
 && rm -rf /var/lib/apt/lists/*
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
WORKDIR /tmp/duckdb
RUN case "$TARGETARCH" in \
      amd64) sha="$DUCKDB_SHA256_AMD64" ;; \
      arm64) sha="$DUCKDB_SHA256_ARM64" ;; \
      *) echo "ERROR: [docker] no DuckDB CLI for $TARGETARCH" >&2; exit 1 ;; \
    esac \
 && curl -fsSLo duckdb.tar.gz "https://github.com/duckdb/duckdb/releases/download/v${DUCKDB_VERSION}/duckdb-cli-linux-${TARGETARCH}.tar.gz" \
 && curl -fsSLo LICENSE "https://raw.githubusercontent.com/duckdb/duckdb/v${DUCKDB_VERSION}/LICENSE" \
 && printf '%s  duckdb.tar.gz\n%s  LICENSE\n' "$sha" "$DUCKDB_LICENSE_SHA256" | sha256sum --check --strict - \
 && tar -xzf duckdb.tar.gz duckdb

FROM debian:trixie-slim@sha256:a29215f6a35e51e22adffa17f89e9d2ef06214e64a2bad10d765c46aea49f11f
ARG VERSION
LABEL org.opencontainers.image.title="vep-rs" \
      org.opencontainers.image.description="Ensembl Variant Effect Predictor in Rust: vep, vep-cache-builder, vep-cache-converter, bgzip, tabix and the duckdb CLI for Parquet output" \
      org.opencontainers.image.source="https://github.com/natera-open-source/vep-rs" \
      org.opencontainers.image.documentation="https://natera-open-source.github.io/vep-rs/getting-started/installation.html#container-image" \
      org.opencontainers.image.licenses="Apache-2.0" \
      org.opencontainers.image.version="${VERSION}"
# ca-certificates: vep-cache-builder downloads from Ensembl over HTTPS.
# libssl3t64: vep-cache-builder links OpenSSL 3. libstdc++6: duckdb.
# tabix: the bgzip and tabix binaries, for indexing --compress_output bgzip
# VCF output.
# hadolint ignore=DL3008
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates libssl3t64 libstdc++6 tabix zlib1g \
 && rm -rf /var/lib/apt/lists/* \
 && useradd --create-home --uid 1000 --user-group vep
# eMudhra (October 2026): its intermediates were issued from near-zero-entropy
# keys, so anyone on-path can forge a certificate for any domain, and OpenSSL
# clients never consult revocation. Distrust the two roots those chains hang
# from until trixie ships ca-certificates 20250419+deb13u1 (Debian bug #1151413,
# https://emudhra.we-have-a-problem.com/). Drop this once it lands: #69.
RUN sed -i \
      -e 's|^mozilla/emSign_ECC_Root_CA_-_G3.crt|!&|' \
      -e 's|^mozilla/Comodo_AAA_Services_root.crt|!&|' \
      /etc/ca-certificates.conf \
 && update-ca-certificates
# Gate: the store must never (re)activate those roots, workaround or not.
RUN ! grep -qE '^(mozilla/emSign_ECC_Root_CA_-_G3|mozilla/Comodo_AAA_Services_root)\.crt$' /etc/ca-certificates.conf
COPY --from=duckdb /tmp/duckdb/duckdb /usr/local/bin/duckdb
COPY --from=duckdb /tmp/duckdb/LICENSE /usr/share/doc/duckdb/LICENSE
COPY LICENSE NOTICE /usr/share/doc/vep-rs/
COPY vep vep-cache-builder vep-cache-converter /usr/local/bin/
USER 1000:1000
WORKDIR /data
CMD ["vep", "--help"]
