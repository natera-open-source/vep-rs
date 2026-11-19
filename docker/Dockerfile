# The vep-rs runtime image. It packages the binaries staged into the build
# context and never compiles them: a release's attested ones
# (scripts/docker/stage_release.sh) or, on a pull request, the tree's own
# build by the release's recipe (scripts/docker/stage_tree.sh). Built by
# .github/workflows/docker.yml, once per platform on a runner of that
# architecture; only an image of a release's binaries is pushed.
#
# From either source the binaries need the glibc floor
# .github/workflows/release.yml builds them for (GLIBC_FLOOR) and
# vep-cache-builder OpenSSL 3; any current Debian release meets both.

FROM debian:trixie-slim@sha256:f6b08b09a05725f9431502d3b5714f5c4576e82c1a84fa1debe0c24a7e0c3ec4
ARG VERSION
LABEL org.opencontainers.image.title="vep-rs" \
      org.opencontainers.image.description="Ensembl Variant Effect Predictor in Rust: vep, vep-cache-builder, vep-cache-converter, bgzip, tabix and ps" \
      org.opencontainers.image.source="https://github.com/natera-open-source/vep-rs" \
      org.opencontainers.image.documentation="https://natera-open-source.github.io/vep-rs/getting-started/installation.html#container-image" \
      org.opencontainers.image.licenses="Apache-2.0" \
      org.opencontainers.image.version="${VERSION}"
# ca-certificates: vep-cache-builder downloads from Ensembl over HTTPS.
# libssl3t64: vep-cache-builder links OpenSSL 3.
# procps: ps, which Nextflow's task wrapper calls inside the container.
# tabix: the bgzip and tabix binaries, for indexing --compress_output bgzip
# VCF output.
# hadolint ignore=DL3008
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates libssl3t64 procps tabix zlib1g \
 && rm -rf /var/lib/apt/lists/* \
 && useradd --create-home --uid 1000 --user-group vep
# Gate: the store must never activate these two roots, whichever
# ca-certificates version apt installs. test -s first: `!` also inverts
# grep's exit 2, so the bare negation passes when the file is missing.
RUN test -s /etc/ca-certificates.conf \
 && ! grep -qE '^(mozilla/emSign_ECC_Root_CA_-_G3|mozilla/Comodo_AAA_Services_root)\.crt$' /etc/ca-certificates.conf
COPY LICENSE NOTICE /usr/share/doc/vep-rs/
COPY vep vep-cache-builder vep-cache-converter /usr/local/bin/
USER 1000:1000
WORKDIR /data
CMD ["vep", "--help"]
