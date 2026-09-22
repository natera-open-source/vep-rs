# Run the full test suite
test:
    cargo test --workspace

# Build a release binary
build:
    cargo build --release

# Lint and format check
check:
    cargo clippy --workspace -- -D warnings
    cargo fmt --all -- --check

# Build release and run a small benchmark on the committed golden corpus
bench: build
    @echo "Benchmarking the 1,008-record GRCh37 golden corpus..."
    target/release/vep \
        -i tests/golden/115/GRCh37/variants.vcf \
        -o /dev/null \
        --offline --json_cache tests/golden/115/GRCh37/json_cache --assembly GRCh37 \
        --fork 1 --buffer_size 5000 --force --no_stats --quiet

# Run smoke concordance test (5000 variants). The harness needs the Perl VEP cache
# root, a converted JSON cache and an indexed FASTA; pass them through the variables.
concordance-smoke perl_cache json_cache fasta:
    scripts/concordance/run_concordance.sh --mode smoke --smoke-variants 5000 \
        --perl-cache-dir {{perl_cache}} --json-cache-dir {{json_cache}} --fasta {{fasta}}

# Generate and open workspace docs
doc:
    cargo doc --workspace --no-deps --open
