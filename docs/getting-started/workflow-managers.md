# Workflow managers

Run vep-rs in place of Ensembl VEP in a pipeline that calls it through the nf-core module, Ensembl's own Nextflow pipeline or the Snakemake wrapper.

Each section names the command the wrapper issues, what vep-rs does with every flag of it, and the one or two lines that change.

## The nf-core module

To run the module on vep-rs, override its container in the pipeline's Nextflow configuration:

```groovy
process {
    withName: 'ENSEMBLVEP_VEP' {
        container = 'ghcr.io/natera-open-source/vep-rs:X.Y.Z'
        ext.args  = '--vcf --everything'
    }
}
```

The nf-core module `ensemblvep/vep` runs one command per sample:

```bash
vep -i <sample>.vcf -o <prefix>.vcf.gz <ext.args> --compress_output bgzip \
  --fasta <fasta> --assembly <assembly> --species <species> \
  --cache --cache_version <version> --dir_cache <cache> --fork <cpus>
```

Every flag of that line is honoured by vep-rs: `--dir_cache` with `--species` and `--cache_version` resolves the cache in Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout, `--compress_output bgzip` writes the BGZF stream `tabix` indexes, and `--fork` is the thread count.

The cache the module stages must be a vep-rs cache: `vep-cache-builder --auto <species> <assembly> --release <version> --dir_cache <cache>` writes one in that layout, where the module's own cache input is Ensembl's Storable archive, which `vep` does not read.

The module's `ext.args` carry the output format and the annotation flags as for Ensembl VEP; a `--gtf` input in place of the cache is refused, because vep-rs reads a GTF only as a cache built from it ([A GTF or GFF3 as the transcript model](../cache-from-gtf-or-gff3.md)).

The module's `versions` outputs run `vep --help`, `tabix -h` and `perl -MMath::CDF` inside the container and read a version from each; the vep-rs image answers the `tabix` line alone. `vep --help` prints no `ensembl-vep` line, so the module's `ensemblvep` version reads empty. The image carries no `Math::CDF`, so the `perl -MMath::CDF` probe fails, and Nextflow fails the task with `Unable to evaluate output` after the annotation has been written. A copy of the module patched with `nf-core modules patch ensemblvep/vep` drops the `versions_perlmathcdf` output, and reads the `ensemblvep` version from `vep --version` or drops that output too.

## Ensembl's Nextflow pipeline

Ensembl VEP ships a Nextflow pipeline under `nextflow/`. It splits a VCF and annotates each part with `vep -i <part> -o out.vcf --vcf --config <vep.ini>`. When `--filters` is given, it filters the part before sorting, compressing and indexing it. The filter command is `filter_vep --format vcf -i vep.raw.vcf -o out.vcf --force_overwrite --only_matched -filter "<expression>"`.

vep-rs reads the same `--config` file, one flag per line as `<name> <value>` without the dashes, and `vep filter` takes `filter_vep`'s options under the same spellings ([vep filter](../cli/filter.md)). Two lines change: the pipeline's `vep` label takes the vep-rs image as its container (`withLabel: vep { container = 'ghcr.io/natera-open-source/vep-rs:X.Y.Z' }`) and its `filter_vep` command becomes `vep filter`, since the image carries no executable of that name. The pipeline's `--vep_version` option, which picks an `ensemblorg/ensembl-vep` tag from Docker Hub, is left unset once the container is overridden.

## The Snakemake wrapper

The snakemake-wrappers rule `bio/vep/annotate` pipes the calls through `bcftools view` into:

```bash
vep <extra> --fork <threads> --format vcf --vcf \
  --offline --cache --dir_cache <cache> --cache_version <release> \
  --species <species> --assembly <build> --fasta <fasta> \
  --dir_plugins <plugins> --plugin <Name>,<file> \
  --output_file STDOUT --warning_file STDERR --safe --stats_file <stats>
```

and back through `bcftools view` into the compressed output. vep-rs honours each flag: `--output_file STDOUT` writes the annotation on standard output, `--warning_file STDERR` keeps the warnings on standard error, `--safe` makes a plugin that cannot be loaded an error, and `--stats_file` names the statistics report.

The wrapper reads the species, release and assembly from the cache directory's layout, `<cache>/<species>/<release>_<build>`. A cache built with `vep-cache-builder --auto` into the rule's cache directory is therefore found as the wrapper expects, and that build replaces the companion rule `bio/vep/cache`, which downloads Ensembl's Storable archive.

A `--plugin` name among the built-in plugins ([Plugins](../plugins.md)) runs built in with the data file the rule names. The Perl modules the `bio/vep/plugins` rule downloads into `--dir_plugins` are not loaded, so a name that is neither built in nor a shared library there is refused under the wrapper's `--safe`.

The wrapper runs in a conda environment holding Ensembl VEP; vep-rs has no Bioconda package ([What vep-rs does not do](what-vep-rs-does-not-do.md#bioconda)), so the rule runs with the `vep` binary on `PATH`, or in the vep-rs container, in place of that environment.

## Where next

- [CLI reference](../cli-reference.md): the flags themselves.
- [Cache setup](../cache-setup.md): building the cache a wrapper stages.
- [Plugins](../plugins.md): the built-in plugins a `--plugin` name runs.
- [vep filter](../cli/filter.md): the command that takes `filter_vep`'s options.
