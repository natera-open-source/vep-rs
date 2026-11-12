# Nextflow and nf-core

nf-core pipelines such as sarek and raredisease annotate with the nf-core `ENSEMBLVEP_VEP` module. vep-rs runs that module's command line unchanged: the container image takes the place of Ensembl VEP's, and the module's script, flags and outputs stay as they are. This page states what the swap covers and what it leaves out.

## What the module runs

The module calls `vep` with `-i`, `-o <prefix>.vcf.gz`, the pipeline's own arguments (`task.ext.args`), `--compress_output bgzip`, `--assembly`, `--species`, `--cache --cache_version <version> --dir_cache <cache>` and `--fork <cpus>`, then indexes the output with `tabix`. vep-rs reads that command line as Ensembl VEP does. It writes the VCF Ensembl VEP 116.2 writes for the same records, BGZF-compressed so `tabix` indexes it.

Sarek's default arguments, `--stats_file`, `--vcf`, `--everything`, `--filter_common`, `--per_gene`, `--total_length`, `--offline` and `--format vcf`, run through the same script. The image carries `vep`, `bgzip` and `tabix`, which are the commands the module's script and its version probes call.

## Running a pipeline on vep-rs

Point the module at the image in a Nextflow configuration file, passed with `-c`:

```groovy
process {
    withName: 'ENSEMBLVEP_VEP' {
        container = 'ghcr.io/natera-open-source/vep-rs:0.5.0'
    }
}
```

The cache the pipeline stages, sarek's `--vep_cache` for one, must be a vep-rs cache in Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout. Build it or convert an Ensembl VEP cache into that layout as [Cache setup](../cache-setup.md) describes, and pass the pipeline's cache version and assembly as before.

## What the swap leaves out

These parts of the module do not carry over to vep-rs:

- **The `perl-math-cdf` version probe.** The module records the version of Perl's `Math::CDF` with `perl -MMath::CDF`. The vep-rs image does not carry `Math::CDF`, so Nextflow fails the task with `Unable to evaluate output` after the annotation has been written. Run the module patched without its `versions_perlmathcdf` output (`nf-core modules patch ensemblvep/vep`), or from an image built on vep-rs's that adds `Math::CDF`.
- **The `ensemblvep` version.** The module reads the version from the `ensembl-vep` line of `vep --help`. vep-rs prints no such line, so the pipeline's version report shows an empty `ensemblvep` version.
- **The summary report.** `--stats_file` is accepted and no summary file is written, so the module emits no `report` and MultiQC shows no VEP section.
- **GTF input.** The module's `gtf` input runs `vep --gtf`, which vep-rs does not accept; the module needs a cache.
- **Plugins.** A plugin argument names one of vep-rs's built-in plugins or a dynamic library in `--dir_plugins` ([Plugins](../plugins.md)). A plugin that cannot be loaded is a warning and the run goes on without it, or an error under `--safe`.
- **An nf-core module of its own.** nf-core modules install their tool from Bioconda, and vep-rs has no Bioconda package yet ([What vep-rs does not do](what-vep-rs-does-not-do.md#bioconda)).
