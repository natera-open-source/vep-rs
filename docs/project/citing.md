# Citing vep-rs

If you use vep-rs in published work, cite the paper:

> Porter M, Borkowski R. _vep-rs: high-throughput Rust variant annotation with population-scale concordance to Ensembl VEP._ bioRxiv (2026). [doi:10.64898/2026.09.22.753614](https://doi.org/10.64898/2026.09.22.753614).

The journal DOI will be added here on acceptance. The paper's figures are in [docs/published-figures.md](../published-figures.md) and re-derive from `manuscript/data/*.csv` in the repository ([manuscript/data/README.md](../../manuscript/data/README.md)). The current release's own measurements are in [Current release](../benchmarks/current-release.md#current-release-v050).

The software itself is archived at Zenodo, and that archive is the citation for the code and the released data files:

> Porter M, Borkowski R. _vep-rs_, version 0.1.0. Zenodo (2026). [doi:10.5281/zenodo.22837897](https://doi.org/10.5281/zenodo.22837897)

`CITATION.cff` at the repository root carries the software citation in machine-readable form, with the all-versions DOI ([10.5281/zenodo.22837896](https://doi.org/10.5281/zenodo.22837896)) and the archived version's own. GitHub's "Cite this repository" button and citation managers pick it up directly. The article is added there on acceptance.

## Acknowledgments

vep-rs is an independent Rust implementation of the variant effect prediction algorithm. It is not affiliated with or endorsed by EMBL-EBI or the Ensembl project.

The algorithm is described in:

> McLaren W, Gil L, Hunt SE, et al. _The Ensembl Variant Effect Predictor._ Genome Biology 17, 122 (2016). [doi:10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4)

## License

Apache-2.0. See [LICENSE](../../LICENSE).
