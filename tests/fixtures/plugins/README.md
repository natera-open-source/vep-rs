# Bundled plugin score files

The score files Ensembl VEP's LoFtool and pLI plugins read from the directory
holding their modules when a run names no file, committed here so the tests of
the two built-ins and the plugin corpora read the same bytes Ensembl VEP reads.

| File                        | Content                                                                                                                                                                           | Source and terms                                                                                                           |
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `LoFtool_scores.txt`        | Gene symbol and LoFtool percentile, tab-separated with CRLF line ends, as bundled; the first line is the header `Gene LoFtool_percentile`                                         | [Ensembl VEP_plugins](https://github.com/Ensembl/VEP_plugins), release 116; Apache License 2.0                             |
| `pLI_values.txt`            | Gene symbol and pLI, tab-separated, as bundled; the first line is the header `gene pLI`                                                                                           | [Ensembl VEP_plugins](https://github.com/Ensembl/VEP_plugins), release 116; Apache License 2.0                             |
| `pLI_transcript_values.txt` | Transcript stable id and pLI, tab-separated under the header `transcript pLI`: the chromosome 21 rows of the ExAC release 0.3 constraint table, the id without its version suffix | the ExAC release 0.3 `functional_gene_constraint` table (Lek et al. 2016), redistributed under the ExAC data release terms |

The transcript file is the extract the pLI plugin's documentation describes
(`transcript` and `pLI` columns of the constraint table), cut to one
chromosome so that a corpus of chromosome 21 records exercises the plugin's
transcript mode. The version suffix is dropped because Ensembl VEP's plugin
looks a row up by the transcript's stable id, which carries none.

A corpus under `tests/golden/116/` that runs one of these plugins commits its
own cut of the file under `plugin_data/`, holding the rows of the genes or
transcripts its records reach, so the corpus stays small and Ensembl VEP's run
over the bundled file and vep-rs's run over the cut write the same rows.
