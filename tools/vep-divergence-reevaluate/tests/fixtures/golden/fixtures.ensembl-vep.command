{{ROOT}}/fake-ensembl-vep -i {{ROOT}}/work/still_present_fixture/stage/input/record.vcf -o {{ROOT}}/work/still_present_fixture/out/ensembl-vep/default.txt --offline --cache --dir_cache {{ROOT}}/perl --cache_version 116 --species homo_sapiens --assembly GRCh38 --format vcf --buffer_size 5000 --fork 1 --force_overwrite --quiet --no_stats --minimal --no_intergenic
exit 0
