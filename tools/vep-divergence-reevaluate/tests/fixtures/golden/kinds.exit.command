{{ROOT}}/stand-in/ensembl-vep/vep -i {{ROOT}}/work/exit_kind/stage/input/record.vcf -o {{ROOT}}/work/exit_kind/out/ensembl-vep/default.txt --offline --cache --dir_cache {{ROOT}}/perl --cache_version 116 --species homo_sapiens --assembly GRCh38 --format vcf --buffer_size 5000 --fork 1 --force_overwrite --quiet --no_stats
exit 2
