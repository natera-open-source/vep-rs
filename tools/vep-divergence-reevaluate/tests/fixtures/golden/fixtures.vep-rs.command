{{ROOT}}/fake-vep-rs -i {{ROOT}}/work/still_present_fixture/stage/input/record.vcf -o {{ROOT}}/work/still_present_fixture/out/vep-rs/default.txt --offline --json_cache {{ROOT}}/work/still_present_fixture/stage/json_cache --species homo_sapiens --assembly GRCh38 --format vcf --buffer_size 5000 --force_overwrite --quiet --no_stats --minimal --no_intergenic
exit 0
