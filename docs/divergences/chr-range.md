# A range in `--chr`

A difference in which records are annotated under one flag value; nothing is set aside, since no measured run names the flag.

Ensembl VEP splits the `--chr` value on commas and keeps a record when its chromosome, as the input spells it, equals one of the tokens (`Parser.pm` lines 470-472; `Config.pm` lines 760-762). The token `1-5` of the example its own help gives (`Config.pm` line 121) matches a chromosome named `1-5` and no other, so `--chr 1-5,10,MT` annotates the records of `10` and `MT` alone.

vep-rs reads a token of two non-negative integers joined by `-` as every integer between them, so the same value annotates chromosomes 1 to 5 as well. A chromosome literally named `1-5`, which Ensembl VEP admits under that token, is refused silently here, as a record `--chr` leaves out is in both engines. Every other token is matched exactly in both engines: `21` does not admit `chr21`, and a synonym is not consulted.
