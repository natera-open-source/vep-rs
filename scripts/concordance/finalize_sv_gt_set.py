#!/usr/bin/env python3
"""Finalize a per-VCF SV ground-truth set directory: refresh the canonical-inputs manifest entry
of every regenerated file from the file itself, and write the set's root provenance.json (per-file
digests of every reference output and canonical input, for both assemblies) plus PROVENANCE and
STATUS from the per-file provenance sidecars. Every value is read from the files; nothing is typed.

Usage: finalize_sv_gt_set.py <set-dir> <set-name> <regenerated-basename> <predecessor-set-name> <reason>
                             [--sites-only <assembly>:<basename> ...]

`--sites-only GRCh38:1kg_sv_chr21` names a canonical input rewritten to its eight fixed
columns (sample columns removed) without a reference regeneration: VEP annotates per
position and allele, so its reference output is unchanged, and the entry records the
rewrite. The manifest entry of every regenerated or rewritten input is refreshed from the
file.
"""
import hashlib, json, os, sys, datetime

argv = sys.argv[1:]
sites_only: dict[str, set[str]] = {}
while "--sites-only" in argv:
    i = argv.index("--sites-only")
    asm, base = argv[i + 1].split(":", 1)
    sites_only.setdefault(asm, set()).add(base)
    del argv[i:i + 2]
if len(argv) != 5:
    sys.exit(__doc__)
set_dir, set_name, regen, predecessor, reason = argv

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

def data_lines(p):
    n = 0
    with open(p, "rb") as fh:
        for line in fh:
            if not line.startswith(b"#"):
                n += 1
    return n

def vcf_counts(p):
    headers = variants = 0
    with open(p, "rb") as fh:
        for line in fh:
            if line.startswith(b"#"):
                headers += 1
            else:
                variants += 1
    return headers, variants

per_vcf = {}
sidecars = {}
for asm in ("grch37", "grch38"):
    asm_up = {"grch37": "GRCh37", "grch38": "GRCh38"}[asm]
    gt_dir = os.path.join(set_dir, asm)
    ci_dir = os.path.join(set_dir, "canonical_inputs", asm)
    # canonical-inputs manifest: refresh the regenerated file's entry from the file
    mpath = os.path.join(ci_dir, "manifest.json")
    m = json.load(open(mpath))
    entries = m["files"] if isinstance(m, dict) else m
    refresh = {regen: "regenerated"} | {b: "sites_only" for b in sites_only.get(asm_up, set())}
    hits = {b: 0 for b in refresh}
    for e in entries:
        base = os.path.basename(e["output_file"])[:-4]
        if base in refresh:
            h, v = vcf_counts(os.path.join(ci_dir, base + ".vcf"))
            e["headers_written"], e["variants_written"] = h, v
            e["output_sha256"] = sha256(os.path.join(ci_dir, base + ".vcf"))
            if refresh[base] == "regenerated":
                e["regenerated_by"] = "scripts/validation/generate_sv_test_vcfs.py"
            else:
                e["sample_columns_removed"] = True
            hits[base] += 1
    for base, n in hits.items():
        if n != 1:
            sys.exit(f"ERROR: [finalize_set] {mpath} carries {n} entries for {base}.vcf, expected 1")
    json.dump(m, open(mpath, "w"), indent=2)
    # per-file digests
    per_vcf[asm_up] = {}
    for fn in sorted(os.listdir(gt_dir)):
        if not fn.endswith(".txt") or fn.endswith("_warnings.txt"):
            continue
        base = fn[:-4]
        ci = os.path.join(ci_dir, base + ".vcf")
        if not os.path.isfile(ci):
            sys.exit(f"ERROR: [finalize_set] no canonical input for {asm}/{fn}")
        per_vcf[asm_up][base] = {
            "gt_sha256": sha256(os.path.join(gt_dir, fn)),
            "gt_lines": data_lines(os.path.join(gt_dir, fn)),
            "canonical_input_sha256": sha256(ci),
            "canonical_input_records": vcf_counts(ci)[1],
            "canonical_input_columns": len(open(ci, "rb").read(1 << 26).split(b"#CHROM", 1)[1].split(b"\n", 1)[0].split(b"\t")),
            "regenerated_in_this_set": base == regen,
            "sample_columns_removed_in_this_set": base in sites_only.get(asm_up, set()),
        }
    sc = os.path.join(gt_dir, regen + ".provenance.json")
    sidecars[asm_up] = json.load(open(sc))
    if sidecars[asm_up]["output_sha256"] != per_vcf[asm_up][regen]["gt_sha256"]:
        sys.exit(f"ERROR: [finalize_set] {asm} {regen} sidecar digest differs from the file")

s37, s38 = sidecars["GRCh37"], sidecars["GRCh38"]
for k in ("perl_image", "perl_image_digest", "cache_version", "fork", "buffer_size"):
    if s37[k] != s38[k]:
        sys.exit(f"ERROR: [finalize_set] sidecars disagree on {k}")
now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
root = {
    "set": set_name,
    "generated_at_utc": now,
    "predecessor_set": predecessor,
    "regenerated_files": [regen],
    "canonical_inputs_rewritten_sites_only": sorted(f"{a}:{b}" for a, bs in sites_only.items() for b in sorted(bs)),
    "canonical_inputs_rewritten_note": (
        "sample columns removed with the reference output unchanged: VEP annotates per position and allele "
        "and reads no genotype without --individual, so the rows are those of the genotype-carrying file"
    ) if sites_only else None,
    "reason": reason,
    "perl_docker_image": s37["perl_image"],
    "perl_docker_digest": s37["perl_image_digest"],
    "cache_version": s37["cache_version"],
    "vep_flags": s37["vep_flags"].replace("--assembly GRCh37 ", "--assembly <asm> "),
    "fork": s37["fork"],
    "buffer_size": s37["buffer_size"],
    "regeneration_sidecars": {a: os.path.relpath(os.path.join(set_dir, a.lower(), regen + ".provenance.json"), set_dir) for a in ("GRCh37", "GRCh38")},
    "files_per_assembly": {a: len(v) for a, v in per_vcf.items()},
    "per_vcf": per_vcf,
}
json.dump(root, open(os.path.join(set_dir, "provenance.json"), "w"), indent=2)
tot = {a: sum(v["gt_lines"] for v in per_vcf[a].values()) for a in per_vcf}
with open(os.path.join(set_dir, "PROVENANCE"), "w") as fh:
    fh.write(f"set={set_name}\npredecessor={predecessor}\nregenerated={regen} (both assemblies)\n"
             f"fork={s37['fork']}\nbuffer_size={s37['buffer_size']}\ncache_version={s37['cache_version']}\n"
             f"image={s37['perl_image_digest']}\ngenerated={now}\n"
             f"total_data_lines_grch37={tot['GRCh37']}\ntotal_data_lines_grch38={tot['GRCh38']}\n"
             f"{regen}_data_lines_grch37={per_vcf['GRCh37'][regen]['gt_lines']}\n{regen}_data_lines_grch38={per_vcf['GRCh38'][regen]['gt_lines']}\n"
             f"inputs=canonical_inputs/<asm>/ ({regen}.vcf regenerated; the other {len(per_vcf['GRCh37']) - 1} files copied from the predecessor)\n")
with open(os.path.join(set_dir, "STATUS"), "w") as fh:
    fh.write(f"DONE files={sum(len(per_vcf[a]) for a in per_vcf)} total_grch37={tot['GRCh37']} total_grch38={tot['GRCh38']} {now}\n")
print(f"finalized {set_name}: " + ", ".join(f"{a} {len(per_vcf[a])} files / {tot[a]} rows" for a in per_vcf))
