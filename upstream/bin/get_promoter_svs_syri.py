#!/usr/bin/env python3
"""
Find SVs (SyRI VCF format) where the reference coordinates overlap the promoter
of one gene AND the query coordinates overlap the promoter of a different gene,
on the same chromosome.

Filters:
  - INVAL  (<INVAL>)
  - NOTAL  (<NOTAL>)
  - Cross-chromosome SVs (ChrB base name != CHROM base name, or ChrB == ".")

Haplotype suffixes (_1, _2) are stripped before the same-chromosome check so
that chr_01_1 (hap1) and chr_01_2 (hap2) are treated as the same chromosome.
GFF3 lookups use the full suffixed names directly.

No tabix index needed — reads the VCF line-by-line (plain or bgzipped).

Usage:
    python get_promoter_svs_syri.py annotation.gff3 variants.vcf[.gz] \\
        -o promoter_sv_pairs.tsv --upstream 2000
"""

import sys
import gzip
import argparse
import bisect
from collections import defaultdict


# ---------------------------------------------------------------------------
# GFF3 parsing
# ---------------------------------------------------------------------------

def _parse_attrs(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def parse_all_genes(gff_file: str, upstream_bp: int = 2000, downstream_bp: int = 200) -> dict:
    """
    Parse gene features from GFF3.

    Returns
    -------
    dict  chrom -> sorted list of
          (prom_start, prom_end, gene_id, tss, strand)  [1-based inclusive]
    """
    genes_by_chrom = defaultdict(list)

    with open(gff_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2].lower() != "gene":
                continue

            attrs   = _parse_attrs(fields[8])
            gene_id = attrs.get("ID", "").replace("gene:", "")
            if not gene_id:
                gene_id = attrs.get("Name", "")
            if not gene_id:
                continue

            chrom  = fields[0]
            start  = int(fields[3])   # 1-based
            end    = int(fields[4])   # 1-based
            strand = fields[6]
            tss    = start if strand == "+" else end

            if strand == "+":
                prom_start = max(1, tss - upstream_bp)
                prom_end   = tss + downstream_bp
            else:
                prom_start = max(1, tss - downstream_bp)
                prom_end   = tss + upstream_bp

            genes_by_chrom[chrom].append(
                (prom_start, prom_end, gene_id, tss, strand)
            )

    # Sort each chromosome's list by promoter start; cache start positions for bisect
    result = {}
    for chrom, lst in genes_by_chrom.items():
        sorted_lst = sorted(lst, key=lambda x: x[0])
        prom_starts = [g[0] for g in sorted_lst]
        result[chrom] = (sorted_lst, prom_starts)

    return result


def find_overlapping_genes(genes_by_chrom: dict, chrom: str, start: int, end: int) -> list:
    """
    Return all genes whose promoter overlaps [start, end] (1-based inclusive).

    Returns list of (gene_id, tss, strand, prom_start, prom_end).
    Uses bisect on sorted prom_starts to skip genes where prom_start > end.
    """
    if chrom not in genes_by_chrom:
        return []

    genes, prom_starts = genes_by_chrom[chrom]
    # All genes with prom_start <= end can possibly overlap; skip the rest
    idx = bisect.bisect_right(prom_starts, end)
    hits = []
    for prom_start, prom_end, gene_id, tss, strand in genes[:idx]:
        if prom_end >= start:
            hits.append((gene_id, tss, strand, prom_start, prom_end))
    return hits


# ---------------------------------------------------------------------------
# VCF helpers
# ---------------------------------------------------------------------------

# ALT values to unconditionally skip
_ALT_SKIP = {"<INVAL>", "<NOTAL>"}
# 
_ALT_SKIP = {}

def _is_filtered(alt: str, ref: str, sv_id: str) -> bool:
    """Return True if this record should be excluded."""
    # Skip INVAL and NOTAL
    if alt in _ALT_SKIP:
        return True
    return False


def _parse_info(info_str: str) -> dict:
    info = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            info[k] = v
        else:
            info[part] = True
    return info


def _sv_type(alt: str, sv_id: str) -> str:
    """Derive a human-readable SV type."""
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1]          # e.g. INV, DEL, HDR, INVDP …
    # Literal-sequence variant — infer from ID prefix
    for prefix in ("SNP","SNV", "DEL", "INS", "DUP", "INV", "HDR", "TDM"):
        if sv_id and sv_id.upper().startswith(prefix):
            return prefix
    return "INDEL"


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _base_chrom(chrom: str) -> str:
    """Strip haplotype suffix (_1 / _2) for same-chromosome comparison."""
    if len(chrom) >= 2 and chrom[-1].isdigit() and chrom[-2] == "_":
        return chrom[:-2]
    return chrom


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

def _signed_dist_tss(sv_start: int, sv_end: int, tss: int, strand: str) -> int:
    """
    Signed distance from SV to TSS.
    Positive = upstream of TSS, negative = downstream of TSS (into gene body).
    Zero when TSS falls inside the SV.
    """
    if sv_start <= tss <= sv_end:
        return 0
    if strand == "+":
        if sv_end < tss:
            return tss - sv_end   # SV is upstream
        return -(sv_start - tss)  # SV is downstream
    else:
        if sv_start > tss:
            return sv_start - tss  # SV is upstream
        return -(tss - sv_end)     # SV is downstream

def find_promoter_sv_pairs(
    gff_file: str,
    vcf_file: str,
    upstream_bp: int = 2000,
    ref_downstream_bp: int = 500,
    query_upstream_bp: int = 500,
    query_downstream_bp: int = 500,
    ref_suffix: str = "",
    query_suffix: str = "",
    max_sv_len: int = 2500,
) -> list:
    """
    Main function — returns list of dicts, one per (SV, ref_gene, query_gene) triple.

    The SV must overlap the promoter (upstream_bp) of the ref gene.
    On the query side the SV only needs to fall within query_upstream_bp of any gene's TSS.

    Each dict contains:
        sv_id, sv_type, chrom,
        ref_start, ref_end, query_start, query_end,
        ref_gene, ref_tss, ref_strand, ref_promoter,
        query_gene, query_tss, query_strand, query_promoter
    """
    print(f"[info] Parsing GFF3: {gff_file}", file=sys.stderr)
    ref_genes_by_chrom   = parse_all_genes(gff_file, upstream_bp, ref_downstream_bp)
    query_genes_by_chrom = parse_all_genes(gff_file, query_upstream_bp, query_downstream_bp)
    n_genes = sum(len(genes) for genes, _ in ref_genes_by_chrom.values())
    print(f"[info] {n_genes} genes on {len(ref_genes_by_chrom)} chromosomes.", file=sys.stderr)
    print(f"[info] Ref promoter window : -{upstream_bp:,} / +{ref_downstream_bp:,} bp from TSS", file=sys.stderr)
    print(f"[info] Query promoter window : -{query_upstream_bp:,} / +{query_downstream_bp:,} bp from TSS", file=sys.stderr)
    if ref_suffix or query_suffix:
        print(f"[info] Chrom mapping: VCF CHROM + '{ref_suffix}' -> GFF3 ref; "
              f"VCF ChrB + '{query_suffix}' -> GFF3 query", file=sys.stderr)

    results        = []
    n_total        = 0
    n_filtered_snp = 0
    n_filtered_xc  = 0  # cross-chromosome
    n_filtered_dot = 0  # missing query coords
    n_filtered_len = 0  # too long
    n_hits         = 0

    print(f"[info] Scanning VCF: {vcf_file}", file=sys.stderr)
    with _open(vcf_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom    = fields[0]
            pos      = int(fields[1])    # 1-based
            sv_id    = fields[2]
            ref      = fields[3]
            alt      = fields[4]
            info_str = fields[7]

            n_total += 1

            # ── Filter 1: SNP / INVAL / NOTAL ──────────────────────────────
            if _is_filtered(alt, ref, sv_id):
                n_filtered_snp += 1
                continue

            info    = _parse_info(info_str)
            end     = int(info.get("END", pos))
            chrb    = info.get("ChrB", ".")
            sv_type = _sv_type(alt, sv_id)
            if alt.startswith("<") and alt.endswith(">"):
                # Symbolic allele — use SVLEN if present, else ref span
                sv_len = int(info.get("SVLEN", end - pos))
            else:
                # Literal sequence — insertion length = len(ALT) - len(REF)
                sv_len = abs(len(alt) - len(ref))

            # ── Filter 2: Max variant length (inversions only) ──────────────
            if sv_type == "INV" and sv_len > max_sv_len:
                n_filtered_len += 1
                continue

            # ── Filter 3: Same-chromosome only (ignoring haplotype suffix) ────
            if chrb == "." or _base_chrom(chrb) != _base_chrom(chrom):
                n_filtered_xc += 1
                continue

            startb_str = info.get("StartB", ".")
            endb_str   = info.get("EndB",   ".")

            # ── Filter 4: Must have query coordinates ──────────────────────
            if startb_str == "." or endb_str == ".":
                n_filtered_dot += 1
                continue

            startb = int(startb_str)
            endb   = int(endb_str)
            # Normalize — StartB > EndB for inverted regions
            qstart = min(startb, endb)
            qend   = max(startb, endb)

            # ── Promoter overlaps ──────────────────────────────────────────
            ref_genes   = find_overlapping_genes(ref_genes_by_chrom,   chrom + ref_suffix,   pos,    end)
            query_genes = find_overlapping_genes(query_genes_by_chrom, chrb  + query_suffix, qstart, qend)

            if not ref_genes or not query_genes:
                continue

            # Emit one row per (ref_gene, query_gene) pair, skip self-pairs
            for rg_id, rg_tss, rg_strand, rg_ps, rg_pe in ref_genes:
                ref_dist = _signed_dist_tss(pos, end, rg_tss, rg_strand)
                ref_overlap_len = min(end, rg_pe) - max(pos, rg_ps) + 1
                for qg_id, qg_tss, qg_strand, qg_ps, qg_pe in query_genes:
                    if rg_id == qg_id:
                        continue
                    query_dist = _signed_dist_tss(qstart, qend, qg_tss, qg_strand)
                    query_overlap_len = min(qend, qg_pe) - max(qstart, qg_ps) + 1
                    n_hits += 1
                    results.append({
                        "sv_id":           sv_id,
                        "sv_type":         sv_type,
                        "sv_len":          sv_len if sv_type == "INS" else ref_overlap_len,
                        "query_sv_len":    query_overlap_len,
                        "chrom":           chrom,
                        "ref_start":       pos,
                        "ref_end":         end,
                        "query_start":     startb,
                        "query_end":       endb,
                        "ref_gene":        rg_id,
                        "ref_tss":         rg_tss,
                        "ref_strand":      rg_strand,
                        "ref_promoter":    f"{rg_ps}-{rg_pe}",
                        "ref_dist_tss":    ref_dist,
                        "query_gene":      qg_id,
                        "query_tss":       qg_tss,
                        "query_strand":    qg_strand,
                        "query_promoter":  f"{qg_ps}-{qg_pe}",
                        "query_dist_tss":  query_dist,
                    })

    print(f"[info] VCF records read       : {n_total}",        file=sys.stderr)
    print(f"[info] Filtered (INVAL/NOTAL) : {n_filtered_snp}", file=sys.stderr)
    print(f"[info] Filtered (INV len >{max_sv_len:,}bp) : {n_filtered_len}", file=sys.stderr)
    print(f"[info] Filtered (cross-chrom) : {n_filtered_xc}",  file=sys.stderr)
    print(f"[info] Filtered (no coords)   : {n_filtered_dot}", file=sys.stderr)
    print(f"[info] Promoter-pair hits     : {n_hits}",         file=sys.stderr)
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Report SVs (SyRI VCF) where reference coordinates overlap one gene's "
            "promoter and query coordinates overlap a different gene's promoter "
            "on the same chromosome."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("gff_file", help="GFF3 annotation file")
    parser.add_argument("vcf_file", help="SyRI VCF file (plain .vcf or bgzipped .vcf.gz)")
    parser.add_argument(
        "--upstream", type=int, default=2000,
        help="Promoter window size in bp upstream of TSS (ref side)",
    )
    parser.add_argument(
        "--ref-downstream", type=int, default=500,
        help="Extension in bp downstream of TSS for the ref promoter window (default: 500)",
    )
    parser.add_argument(
        "--query-upstream", type=int, default=10000,
        help="Window size in bp upstream of TSS on the query side (default: 10000)",
    )
    parser.add_argument(
        "--query-downstream", type=int, default=200,
        help="Extension in bp downstream of TSS for the query promoter window (default: 200)",
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output TSV (default: stdout)",
    )
    parser.add_argument(
        "--max-sv-len", type=int, default=2000,
        help="Maximum length in bp (END - POS) for INV variants; longer inversions are excluded (default: 2000)",
    )
    parser.add_argument(
        "--ref-suffix", default="",
        help="Suffix appended to VCF CHROM to match GFF3 chromosome names "
             "(e.g. '_1' when GFF3 has chr_01_1 but VCF has chr_01)",
    )
    parser.add_argument(
        "--query-suffix", default="",
        help="Suffix appended to VCF ChrB to match GFF3 chromosome names "
             "(e.g. '_2' for hap2)",
    )
    args = parser.parse_args()

    rows = find_promoter_sv_pairs(
        args.gff_file, args.vcf_file, args.upstream,
        ref_downstream_bp=args.ref_downstream,
        query_upstream_bp=args.query_upstream,
        query_downstream_bp=args.query_downstream,
        ref_suffix=args.ref_suffix, query_suffix=args.query_suffix,
        max_sv_len=args.max_sv_len,
    )

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        header = [
            "sv_id", "sv_type", "sv_len", "query_sv_len", "chrom",
            "ref_start", "ref_end", "query_start", "query_end",
            "ref_gene", "ref_tss", "ref_strand", "ref_promoter", "ref_dist_tss",
            "query_gene", "query_tss", "query_strand", "query_promoter", "query_dist_tss",
        ]
        print("\t".join(header), file=out)
        for r in rows:
            print(
                f"{r['sv_id']}\t{r['sv_type']}\t{r['sv_len']}\t{r['query_sv_len']}\t{r['chrom']}\t"
                f"{r['ref_start']}\t{r['ref_end']}\t"
                f"{r['query_start']}\t{r['query_end']}\t"
                f"{r['ref_gene']}\t{r['ref_tss']}\t{r['ref_strand']}\t{r['ref_promoter']}\t{r['ref_dist_tss']}\t"
                f"{r['query_gene']}\t{r['query_tss']}\t{r['query_strand']}\t{r['query_promoter']}\t{r['query_dist_tss']}",
                file=out,
            )
    finally:
        if args.output:
            out.close()

    print(f"[info] Done — {len(rows)} rows written.", file=sys.stderr)
