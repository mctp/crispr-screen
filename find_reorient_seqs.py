#!/usr/bin/env python3
"""
find_reorient_seqs.py

Empirically identify R1/R2 search sequences for fastq reorientation by scanning
FASTQ reads for sgRNA sequences.  Makes no assumptions about scaffold, primer, or
backbone sequences.

Adaptive sgRNA selection (pass 1):
  Tries sgRNAs in small batches against R1.  Stops as soon as accumulated hit reads
  reach --n-target.  Gives up after --max-sgrnas tried without reaching the target.
  High-count sgRNAs in R1 are guaranteed to be well-represented in R2 too (same
  insert, same cluster), so a small number of well-chosen sgRNAs is all that is needed.

Extraction (pass 2):
  Scans R1 and R2 once with the top --n-top sgRNAs by hit count.

Tries to find two sequences per orientation:
  seq 1 — sgRNA-adjacent: seq immediately adjacent to the sgRNA hit. approximately 15-20 nucleotides long.
  seq 2 — read-end proximal: seq at --end-offset from the 5' end. approximately 15-20 nucleotides long.
           Falls back to 3'-proximal if seq 1 is within --min-distance of the 5' end. approximately 15-20 nucleotides long.

Usage (inside Docker):
    python find_reorient_seqs.py \\
        --sgrna-file /repo/sgrnas.txt \\
        --fastq-r1 /input/sample_R1.fastq.gz \\
        [--fastq-r2 /input/sample_R2.fastq.gz]
"""

import argparse
import gzip
import random
import sys
from collections import Counter


def revcomp(seq):
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]


def load_sgrnas(sgrna_file, seed):
    seqs = {}
    with open(sgrna_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2 and all(c in 'ACGTacgtNn' for c in parts[1]):
                name, seq = parts[0], parts[1].upper()
            else:
                name, seq = parts[0], parts[0].upper()
            if len(seq) >= 15 and all(c in 'ACGTN' for c in seq):
                seqs[name] = seq
    if not seqs:
        sys.exit(f"ERROR: No valid sgRNA sequences found in {sgrna_file}")
    names = list(seqs)
    random.seed(seed)
    random.shuffle(names)
    return seqs, names          # dict + shuffled order


def iter_fastq(path, max_reads):
    opener = gzip.open if path.endswith('.gz') else open
    count = 0
    with opener(path, 'rt') as f:
        while True:
            if not f.readline():
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            yield seq
            count += 1
            if max_reads and count >= max_reads:
                break


def count_hits_batch(fastq_path, patterns_dict, max_reads):
    """Scan fastq once; return per-name hit counts for patterns_dict."""
    counts = Counter()
    for seq in iter_fastq(fastq_path, max_reads):
        for name, pat in patterns_dict.items():
            if seq.find(pat) > 0:
                counts[name] += 1
                break
    return counts


def consensus(seq_list, min_frac):
    total = len(seq_list)
    if not total:
        return None, 0, 0
    best, count = Counter(seq_list).most_common(1)[0]
    return (best, count, total) if count / total >= min_frac else (None, count, total)


def extract_at(seq, start, length):
    s = seq[start: start + length]
    return s if len(s) == length else None


def med(lst):
    return sorted(lst)[len(lst) // 2] if lst else None


def show(label, seq, count, total, src):
    pct = 100 * count / total if total else 0
    if seq:
        print(f"\n  {label}")
        print(f"    {seq}   ({count}/{total}, {pct:.1f}%)")
        print(f"    source: {src}")
    else:
        print(f"\n  {label}  —  NOT FOUND   (best: {count}/{total})")
        if src:
            print(f"    attempted: {src}")


def main():
    ap = argparse.ArgumentParser(
        description='Find R1/R2 reorientation search sequences from sgRNA library FASTQs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument('--sgrna-file',        required=True)
    ap.add_argument('--fastq-r1',          required=True)
    ap.add_argument('--fastq-r2',          default=None)
    ap.add_argument('--max-reads',         type=int,   default=500000)
    ap.add_argument('--no-search-revcomp', action='store_true')
    ap.add_argument('--n-target',          type=int,   default=500,
                    help='Target number of hit reads to accumulate before stopping sgRNA search')
    ap.add_argument('--max-sgrnas',        type=int,   default=50,
                    help='Give up after trying this many sgRNAs without reaching --n-target')
    ap.add_argument('--batch-size',        type=int,   default=5,
                    help='sgRNAs per R1 scan batch during adaptive pass 1')
    ap.add_argument('--n-top',             type=int,   default=10,
                    help='Top N sgRNAs by hit count to use for sequence extraction')
    ap.add_argument('--adjacent-len',      type=int,   default=20)
    ap.add_argument('--sgrna-len',         type=int,   default=20)
    ap.add_argument('--end-offset',        type=int,   default=20,
                    help='Offset from read end for proximal seqs')
    ap.add_argument('--window-len',        type=int,   default=18)
    ap.add_argument('--min-distance',      type=int,   default=40,
                    help='Min nt separation required between the two output sequences')
    ap.add_argument('--min-frac',          type=float, default=0.25)
    ap.add_argument('--seed',              type=int,   default=42)
    args = ap.parse_args()

    search_revcomp = not args.no_search_revcomp
    alen = args.adjacent_len
    slen = args.sgrna_len
    eoff = args.end_offset
    wlen = args.window_len

    print(f"Loading sgRNAs from {args.sgrna_file} ...", flush=True)
    all_sgrnas, shuffled_names = load_sgrnas(args.sgrna_file, args.seed)
    print(f"  Loaded {len(all_sgrnas)} sgRNAs  |  search_revcomp: {search_revcomp}")

    def make_r1_pat(seq):
        return revcomp(seq) if search_revcomp else seq
    def make_r2_pat(seq):
        return seq if search_revcomp else revcomp(seq)

    # ---- Pass 1: adaptive — find well-represented sgRNAs ----
    print(f"\nPass 1 — finding well-represented sgRNAs in R1 ...", flush=True)
    cumulative_hits = Counter()
    sgrnas_tried = 0
    total_hits = 0

    for batch_start in range(0, len(shuffled_names), args.batch_size):
        if sgrnas_tried >= args.max_sgrnas or total_hits >= args.n_target:
            break
        batch_names = shuffled_names[batch_start: batch_start + args.batch_size]
        batch_pats  = {n: make_r1_pat(all_sgrnas[n]) for n in batch_names}
        batch_hits  = count_hits_batch(args.fastq_r1, batch_pats, args.max_reads)
        cumulative_hits.update(batch_hits)
        sgrnas_tried += len(batch_names)
        total_hits    = sum(cumulative_hits.values())
        print(f"  Tried {sgrnas_tried} sgRNAs, total hits so far: {total_hits}", flush=True)

    if total_hits == 0:
        sys.exit("ERROR: No hits found in R1. Check --sgrna-file and search_revcomp setting.")

    if total_hits < args.n_target:
        print(f"  WARNING: reached --max-sgrnas ({args.max_sgrnas}) with only {total_hits} hits "
              f"(target {args.n_target}). Proceeding with what was found.", flush=True)

    top_names  = [n for n, _ in cumulative_hits.most_common(args.n_top)]
    top_counts = [cumulative_hits[n] for n in top_names]
    print(f"  Using top {len(top_names)} sgRNAs: hit counts {top_counts[0]}–{top_counts[-1]}")

    top_r1_pats = {n: make_r1_pat(all_sgrnas[n]) for n in top_names}
    top_r2_pats = {n: make_r2_pat(all_sgrnas[n]) for n in top_names}

    # ---- Pass 2: extract candidate sequences ----
    r1_adj=[]; r1_adj_pos=[]
    r1_5p=[];  r1_3p=[];    r1_3p_pos=[]
    r1_nomatch_5p=[]          # 5' of R1 reads that don't match any sgRNA (= correct R2 content)
    r2_adj_up=[]; r2_adj_up_pos=[]
    r2_adj_dn=[]; r2_adj_dn_pos=[]
    r2_5p=[];  r2_3p=[];    r2_3p_pos=[]
    r2_correct_3p=[]; r2_correct_3p_pos=[]   # 3' of correct R2 reads
    r1_hits=r1_total=0
    r2_hits=r2_total=0

    print(f"\nPass 2 — extracting sequences from R1 ...", flush=True)
    for seq in iter_fastq(args.fastq_r1, args.max_reads):
        r1_total += 1
        matched = False
        for pat in top_r1_pats.values():
            pos = seq.find(pat)
            if pos > 0:
                if pos >= alen:
                    r1_adj.append(seq[pos - alen: pos])
                    r1_adj_pos.append(pos - alen)
                w = extract_at(seq, eoff, wlen)
                if w:
                    r1_5p.append(w)
                rlen = len(seq)
                s3 = rlen - eoff - wlen
                if s3 >= 0:
                    r1_3p.append(seq[s3: s3 + wlen])
                    r1_3p_pos.append(s3)
                r1_hits += 1
                matched = True
                break
        if not matched:
            # Non-matching R1 reads are correct R2 content; their 5' seq identifies correct R2 reads
            w = extract_at(seq, eoff, wlen)
            if w:
                r1_nomatch_5p.append(w)
    print(f"  R1 pass 2: {r1_total} reads, {r1_hits} hits ({100*r1_hits/r1_total:.1f}%)")

    r1_nomatch_5p_seq, r1_nomatch_5p_n, r1_nomatch_5p_t = consensus(r1_nomatch_5p, args.min_frac)

    if args.fastq_r2:
        print(f"\nPass 2 — extracting sequences from R2 ...", flush=True)
        for seq in iter_fastq(args.fastq_r2, args.max_reads):
            r2_total += 1
            for pat in top_r2_pats.values():
                pos = seq.find(pat)
                if pos > 0:
                    if pos >= alen:
                        r2_adj_up.append(seq[pos - alen: pos])
                        r2_adj_up_pos.append(pos - alen)
                    dn = pos + slen
                    w = extract_at(seq, dn, alen)
                    if w:
                        r2_adj_dn.append(w)
                        r2_adj_dn_pos.append(dn)
                    w = extract_at(seq, eoff, wlen)
                    if w:
                        r2_5p.append(w)
                    rlen = len(seq)
                    s3 = rlen - eoff - wlen
                    if s3 >= 0:
                        r2_3p.append(seq[s3: s3 + wlen])
                        r2_3p_pos.append(s3)
                    r2_hits += 1
                    break
            # Identify correct R2 reads by matching r1_nomatch_5p_seq near 5' end;
            # collect their 3' end sequence
            if r1_nomatch_5p_seq and r1_nomatch_5p_seq in seq[:eoff + wlen + 5]:
                rlen = len(seq)
                s3 = rlen - eoff - wlen
                if s3 >= 0:
                    r2_correct_3p.append(seq[s3: s3 + wlen])
                    r2_correct_3p_pos.append(s3)
        print(f"  R2 pass 2: {r2_total} reads, {r2_hits} hits ({100*r2_hits/r2_total:.1f}%)")

    # ---- Consensus ----
    r1_adj_seq,r1_adj_n,r1_adj_t = consensus(r1_adj,    args.min_frac)
    r1_5p_seq, r1_5p_n, r1_5p_t  = consensus(r1_5p,    args.min_frac)
    r1_3p_seq, r1_3p_n, r1_3p_t  = consensus(r1_3p,    args.min_frac)
    r2_up_seq, r2_up_n, r2_up_t  = consensus(r2_adj_up, args.min_frac)
    r2_dn_seq, r2_dn_n, r2_dn_t  = consensus(r2_adj_dn, args.min_frac)
    r2_5p_seq, r2_5p_n, r2_5p_t  = consensus(r2_5p,    args.min_frac)
    r2_3p_seq, r2_3p_n, r2_3p_t  = consensus(r2_3p,    args.min_frac)
    r2_c3p_seq, r2_c3p_n, r2_c3p_t = consensus(r2_correct_3p, args.min_frac)

    r1_adj_mpos    = med(r1_adj_pos)
    r1_3p_mpos     = med(r1_3p_pos)
    r2_up_mpos     = med(r2_adj_up_pos)
    r2_dn_mpos     = med(r2_adj_dn_pos)
    r2_3p_mpos     = med(r2_3p_pos)
    r2_c3p_mpos    = med(r2_correct_3p_pos)

    # ---- Select R1 pair ----
    r1s1, r1s1_n, r1s1_t = r1_adj_seq, r1_adj_n, r1_adj_t
    r1s1_pos = r1_adj_mpos
    r1s1_src = f"upstream of sgRNA hit in R1 (~pos {r1s1_pos})"

    if r1s1_pos is not None and abs(r1s1_pos - eoff) >= args.min_distance:
        r1s2, r1s2_n, r1s2_t = r1_5p_seq, r1_5p_n, r1_5p_t
        r1s2_src = f"5' end of read (~pos {eoff})"
    else:
        r1s2, r1s2_n, r1s2_t = r1_3p_seq, r1_3p_n, r1_3p_t
        r1s2_src = f"3' end of read (~pos {r1_3p_mpos}) [5' end too close to seq 1]"

    # ---- Select R2 pair ----
    LOW_R2_HITS = 200
    r2s1, r2s1_n, r2s1_t, r2s1_pos, r2s1_src = None, 0, 0, None, None
    r2s2, r2s2_n, r2s2_t, r2s2_src = None, 0, 0, ""

    if r2_hits >= LOW_R2_HITS:
        # Enough sgRNA hits in R2 — use sgRNA-adjacent sequences
        candidates = []
        if r2_up_seq:
            candidates.append((r2_up_n, r2_up_seq, r2_up_n, r2_up_t, r2_up_mpos, "adjacent to sgRNA in R2"))
        if r2_dn_seq:
            candidates.append((r2_dn_n, r2_dn_seq, r2_dn_n, r2_dn_t, r2_dn_mpos, "adjacent to sgRNA in R2"))
        if candidates:
            candidates.sort(reverse=True)
            _, r2s1, r2s1_n, r2s1_t, r2s1_pos, r2s1_src = candidates[0]
            if r2s1_pos is not None:
                r2s1_src += f" (~pos {r2s1_pos})"
        use_5p = (r2s1_pos is None or abs(r2s1_pos - eoff) >= args.min_distance)
        if use_5p and r2_5p_seq:
            r2s2, r2s2_n, r2s2_t = r2_5p_seq, r2_5p_n, r2_5p_t
            r2s2_src = f"5' end of R2 hits (~pos {eoff})"
        elif r2_3p_seq:
            r2s2, r2s2_n, r2s2_t = r2_3p_seq, r2_3p_n, r2_3p_t
            r2s2_src = f"3' end of R2 hits (~pos {r2_3p_mpos})"
    else:
        # Too few sgRNA hits in R2 (amplicon may be too long for sgRNA to fall in read).
        # Derive R2 seqs from non-matching R1 reads, which carry correct R2 content.
        if r1_nomatch_5p_seq:
            r2s1     = r1_nomatch_5p_seq
            r2s1_n   = r1_nomatch_5p_n
            r2s1_t   = r1_nomatch_5p_t
            r2s1_pos = eoff
            r2s1_src = f"5' end of R2 reads (~pos {eoff}) [from non-matching R1 reads]"
        if r2_c3p_seq:
            r2s2     = r2_c3p_seq
            r2s2_n   = r2_c3p_n
            r2s2_t   = r2_c3p_t
            r2s2_src = f"3' end of R2 reads (~pos {r2_c3p_mpos}) [from correct R2 reads]"

    # Warn if both R2 seqs are close together
    if r2s1_pos is not None and r2s2 and r2s2_src:
        r2s2_pos_est = r2_c3p_mpos if r2_hits < LOW_R2_HITS else eoff
        if abs(r2s1_pos - r2s2_pos_est) < args.min_distance:
            print(f"\n  WARNING: R2 seq 1 (~pos {r2s1_pos}) and seq 2 (~pos {r2s2_pos_est}) "
                  f"are within {args.min_distance} nt of each other")

    # ---- Output ----
    sep = '=' * 60
    print(f"\n{sep}")
    print("\nR1 results:")
    show("seq 1 — sgRNA-adjacent",  r1s1, r1s1_n, r1s1_t, r1s1_src)
    show("seq 2 — read-end",        r1s2, r1s2_n, r1s2_t, r1s2_src)

    if args.fastq_r2:
        low = r2_hits < LOW_R2_HITS
        print(f"\nR2 results  (sgRNA hits in R2: {r2_hits}"
              f"{', LOW — using fallback derivation from R1' if low else ''}):")
        show("seq 1 — 5' end of read",   r2s1, r2s1_n, r2s1_t, r2s1_src or "")
        show("seq 2 — 3' end of read",   r2s2, r2s2_n, r2s2_t, r2s2_src)
    else:
        print("\nR2: (--fastq-r2 not provided)")

    print(f"\n{sep}")
    print("Suggested config.yaml block:\n")
    print(f"r1_seqs:")
    print(f"  - {r1s1 or '# NOT FOUND'}   # {r1s1_src or ''}")
    print(f"  - {r1s2 or '# NOT FOUND'}   # {r1s2_src or ''}")
    print(f"r2_seqs:")
    print(f"  - {r2s1 or '# NOT FOUND'}   # {r2s1_src or ''}")
    print(f"  - {r2s2 or '# NOT FOUND'}   # {r2s2_src or ''}")


if __name__ == '__main__':
    main()
