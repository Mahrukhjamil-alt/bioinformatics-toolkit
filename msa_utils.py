"""
Simple center-star Multiple Sequence Alignment.

This is a lightweight heuristic MSA built from repeated pairwise alignments
(no external binaries like Clustal Omega / MUSCLE required). It picks the
sequence with the highest total pairwise similarity as the "center," aligns
every other sequence to it, then merges all the alignments into one set of
MSA columns. Good enough for small demo/teaching alignments — for
publication-quality MSA, use a dedicated tool like Clustal Omega or MAFFT.
"""


def center_star_msa(named_seqs: dict, aligner) -> dict:
    names = list(named_seqs.keys())
    seqs = [named_seqs[n] for n in names]
    n = len(seqs)

    if n < 2:
        return dict(named_seqs)

    totals = [0.0] * n
    for i in range(n):
        for j in range(i + 1, n):
            score = aligner.align(seqs[i], seqs[j]).score
            totals[i] += score
            totals[j] += score
    center = totals.index(max(totals))
    center_seq = seqs[center]
    L = len(center_seq)

    insert_counts = [0] * (L + 1)
    match_residue = {}
    per_seq_inserts = {}

    for i in range(n):
        if i == center:
            continue
        aln = aligner.align(center_seq, seqs[i])[0]
        idx = aln.indices

        matches = {}
        inserts = {}
        cur_center_pos = 0
        j = 0
        ncols = idx.shape[1]
        while j < ncols:
            c, o = idx[0][j], idx[1][j]
            if c == -1:
                run = []
                while j < ncols and idx[0][j] == -1:
                    run.append(seqs[i][idx[1][j]])
                    j += 1
                inserts.setdefault(cur_center_pos, []).append("".join(run))
                insert_counts[cur_center_pos] = max(
                    insert_counts[cur_center_pos], len("".join(run))
                )
            else:
                matches[c] = seqs[i][o] if o != -1 else "-"
                cur_center_pos = c + 1
                j += 1

        match_residue[i] = matches
        per_seq_inserts[i] = inserts

    final = {name: [] for name in names}
    for k in range(L + 1):
        width = insert_counts[k]
        if width > 0:
            for i, name in enumerate(names):
                if i == center:
                    final[name].append("-" * width)
                else:
                    run = "".join(per_seq_inserts.get(i, {}).get(k, []))
                    final[name].append(run.ljust(width, "-"))
        if k < L:
            for i, name in enumerate(names):
                if i == center:
                    final[name].append(center_seq[k])
                else:
                    final[name].append(match_residue.get(i, {}).get(k, "-"))

    return {name: "".join(cols) for name, cols in final.items()}