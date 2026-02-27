"""
Analyze human promoter sequences using PyScale propensity profiles and Lyapunov estimator.
Biological Case Study: TATA-box vs CpG island promoters

Authors: Yegi Esfandiari, Joshua Thomas Oludare, Suresh Chattopadhyay
"""

import os
import numpy as np


def return_lyapunov(array, periodicity, width,
                    m=2, tau=1, k_max=30, theiler=0, eps=0.0,
                    fit_start=1, fit_end=10,
                    k_neighbors=15, auto_fit=True, r2_min=0.98, min_fit_len=6):
    """
    Estimate the largest Lyapunov exponent using the Rosenstein algorithm
    with K-nearest neighbor averaging and optional automatic fit-region selection.
    """
    x = np.asarray(array, dtype=float)

    if (not periodicity) and (width != 'avg'):
        half_width = int(width) // 2
        if len(x) > 2 * half_width:
            x = x[half_width:-half_width]

    N = len(x)
    if N < 20:
        return (float('nan'), "unavailable: series too short")

    if theiler is None or theiler == 0:
        theiler = m * tau

    eps_cutoff = None if (eps is None or eps <= 0.0) else float(eps)
    k_neighbors = max(1, min(k_neighbors, 30))

    if m < 1 or tau < 1:
        return (float('nan'), "unavailable: m and tau must be >=1")

    M = N - (m - 1) * tau
    if M <= 2:
        return (float('nan'), "unavailable: series too short for embedding")

    if M <= k_max + 1:
        k_max = M - 2
    if k_max < 2:
        return (float('nan'), "unavailable: k_max too large for series length")

    idx = np.arange(M)[:, None] + np.arange(m)[None, :] * tau
    X = x[idx]

    nn = np.full((M, k_neighbors), -1, dtype=int)

    for i in range(M):
        d = np.linalg.norm(X - X[i], axis=1)
        lo = max(0, i - theiler)
        hi = min(M, i + theiler + 1)
        d[lo:hi] = np.inf

        finite_mask = np.isfinite(d)
        if not finite_mask.any():
            continue

        sorted_idxs = np.argsort(d)

        count = 0
        for j in sorted_idxs:
            if count >= k_neighbors:
                break
            if not np.isfinite(d[j]):
                break
            if eps_cutoff is not None and d[j] > eps_cutoff:
                continue
            nn[i, count] = j
            count += 1

    valid = np.where(nn[:, 0] >= 0)[0]
    if len(valid) < max(5, M // 20):
        return (float('nan'), "unavailable: too few neighbor pairs")

    small = 1e-12
    y = np.full(k_max + 1, np.nan)

    for k in range(k_max + 1):
        log_divs = []
        for i in valid:
            neighbors = nn[i, :]
            neighbors = neighbors[neighbors >= 0]

            if len(neighbors) == 0:
                continue

            dists = []
            for j in neighbors:
                if i + k < M and j + k < M:
                    dist = np.linalg.norm(X[i + k] - X[j + k])
                    dists.append(max(dist, small))

            if len(dists) > 0:
                mean_dist = np.mean(dists)
                log_divs.append(np.log(mean_dist))

        if len(log_divs) < 5:
            break
        y[k] = float(np.mean(log_divs))

    if auto_fit:
        best_fit = None
        best_score = -1

        for a in range(k_max + 1):
            for b in range(a + min_fit_len, k_max + 1):
                segment = y[a:b + 1]
                if not np.isfinite(segment).all():
                    continue
                if len(segment) < min_fit_len:
                    continue

                t = np.arange(a, b + 1)
                slope, intercept = np.polyfit(t, segment, 1)

                if slope <= 0:
                    continue

                y_pred = slope * t + intercept
                ss_res = np.sum((segment - y_pred) ** 2)
                ss_tot = np.sum((segment - np.mean(segment)) ** 2)
                if ss_tot < 1e-12:
                    continue
                r2 = 1 - ss_res / ss_tot

                if r2 < r2_min:
                    continue

                seg_len = b - a + 1
                if seg_len > best_score or (seg_len == best_score and best_fit is not None and r2 > best_fit[2]):
                    best_fit = (a, b, r2, slope)
                    best_score = seg_len

        if best_fit is None:
            return (float('nan'), f"unavailable: no linear region (r2_min={r2_min}, min_len={min_fit_len})")

        fit_start, fit_end, r2, slope = best_fit
        meta = f"method=rosenstein_knn m={m} tau={tau} theiler={theiler} k={k_neighbors} k_max={k_max} fit=[{fit_start},{fit_end}] r2={r2:.3f}"
        return (float(slope), meta)

    else:
        if fit_start < 0 or fit_end > k_max or fit_start >= fit_end:
            return (float('nan'), "unavailable: invalid fit range")

        segment = y[fit_start:fit_end + 1]
        if not np.isfinite(segment).all():
            return (float('nan'), "unavailable: insufficient divergence data in fit range")

        t = np.arange(fit_start, fit_end + 1)
        slope, intercept = np.polyfit(t, segment, 1)

        y_pred = slope * t + intercept
        ss_res = np.sum((segment - y_pred) ** 2)
        ss_tot = np.sum((segment - np.mean(segment)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0

        meta = f"method=rosenstein_knn m={m} tau={tau} theiler={theiler} k={k_neighbors} k_max={k_max} fit=[{fit_start},{fit_end}] r2={r2:.3f}"
        return (float(slope), meta)


DNA_SCALES = {
    "twist": {"A": 34.4, "T": 34.4, "G": 33.6, "C": 33.6},
    "rise": {"A": 3.32, "T": 3.32, "G": 3.38, "C": 3.38},
    "roll": {"A": 0.0, "T": 0.0, "G": 4.0, "C": 4.0},
    "tilt": {"A": 0.0, "T": 0.0, "G": 1.0, "C": 1.0},
    "gc_content": {"A": 0, "T": 0, "G": 1, "C": 1},
    "purine": {"A": 1, "T": 0, "G": 1, "C": 0},
    "amino": {"A": 1, "T": 0, "G": 0, "C": 1},
    "weak_strong": {"A": 0, "T": 0, "G": 1, "C": 1},
}


def read_fasta(filepath):
    """Read sequences from FASTA file."""
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_name:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences


def sequence_to_profile(sequence, scale_name):
    """Convert DNA sequence to numerical profile using specified scale."""
    scale = DNA_SCALES[scale_name]
    profile = []
    for base in sequence.upper():
        if base in scale:
            profile.append(scale[base])
        else:
            profile.append(0)
    return np.array(profile)


def moving_average(profile, window_size=23, step=1):
    """Apply moving average smoothing."""
    if window_size >= len(profile):
        return np.array([np.mean(profile)])
    
    result = []
    for i in range(0, len(profile) - window_size + 1, step):
        result.append(np.mean(profile[i:i + window_size]))
    return np.array(result)


def analyze_sequence(name, sequence, scale_name, window_size=23):
    """Analyze a single sequence and return Lyapunov estimate."""
    profile = sequence_to_profile(sequence, scale_name)
    smoothed = moving_average(profile, window_size)
    
    lyap_val, lyap_meta = return_lyapunov(
        array=smoothed,
        periodicity=False,
        width=window_size,
        m=2,
        tau=1,
        k_max=30,
        theiler=0,
        eps=0.0,
        fit_start=1,
        fit_end=10,
        k_neighbors=15,
        auto_fit=True,
        r2_min=0.98,
        min_fit_len=6
    )
    
    return {
        'name': name,
        'length': len(sequence),
        'profile_length': len(smoothed),
        'lyapunov': lyap_val,
        'metadata': lyap_meta,
        'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    }


def main():
    fasta_dir = "fasta"
    
    print("=" * 80)
    print("BIOLOGICAL CASE STUDY: Human Promoter Analysis")
    print("Using PyScale Propensity Profiles and Lyapunov Estimator")
    print("=" * 80)
    
    tata_seqs = read_fasta(f"{fasta_dir}/tata_promoters.fa")
    cpg_seqs = read_fasta(f"{fasta_dir}/cpg_promoters.fa")
    
    print(f"\nLoaded {len(tata_seqs)} TATA-box promoters")
    print(f"Loaded {len(cpg_seqs)} CpG island promoters")
    
    scales_to_test = ["twist", "gc_content", "purine"]
    
    results = {scale: {"TATA": [], "CpG": []} for scale in scales_to_test}
    
    for scale_name in scales_to_test:
        print(f"\n{'='*80}")
        print(f"ANALYSIS WITH '{scale_name.upper()}' SCALE")
        print("="*80)
        
        print("\n--- TATA-box Promoters ---")
        for name, seq in tata_seqs.items():
            result = analyze_sequence(name, seq, scale_name)
            results[scale_name]["TATA"].append(result)
            lyap_str = f"{result['lyapunov']:.4f}" if not np.isnan(result['lyapunov']) else "unavailable"
            print(f"  {name:<12}: LLE={lyap_str:<12} GC={result['gc_content']:.1f}%")
        
        print("\n--- CpG Island Promoters ---")
        for name, seq in cpg_seqs.items():
            result = analyze_sequence(name, seq, scale_name)
            results[scale_name]["CpG"].append(result)
            lyap_str = f"{result['lyapunov']:.4f}" if not np.isnan(result['lyapunov']) else "unavailable"
            print(f"  {name:<12}: LLE={lyap_str:<12} GC={result['gc_content']:.1f}%")
    
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    for scale_name in scales_to_test:
        print(f"\n--- {scale_name.upper()} Scale ---")
        
        tata_lyaps = [r['lyapunov'] for r in results[scale_name]["TATA"] 
                      if not np.isnan(r['lyapunov'])]
        cpg_lyaps = [r['lyapunov'] for r in results[scale_name]["CpG"] 
                     if not np.isnan(r['lyapunov'])]
        
        if tata_lyaps:
            print(f"  TATA promoters (n={len(tata_lyaps)}): "
                  f"mean={np.mean(tata_lyaps):.4f}, std={np.std(tata_lyaps):.4f}")
        else:
            print(f"  TATA promoters: No valid LLE estimates")
        
        if cpg_lyaps:
            print(f"  CpG promoters (n={len(cpg_lyaps)}):  "
                  f"mean={np.mean(cpg_lyaps):.4f}, std={np.std(cpg_lyaps):.4f}")
        else:
            print(f"  CpG promoters: No valid LLE estimates")
        
        if tata_lyaps and cpg_lyaps:
            diff = np.mean(tata_lyaps) - np.mean(cpg_lyaps)
            print(f"  Difference (TATA - CpG): {diff:+.4f}")
    
    print("\n" + "="*80)
    print("GC CONTENT COMPARISON")
    print("="*80)
    
    tata_gc = [r['gc_content'] for r in results[scales_to_test[0]]["TATA"]]
    cpg_gc = [r['gc_content'] for r in results[scales_to_test[0]]["CpG"]]
    
    print(f"  TATA promoters: mean GC = {np.mean(tata_gc):.1f}% (range: {min(tata_gc):.1f}-{max(tata_gc):.1f}%)")
    print(f"  CpG promoters:  mean GC = {np.mean(cpg_gc):.1f}% (range: {min(cpg_gc):.1f}-{max(cpg_gc):.1f}%)")
    
    print("\n" + "="*80)
    print("DETAILED RESULTS TABLE")
    print("="*80)
    print(f"\n{'Promoter':<12} {'Category':<8} {'GC%':<8} ", end="")
    for scale in scales_to_test:
        print(f"{scale[:8]:<12} ", end="")
    print()
    print("-" * 80)
    
    all_results = []
    for category in ["TATA", "CpG"]:
        for i, result in enumerate(results[scales_to_test[0]][category]):
            name = result['name']
            gc = result['gc_content']
            print(f"{name:<12} {category:<8} {gc:<8.1f} ", end="")
            
            row = {'name': name, 'category': category, 'gc': gc}
            for scale in scales_to_test:
                lyap = results[scale][category][i]['lyapunov']
                row[scale] = lyap
                if np.isnan(lyap):
                    print(f"{'N/A':<12} ", end="")
                else:
                    print(f"{lyap:<12.4f} ", end="")
            print()
            all_results.append(row)
    
    with open("analysis_results.csv", "w") as f:
        headers = ["Promoter", "Category", "GC%"] + scales_to_test
        f.write(",".join(headers) + "\n")
        for row in all_results:
            line = [row['name'], row['category'], f"{row['gc']:.2f}"]
            for scale in scales_to_test:
                val = row[scale]
                line.append(f"{val:.4f}" if not np.isnan(val) else "NA")
            f.write(",".join(line) + "\n")
    
    print("\nResults saved to analysis_results.csv")


if __name__ == "__main__":
    main()
