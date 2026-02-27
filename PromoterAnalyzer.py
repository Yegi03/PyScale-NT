#!/usr/bin/env python3
"""
PromoterAnalyzer - Command-line tool for analyzing promoter sequences.
Analyzes TATA motifs, GC content, and Lyapunov complexity of DNA sequences.

Usage:
    python PromoterAnalyzer.py -f <fasta_file>
    python PromoterAnalyzer.py -d <directory>
    python PromoterAnalyzer.py --interactive

Authors: Yegi Esfandiari, Joshua Thomas Oludare, Suresh Chattopadhyay
License: MIT
"""

import numpy as np
import re
import os
import sys
import argparse


DNA_SCALES = {
    "twist": {"A": 34.4, "T": 34.4, "G": 33.6, "C": 33.6},
    "rise": {"A": 3.32, "T": 3.32, "G": 3.38, "C": 3.38},
    "roll": {"A": 0.0, "T": 0.0, "G": 4.0, "C": 4.0},
    "gc_content": {"A": 0, "T": 0, "G": 1, "C": 1},
    "purine": {"A": 1, "T": 0, "G": 1, "C": 0},
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
            elif line:
                current_seq.append(line.upper())
        if current_name:
            sequences[current_name] = ''.join(current_seq)
    return sequences


def calculate_gc(seq):
    """Calculate GC content percentage."""
    if len(seq) == 0:
        return 0
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100


def find_tata_motifs(seq, start, end, tss):
    """Find TATA-like motifs in a region."""
    region = seq[start:end]
    patterns = [
        (r'TATAAA', 'canonical'),
        (r'TATAA[AT]', 'TATAA[A/T]'),
        (r'[TC]ATA[AT]A', '[T/C]ATA[A/T]A'),
        (r'TATA', 'TATA core'),
    ]
    found = []
    for pattern, ptype in patterns:
        for match in re.finditer(pattern, region):
            abs_pos = start + match.start()
            found.append({
                'match': match.group(),
                'type': ptype,
                'position': abs_pos,
                'upstream': tss - abs_pos
            })
    return found


def sequence_to_profile(seq, scale_name):
    """Convert sequence to numerical profile."""
    scale = DNA_SCALES.get(scale_name, DNA_SCALES["gc_content"])
    return np.array([scale.get(b, 0) for b in seq.upper()])


def moving_average(profile, window_size=23):
    """Apply moving average smoothing."""
    if window_size >= len(profile):
        return np.array([np.mean(profile)])
    result = []
    for i in range(len(profile) - window_size + 1):
        result.append(np.mean(profile[i:i + window_size]))
    return np.array(result)


def return_lyapunov(array, m=2, tau=1, k_max=30, k_neighbors=15, window_size=23):
    """Estimate Lyapunov exponent using Rosenstein algorithm with KNN."""
    x = np.asarray(array, dtype=float)
    
    half_width = window_size // 2
    if len(x) > 2 * half_width:
        x = x[half_width:-half_width]
    
    N = len(x)
    
    if N < 50:
        return (float('nan'), "series too short")
    
    M = N - (m - 1) * tau
    if M <= k_max + 5:
        k_max = max(10, M - 5)
    
    if k_max < 6:
        return (float('nan'), "series too short for analysis")
    
    idx = np.arange(M)[:, None] + np.arange(m)[None, :] * tau
    X = x[idx]
    
    theiler = m * tau
    nn = np.full((M, k_neighbors), -1, dtype=int)
    
    for i in range(M):
        d = np.linalg.norm(X - X[i], axis=1)
        lo, hi = max(0, i - theiler), min(M, i + theiler + 1)
        d[lo:hi] = np.inf
        
        sorted_idx = np.argsort(d)
        count = 0
        for j in sorted_idx:
            if count >= k_neighbors:
                break
            if np.isfinite(d[j]):
                nn[i, count] = j
                count += 1
    
    valid = np.where(nn[:, 0] >= 0)[0]
    if len(valid) < 10:
        return (float('nan'), "too few neighbors")
    
    small = 1e-12
    y = np.full(k_max + 1, np.nan)
    
    for k in range(k_max + 1):
        log_divs = []
        for i in valid:
            neighbors = nn[i, nn[i] >= 0]
            dists = []
            for j in neighbors:
                if i + k < M and j + k < M:
                    dist = np.linalg.norm(X[i + k] - X[j + k])
                    dists.append(max(dist, small))
            if dists:
                log_divs.append(np.log(np.mean(dists)))
        if len(log_divs) >= 5:
            y[k] = np.mean(log_divs)
    
    best_fit = None
    best_score = -1
    min_fit_len = 6
    r2_min = 0.98
    
    for a in range(k_max + 1):
        for b in range(a + min_fit_len, k_max + 1):
            segment = y[a:b + 1]
            if not np.isfinite(segment).all():
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
            
            if r2 >= r2_min:
                seg_len = b - a + 1
                if seg_len > best_score:
                    best_fit = (slope, r2)
                    best_score = seg_len
    
    if best_fit:
        return (best_fit[0], f"r2={best_fit[1]:.3f}")
    return (float('nan'), "no linear region")


def analyze_sequence(name, seq, tss=500, scales=None, window_size=23):
    """Analyze a single sequence."""
    if scales is None:
        scales = ['twist', 'gc_content', 'purine']
    
    result = {
        'name': name,
        'length': len(seq),
        'gc': calculate_gc(seq),
        'tata_motifs': [],
        'lyapunov': {}
    }
    
    if tss <= len(seq):
        region_start = max(0, tss - 45)
        region_end = tss - 15
        result['tata_motifs'] = find_tata_motifs(seq, region_start, region_end + 10, tss)
        result['tata_region'] = seq[tss-40:tss-20] if tss >= 40 else ""
        result['tata_region_gc'] = calculate_gc(result['tata_region'])
    
    for scale in scales:
        profile = sequence_to_profile(seq, scale)
        smoothed = moving_average(profile, window_size)
        lyap_val, meta = return_lyapunov(smoothed, window_size=window_size)
        result['lyapunov'][scale] = {'value': lyap_val, 'meta': meta}
    
    return result


def print_result(result):
    """Print analysis result for one sequence."""
    print(f"\n{'─' * 60}")
    print(f"Sequence: {result['name']}")
    print(f"Length: {result['length']} bp")
    print(f"GC Content: {result['gc']:.1f}%")
    
    if result.get('tata_region'):
        print(f"TATA region (-40 to -20): {result['tata_region']}")
        print(f"TATA region GC: {result['tata_region_gc']:.1f}%")
    
    if result['tata_motifs']:
        print("TATA motifs found:")
        seen = set()
        for m in result['tata_motifs']:
            key = (m['match'], m['position'])
            if key not in seen:
                print(f"  {m['match']} ({m['type']}) at {m['upstream']}bp upstream")
                seen.add(key)
    else:
        print("No TATA motif found")
    
    print("Lyapunov estimates:")
    for scale, data in result['lyapunov'].items():
        val = data['value']
        val_str = f"{val:.4f}" if not np.isnan(val) else "N/A"
        print(f"  {scale}: {val_str}")


def print_summary(results):
    """Print summary statistics."""
    print(f"\n{'═' * 60}")
    print("SUMMARY")
    print('═' * 60)
    
    gc_vals = [r['gc'] for r in results]
    print(f"\nGC Content: mean={np.mean(gc_vals):.1f}%, range={min(gc_vals):.1f}-{max(gc_vals):.1f}%")
    
    tata_count = sum(1 for r in results if r['tata_motifs'])
    print(f"Sequences with TATA: {tata_count}/{len(results)}")
    
    if results and results[0]['lyapunov']:
        for scale in results[0]['lyapunov'].keys():
            vals = [r['lyapunov'][scale]['value'] for r in results]
            valid = [v for v in vals if not np.isnan(v)]
            if valid:
                print(f"LLE ({scale}): mean={np.mean(valid):.4f}, std={np.std(valid):.4f}, n={len(valid)}")


def interactive_mode():
    """Run in interactive mode."""
    print("\n" + "═" * 60)
    print("PROMOTER ANALYZER - Interactive Mode")
    print("═" * 60)
    
    sequences = {}
    
    while True:
        print("\nOptions:")
        print("  1. Load FASTA file")
        print("  2. Load directory of FASTA files")
        print("  3. Enter sequence manually")
        print("  4. Run analysis")
        print("  5. Exit")
        
        choice = input("\nChoice: ").strip()
        
        if choice == '1':
            path = input("Enter FASTA file path: ").strip()
            if os.path.exists(path):
                sequences.update(read_fasta(path))
                print(f"Loaded. Total sequences: {len(sequences)}")
            else:
                print("File not found!")
        
        elif choice == '2':
            path = input("Enter directory path: ").strip()
            if os.path.isdir(path):
                for f in os.listdir(path):
                    if f.endswith(('.fa', '.fasta', '.fna')):
                        sequences.update(read_fasta(os.path.join(path, f)))
                print(f"Loaded. Total sequences: {len(sequences)}")
            else:
                print("Directory not found!")
        
        elif choice == '3':
            name = input("Sequence name: ").strip()
            print("Enter sequence (press Enter twice when done):")
            lines = []
            while True:
                line = input()
                if not line:
                    break
                lines.append(line.upper())
            if lines:
                sequences[name] = ''.join(lines)
                print(f"Added. Total sequences: {len(sequences)}")
        
        elif choice == '4':
            if not sequences:
                print("No sequences loaded!")
                continue
            
            tss = input("TSS position [500]: ").strip()
            tss = int(tss) if tss else 500
            
            print("\n" + "═" * 60)
            print("ANALYSIS RESULTS")
            print("═" * 60)
            
            results = []
            for name, seq in sequences.items():
                result = analyze_sequence(name, seq, tss)
                results.append(result)
                print_result(result)
            
            print_summary(results)
        
        elif choice == '5':
            print("Goodbye!")
            break
        
        else:
            print("Invalid choice!")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze promoter sequences for TATA motifs, GC content, and Lyapunov complexity.')
    parser.add_argument('-f', '--file', help='Input FASTA file')
    parser.add_argument('-d', '--directory', help='Directory containing FASTA files')
    parser.add_argument('--tss', type=int, default=500, help='TSS position in sequence (default: 500)')
    parser.add_argument('--scales', nargs='+', default=['twist', 'gc_content', 'purine'],
                        help='Propensity scales to use')
    parser.add_argument('--window', type=int, default=23, help='Window size for averaging')
    parser.add_argument('-i', '--interactive', action='store_true', help='Run in interactive mode')
    parser.add_argument('-o', '--output', help='Output file (optional)')
    
    args = parser.parse_args()
    
    if args.interactive:
        interactive_mode()
        return
    
    if not args.file and not args.directory:
        parser.print_help()
        print("\nExample usage:")
        print("  python PromoterAnalyzer.py -f sequences.fa")
        print("  python PromoterAnalyzer.py -d fasta_directory/")
        print("  python PromoterAnalyzer.py --interactive")
        return
    
    sequences = {}
    
    if args.file:
        if os.path.exists(args.file):
            sequences.update(read_fasta(args.file))
        else:
            print(f"Error: File not found: {args.file}")
            return
    
    if args.directory:
        if os.path.isdir(args.directory):
            for f in os.listdir(args.directory):
                if f.endswith(('.fa', '.fasta', '.fna')):
                    sequences.update(read_fasta(os.path.join(args.directory, f)))
        else:
            print(f"Error: Directory not found: {args.directory}")
            return
    
    if not sequences:
        print("No sequences loaded!")
        return
    
    print("═" * 60)
    print("PROMOTER SEQUENCE ANALYSIS")
    print("═" * 60)
    print(f"Sequences: {len(sequences)}")
    print(f"TSS position: {args.tss}")
    print(f"Scales: {', '.join(args.scales)}")
    
    results = []
    for name, seq in sequences.items():
        result = analyze_sequence(name, seq, args.tss, args.scales, args.window)
        results.append(result)
        print_result(result)
    
    print_summary(results)
    
    if args.output:
        with open(args.output, 'w') as f:
            f.write("Promoter,Length,GC%,HasTATA")
            for scale in args.scales:
                f.write(f",LLE_{scale}")
            f.write("\n")
            
            for r in results:
                has_tata = "Yes" if r['tata_motifs'] else "No"
                f.write(f"{r['name']},{r['length']},{r['gc']:.2f},{has_tata}")
                for scale in args.scales:
                    val = r['lyapunov'][scale]['value']
                    f.write(f",{val:.4f}" if not np.isnan(val) else ",NA")
                f.write("\n")
        print(f"\nResults saved to: {args.output}")


if __name__ == "__main__":
    main()
