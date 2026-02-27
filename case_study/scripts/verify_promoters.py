"""
Verify TATA motifs and GC content in downloaded promoter sequences.
TSS is at position 500 (0-indexed), so TATA box should be ~30bp upstream at positions 465-475.

Authors: Yegi Esfandiari, Joshua Thomas Oludare, Suresh Chattopadhyay
"""

import re


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
                header = line[1:]
                current_name = header.split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_name:
            sequences[current_name] = ''.join(current_seq)
    return sequences


def find_tata_motifs(sequence, start=460, end=485):
    """Find TATA-like motifs in a region of the sequence."""
    region = sequence[start:end]
    tata_patterns = [
        (r'TATAAA', 'canonical TATAAA'),
        (r'TATAA[AT]', 'TATAA[A/T]'),
        (r'TATA[AT]A', 'TATA[A/T]A'),
        (r'[TC]ATA[AT]A', '[T/C]ATA[A/T]A'),
        (r'TATA', 'TATA core'),
    ]
    
    found = []
    for pattern, name in tata_patterns:
        for match in re.finditer(pattern, region):
            abs_pos = start + match.start()
            found.append({
                'pattern': name,
                'match': match.group(),
                'position': abs_pos,
                'upstream_of_tss': 500 - abs_pos
            })
    return found


def calculate_gc_content(sequence, start=None, end=None):
    """Calculate GC content of sequence or a region."""
    if start is not None and end is not None:
        seq = sequence[start:end]
    else:
        seq = sequence
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0


def analyze_promoter_region(sequence, tss_pos=500):
    """Analyze the core promoter region around TSS."""
    upstream_100 = sequence[tss_pos-100:tss_pos]
    downstream_50 = sequence[tss_pos:tss_pos+50]
    tata_region = sequence[tss_pos-40:tss_pos-20]
    
    return {
        'upstream_100_gc': calculate_gc_content(upstream_100),
        'downstream_50_gc': calculate_gc_content(downstream_50),
        'tata_region_gc': calculate_gc_content(tata_region),
        'tata_region_seq': tata_region,
        'full_gc': calculate_gc_content(sequence)
    }


def main():
    print("=" * 80)
    print("PROMOTER VERIFICATION: TATA Motifs and GC Content")
    print("=" * 80)
    print("\nSequence structure: 600bp (-500 to +100 relative to TSS)")
    print("TSS position in sequence: 500 (0-indexed)")
    print("Expected TATA box location: ~30bp upstream of TSS (positions 465-475)")
    print()
    
    tata_seqs = read_fasta("../fasta/tata_promoters.fa")
    cpg_seqs = read_fasta("../fasta/cpg_promoters.fa")
    
    print("=" * 80)
    print("TATA-BOX PROMOTERS - Checking for TATA motifs")
    print("=" * 80)
    
    for name, seq in tata_seqs.items():
        print(f"\n--- {name} ---")
        print(f"Sequence length: {len(seq)} bp")
        
        analysis = analyze_promoter_region(seq)
        print(f"Full sequence GC: {analysis['full_gc']:.1f}%")
        print(f"TATA region (-40 to -20): {analysis['tata_region_seq']}")
        print(f"TATA region GC: {analysis['tata_region_gc']:.1f}%")
        
        motifs = find_tata_motifs(seq, start=455, end=485)
        if motifs:
            print("TATA motifs found:")
            for m in motifs:
                print(f"  - {m['match']} ({m['pattern']}) at position {m['position']} "
                      f"({m['upstream_of_tss']}bp upstream of TSS)")
        else:
            print("  No canonical TATA motif found in expected region")
            print(f"  Region 455-485: {seq[455:485]}")
    
    print("\n" + "=" * 80)
    print("CpG ISLAND PROMOTERS - Checking GC content")
    print("=" * 80)
    
    for name, seq in cpg_seqs.items():
        print(f"\n--- {name} ---")
        print(f"Sequence length: {len(seq)} bp")
        
        analysis = analyze_promoter_region(seq)
        print(f"Full sequence GC: {analysis['full_gc']:.1f}%")
        print(f"TATA region (-40 to -20): {analysis['tata_region_seq']}")
        print(f"TATA region GC: {analysis['tata_region_gc']:.1f}%")
        
        motifs = find_tata_motifs(seq, start=455, end=485)
        if motifs:
            print("TATA-like motifs found (unexpected for CpG promoters):")
            for m in motifs:
                print(f"  - {m['match']} ({m['pattern']}) at position {m['position']}")
        else:
            print("  No TATA motif found (expected for CpG island promoters)")
    
    print("\n" + "=" * 80)
    print("SUMMARY COMPARISON")
    print("=" * 80)
    
    tata_gc = [calculate_gc_content(seq) for seq in tata_seqs.values()]
    cpg_gc = [calculate_gc_content(seq) for seq in cpg_seqs.values()]
    
    tata_region_gc_tata = [calculate_gc_content(seq, 460, 480) for seq in tata_seqs.values()]
    tata_region_gc_cpg = [calculate_gc_content(seq, 460, 480) for seq in cpg_seqs.values()]
    
    print(f"\nFull Sequence GC Content:")
    print(f"  TATA promoters: mean={sum(tata_gc)/len(tata_gc):.1f}%, "
          f"range={min(tata_gc):.1f}-{max(tata_gc):.1f}%")
    print(f"  CpG promoters:  mean={sum(cpg_gc)/len(cpg_gc):.1f}%, "
          f"range={min(cpg_gc):.1f}-{max(cpg_gc):.1f}%")
    
    print(f"\nTATA Region (-40 to -20) GC Content:")
    print(f"  TATA promoters: mean={sum(tata_region_gc_tata)/len(tata_region_gc_tata):.1f}%")
    print(f"  CpG promoters:  mean={sum(tata_region_gc_cpg)/len(tata_region_gc_cpg):.1f}%")
    
    print("\n" + "=" * 80)
    print("VERIFICATION STATUS")
    print("=" * 80)
    print("\nExpected patterns:")
    print("  - TATA promoters should have TATAAA or variant ~30bp upstream of TSS")
    print("  - CpG promoters should have higher GC content, especially in core region")
    print("  - CpG promoters typically lack TATA box")


if __name__ == "__main__":
    main()
