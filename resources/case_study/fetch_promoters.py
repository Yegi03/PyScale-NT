"""
Fetch human promoter sequences from UCSC Genome Browser for biological case study.
Data source: EPDnew Human (Homo sapiens, hg38)

Authors: Yegi Esfandiari, Joshua Thomas Oludare, Suresh Chattopadhyay
"""

import urllib.request
import json
import os

# Promoter definitions from EPDnew (hg38)
# Format: gene_name: (chromosome, tss_position, strand, has_tata, category)
PROMOTERS = {
    # TATA-box containing promoters (canonical TATA genes)
    "HBB_1":    ("chr11", 5227071, "-", True, "TATA"),
    "TBP_2":    ("chr6", 170553850, "+", True, "TATA"),
    "MYC_1":    ("chr8", 127736231, "+", True, "TATA"),
    "FOS_1":    ("chr14", 75278828, "+", True, "TATA"),
    "PCNA_1":   ("chr20", 5119955, "-", False, "TATA"),  # Note: EPD shows no TATA
    
    # Housekeeping/CpG island promoters (typically TATA-less)
    "GAPDH_1":  ("chr12", 6534517, "+", True, "CpG"),   # Note: EPD shows TATA
    "ACTB_1":   ("chr7", 5530601, "-", True, "CpG"),    # Note: EPD shows TATA
    "EEF1A1_2": ("chr6", 73521585, "-", False, "CpG"),  # Using promoter 2 (no TATA)
    "RPS6_1":   ("chr9", 19380236, "-", False, "CpG"),
    "UBB_2":    ("chr17", 16380779, "+", False, "CpG"),
}

# Region to extract: -500 to +100 relative to TSS (600 bp total)
UPSTREAM = 500
DOWNSTREAM = 100


def fetch_sequence_ucsc(chrom, start, end, strand):
    """Fetch DNA sequence from UCSC Genome Browser API (hg38)."""
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38&chrom={chrom}&start={start}&end={end}"
    
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            data = json.loads(response.read().decode())
            seq = data.get("dna", "")
            
            if strand == "-":
                seq = reverse_complement(seq)
            
            return seq.upper()
    except Exception as e:
        print(f"Error fetching {chrom}:{start}-{end}: {e}")
        return None


def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G",
                  "a": "t", "t": "a", "g": "c", "c": "g", "N": "N", "n": "n"}
    return "".join(complement.get(base, "N") for base in reversed(seq))


def get_promoter_region(chrom, tss, strand, upstream=500, downstream=100):
    """Calculate genomic coordinates for promoter region."""
    if strand == "+":
        start = tss - upstream
        end = tss + downstream
    else:
        start = tss - downstream
        end = tss + upstream
    return start, end


def main():
    output_dir = "fasta"
    os.makedirs(output_dir, exist_ok=True)
    
    all_sequences = []
    tata_sequences = []
    cpg_sequences = []
    
    print("Fetching promoter sequences from UCSC (hg38)...")
    print("=" * 60)
    
    for name, (chrom, tss, strand, has_tata, category) in PROMOTERS.items():
        start, end = get_promoter_region(chrom, tss, strand, UPSTREAM, DOWNSTREAM)
        
        print(f"\n{name}:")
        print(f"  Location: {chrom}:{tss} ({strand})")
        print(f"  Region: {chrom}:{start}-{end}")
        print(f"  Category: {category}, EPD TATA: {'Yes' if has_tata else 'No'}")
        
        seq = fetch_sequence_ucsc(chrom, start, end, strand)
        
        if seq:
            print(f"  Length: {len(seq)} bp")
            print(f"  Sequence (first 60bp): {seq[:60]}...")
            
            gene = name.split("_")[0]
            header = f">{name} | {gene} | {chrom}:{start}-{end}({strand}) | TSS:{tss} | {category}"
            
            fasta_entry = f"{header}\n{seq}\n"
            all_sequences.append(fasta_entry)
            
            if category == "TATA":
                tata_sequences.append(fasta_entry)
            else:
                cpg_sequences.append(fasta_entry)
            
            with open(f"{output_dir}/{name}.fa", "w") as f:
                f.write(fasta_entry)
        else:
            print(f"  ERROR: Could not fetch sequence")
    
    with open(f"{output_dir}/all_promoters.fa", "w") as f:
        f.writelines(all_sequences)
    
    with open(f"{output_dir}/tata_promoters.fa", "w") as f:
        f.writelines(tata_sequences)
    
    with open(f"{output_dir}/cpg_promoters.fa", "w") as f:
        f.writelines(cpg_sequences)
    
    print("\n" + "=" * 60)
    print(f"Created {len(all_sequences)} individual FASTA files")
    print(f"Created combined files:")
    print(f"  - all_promoters.fa ({len(all_sequences)} sequences)")
    print(f"  - tata_promoters.fa ({len(tata_sequences)} sequences)")
    print(f"  - cpg_promoters.fa ({len(cpg_sequences)} sequences)")
    
    print("\n" + "=" * 60)
    print("SUMMARY OF EPD ANNOTATIONS:")
    print("-" * 60)
    print(f"{'Promoter':<12} {'Category':<8} {'EPD TATA':<10} {'Note'}")
    print("-" * 60)
    for name, (chrom, tss, strand, has_tata, category) in PROMOTERS.items():
        expected = (category == "TATA")
        match = "OK" if has_tata == expected else "MISMATCH"
        tata_str = "Yes" if has_tata else "No"
        print(f"{name:<12} {category:<8} {tata_str:<10} {match}")


if __name__ == "__main__":
    main()
