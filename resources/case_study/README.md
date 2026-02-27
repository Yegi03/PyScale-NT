# Biological Case Study: Human Promoter Analysis

Analysis of TATA-box vs CpG island promoters using PyScale propensity profiles and Lyapunov exponent estimation.

## Data Source

- **Database**: EPDnew (Eukaryotic Promoter Database)
- **Organism**: Homo sapiens (Human)
- **Assembly**: GRCh38/hg38
- **URL**: https://epd.expasy.org/epd/

## Promoters Analyzed

### TATA-box Promoters (n=5)
Genes with canonical TATA boxes in their core promoter region:

| Gene | EPD ID | Chromosome | TSS Position | Strand | EPD TATA |
|------|--------|------------|--------------|--------|----------|
| HBB | HBB_1 | chr11 | 5227071 | - | Yes |
| TBP | TBP_2 | chr6 | 170553850 | + | Yes |
| MYC | MYC_1 | chr8 | 127736231 | + | Yes |
| FOS | FOS_1 | chr14 | 75278828 | + | Yes |
| PCNA | PCNA_1 | chr20 | 5119955 | - | No* |

*Note: PCNA is traditionally classified as TATA-containing, but EPD computational analysis doesn't detect a canonical TATA box.

### CpG Island Promoters (n=5)
Housekeeping genes whose promoters are embedded in CpG islands:

| Gene | EPD ID | Chromosome | TSS Position | Strand | EPD TATA |
|------|--------|------------|--------------|--------|----------|
| GAPDH | GAPDH_1 | chr12 | 6534517 | + | Yes* |
| ACTB | ACTB_1 | chr7 | 5530601 | - | Yes* |
| EEF1A1 | EEF1A1_2 | chr6 | 73521585 | - | No |
| RPS6 | RPS6_1 | chr9 | 19380236 | - | No |
| UBB | UBB_2 | chr17 | 16380779 | + | No |

*Note: GAPDH and ACTB are housekeeping genes but EPD detects TATA-like sequences.

## Sequence Details

- **Region**: -500 to +100 bp relative to TSS (600 bp total)
- **Format**: FASTA

## Files

- `fetch_promoters.py` - Script to download sequences from UCSC Genome Browser API
- `analyze_promoters.py` - Analysis script using PyScale Lyapunov estimator
- `analysis_results.csv` - Results table
- `fasta/` - Directory containing FASTA files:
  - Individual files: `HBB_1.fa`, `TBP_2.fa`, etc.
  - `all_promoters.fa` - All 10 sequences
  - `tata_promoters.fa` - TATA-box promoters only
  - `cpg_promoters.fa` - CpG island promoters only

## Results Summary

### GC Content
- **TATA promoters**: mean 54.5% (range: 42.5-65.2%)
- **CpG promoters**: mean 60.1% (range: 53.5-69.8%)

### Lyapunov Exponent (LLE)

| Scale | TATA (mean ± std) | CpG (mean ± std) | Difference |
|-------|-------------------|------------------|------------|
| twist | 0.027 ± 0.018 | 0.033 ± 0.006 | -0.006 |
| gc_content | 0.028 ± 0.018 | 0.034 ± 0.008 | -0.006 |
| purine | 0.032 ± 0.011 | 0.031 ± 0.014 | +0.001 |

### Key Observations

1. CpG island promoters show higher GC content (expected for CpG islands)
2. CpG promoters have slightly higher and more consistent LLE values with structural scales
3. TATA promoters show more variability (higher std), suggesting more diverse sequence patterns
4. HBB promoter is an outlier with notably high LLE (~0.06)
5. The purine scale shows no significant difference between promoter classes

## Usage

```bash
# Fetch sequences (requires internet)
python fetch_promoters.py

# Run analysis
python analyze_promoters.py
```

## Citation

EPDnew database:
> Meylan P, Dreos R, Ambrosini G, Groux R, Bucher P. EPD in 2020: enhanced data visualization and extension to ncRNA promoters. Nucleic Acids Res. 2020;48(D1):D65-D69.
