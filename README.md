# PyScale

A curated nucleotide propensity-scale resource with reproducible profile analysis and Lyapunov-based complexity estimation.

## Overview

PyScale converts nucleotide sequences into numerical propensity profiles using curated scales, enabling alignment-free sequence analysis. The toolkit supports:

- **267 nucleotide propensity scales** (243 mono-, 46 di-, 17 trinucleotide)
- **Four profiling paradigms**: window averaging, periodicity, flip-reading, complementary strand
- **Lyapunov exponent estimation** using an enhanced Rosenstein algorithm with K-nearest neighbor averaging and automatic fit-region selection
- **Promoter sequence analysis** with TATA motif detection and GC content calculation

## Requirements

- Python 3.6+
- numpy
- matplotlib
- tkinter (for GUI versions)

```bash
pip install -r requirements.txt
```

## Directory Structure

```
PyScale/
├── PyScaleTERM.py          # Terminal/command-line version
├── PyScaleGUI.py           # Graphical user interface version
├── PromoterAnalyzer.py     # Promoter analysis CLI tool
├── PromoterAnalyzerGUI.py  # Promoter analysis GUI tool
├── README.md               # This file
├── requirements.txt        # Python dependencies
├── test_lyap.py            # Lyapunov estimator validation script
└── resources/
    ├── help.txt            # CLI help documentation
    ├── data_files/         # Input nucleotide sequences (FASTA format)
    ├── scales/             # Propensity scale definitions
    ├── output/             # Generated output files (runtime)
    └── case_study/         # Biological case study (human promoters)
        ├── fasta/          # Promoter FASTA sequences
        ├── analyze_promoters.py
        ├── fetch_promoters.py
        ├── verify_promoters.py
        ├── generate_case_study_figures.py
        └── analysis_results.csv
```

## Usage

### Terminal Version (PyScaleTERM.py)

```bash
python3 PyScaleTERM.py [options]
```

**Basic Options:**

| Flag | Description |
|------|-------------|
| `-h` | Display help menu |
| `-c` | Complementary strand conversion |
| `-r` | Flip-read (reverse sequence) |
| `-p` | Periodicity assumption |
| `-ao` | All output combinations |
| `-f [path]` | Input file location |
| `-o [path]` | Output file location |
| `-s [path]` | Scale file location |
| `-ws [int]` | Window size (odd, default: 7) |
| `-ss [int]` | Step size (even, default: 2) |
| `-wa` | Window averaging mode |

**Lyapunov Options:**

| Flag | Description |
|------|-------------|
| `--lyap_m [int]` | Embedding dimension (default: 2) |
| `--lyap_tau [int]` | Time delay (default: 1) |
| `--lyap_kmax [int]` | Max divergence steps (default: 30) |
| `--lyap_theiler [int]` | Theiler window, 0=auto (default: 0) |
| `--lyap_eps [float]` | Distance cutoff, 0=disabled (default: 0.0) |
| `--lyap_knn [int]` | K-nearest neighbors (default: 15) |
| `--lyap_autofit` | Automatic fit-region selection |
| `--lyap_r2min [float]` | Min R² for auto-fit (default: 0.98) |
| `--lyap_minlen [int]` | Min fit segment length (default: 6) |

### GUI Version (PyScaleGUI.py)

```bash
python3 PyScaleGUI.py
```

The GUI provides interactive controls for all options and an integrated plotter for comparing multiple sequence profiles.

### Promoter Analyzer (PromoterAnalyzer.py)

A specialized tool for analyzing promoter sequences with TATA motif detection and complexity estimation.

```bash
# Analyze a single FASTA file
python3 PromoterAnalyzer.py -f sequence.fa

# Analyze a directory of FASTA files
python3 PromoterAnalyzer.py -d resources/case_study/fasta/

# Interactive mode
python3 PromoterAnalyzer.py -i

# Export results to CSV
python3 PromoterAnalyzer.py -f sequence.fa -o results.csv
```

**Options:**

| Flag | Description |
|------|-------------|
| `-f [path]` | Input FASTA file |
| `-d [path]` | Input directory with FASTA files |
| `--tss [int]` | TSS position in sequence (default: 500) |
| `--scales [list]` | Scales to analyze (default: twist, gc_content, purine) |
| `--window [int]` | Window size for smoothing (default: 23) |
| `-i` | Interactive mode |
| `-o [path]` | Output CSV file |

**Features:**
- TATA motif detection (TATAAA and variants) in the -40 to -20 region upstream of TSS
- GC content calculation (overall and TATA region)
- Lyapunov exponent estimation with multiple propensity scales
- Batch processing of multiple sequences

## Profiling Paradigms

1. **Window averaging** - Smooths profiles using moving window averaging with configurable sizes
2. **Periodicity** - Assumes sequence repeats, wrapping windows around boundaries
3. **Flip-reading** - Reverses nucleotide sequence order for orientation analysis
4. **Complementarity** - Converts to complementary strand (A↔T, C↔G)

## Lyapunov Exponent Estimation

PyScale implements a Rosenstein-style largest Lyapunov exponent (LLE) estimator with enhancements for robustness:

- **K-nearest neighbor averaging** (default K=15): Reduces noise by averaging divergence across multiple neighbors
- **Automatic fit-region selection**: Objectively identifies the optimal linear region based on R² and segment length criteria
- **Theiler window exclusion**: Prevents selection of temporally correlated neighbors
- **Edge trimming**: Removes boundary artifacts from moving average before analysis

### Validation

The estimator has been validated against the chaotic logistic map (r=4.0):
- Theoretical λ = ln(2) ≈ 0.693
- Observed λ = 0.638 ± 0.003 (error < 8%)

For non-chaotic signals, the algorithm correctly returns "unavailable".

### Output Format

```
# lyapunov = 0.6382
# method=rosenstein_knn m=2 tau=1 theiler=2 k=15 k_max=30 eps=None fit=[0,9] r2=0.986
0 0.4567
1 0.5234
...
```

## Testing

Run the Lyapunov estimator validation:

```bash
python3 test_lyap.py
```

## Biological Case Study

The `resources/case_study/` directory contains an analysis of human promoter sequences comparing TATA-box and CpG island promoters.

**Data Source**: EPDnew (Eukaryotic Promoter Database) - Human (hg38)

**Sequences**: 600bp regions (-500 to +100 relative to TSS)

**Promoters analyzed**:
- **TATA-box genes**: HBB, TBP, MYC, FOS, PCNA
- **CpG island genes**: GAPDH, ACTB, EEF1A1, RPS6, UBB

**Key findings**:

| Metric | TATA-box | CpG Island |
|--------|----------|------------|
| Mean GC content | 54.5% | 60.1% |
| Mean LLE (twist) | 0.027 | 0.033 |
| LLE variability | Higher (std 0.018) | Lower (std 0.006) |

- CpG island promoters show higher and more uniform GC content
- CpG promoters have slightly higher and more consistent LLE values
- TATA promoters show more variability, suggesting diverse sequence patterns

**Run the analysis:**

```bash
cd resources/case_study
python3 analyze_promoters.py
python3 generate_case_study_figures.py
```

## Input/Output Formats

**Input (FASTA format):**
```
>SequenceName
ATCGATCGATCG...
```

**Scale file format:**
```
>ScaleName
A 1.23
T 0.45
C 2.34
G 1.56
```

## Authors

- Joshua Davies
- Yeganeh Abdollahinejad
- Darren R. Flower
- Amit K. Chattopadhyay

## Citation

If you use PyScale in your research, please cite:

> Davies J, Abdollahinejad Y, Flower DR, Chattopadhyay AK. PyScale: A Curated Nucleotide Propensity-Scale Resource with Reproducible Profile Analysis and Robust Lyapunov-Based Complexity Estimation. [Journal], [Year].

## References

- Rosenstein MT, Collins JJ, de Luca CJ. (1993). A practical method for calculating largest Lyapunov exponents from small data sets. Physica D 65(1-2): 117-134.
- Meylan P, Dreos R, Ambrosini G, Groux R, Bucher P. (2020). EPD in 2020: enhanced data visualization and extension to ncRNA promoters. Nucleic Acids Research 48(D1): D65-D69.

## License

MIT License

## Repository

https://github.com/Yegi03/PyScale-NT
