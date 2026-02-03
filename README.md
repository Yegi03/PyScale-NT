# PyScale

A curated nucleotide propensity-scale resource with reproducible profile analysis and Lyapunov-based complexity estimation.

## Overview

PyScale converts nucleotide sequences into numerical propensity profiles using curated scales, enabling alignment-free sequence analysis. The toolkit supports:

- **267 nucleotide propensity scales** (243 mono-, 46 di-, 17 trinucleotide)
- **Four profiling paradigms**: window averaging, periodicity, flip-reading, complementary strand
- **Lyapunov exponent estimation** using an enhanced Rosenstein algorithm with K-nearest neighbor averaging and automatic fit-region selection

## Requirements

- Python 3.6+
- numpy
- matplotlib
- tkinter (for GUI version)

```bash
pip install -r requirements.txt
```

## Directory Structure

```
PyScale/
├── PyScaleTERM.py      # Terminal/command-line version
├── PyScaleGUI.py       # Graphical user interface version
├── README.md           # This file
├── requirements.txt    # Python dependencies
├── test_lyap.py        # Lyapunov estimator validation script
└── resources/
    ├── help.txt        # CLI help documentation
    ├── data_files/     # Input nucleotide sequences (FASTA format)
    └── scales/         # Propensity scale definitions
```

## Usage

### Terminal Version

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

### GUI Version

```bash
python3 PyScaleGUI.py
```

The GUI provides interactive controls for all options and an integrated plotter for comparing multiple sequence profiles.

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

## Input/Output Formats

**Input (FASTA-like):**
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

## Citation

If you use PyScale in your research, please cite:

> Davies J, Abdollahinejad Y, Flower DR, Chattopadhyay AK. PyScale: A Curated Nucleotide Propensity-Scale Resource with Reproducible Profile Analysis and Robust Lyapunov-Based Complexity Estimation. [Journal], [Year].

## References

- Rosenstein MT, Collins JJ, de Luca CJ. (1993). A practical method for calculating largest Lyapunov exponents from small data sets. Physica D 65(1-2): 117-134.

## Acknowledgements

This work has been supported by Aston University and the Nuffield Foundation.

## License

MIT License

Copyright (c) 2024 Joshua Davies, Yeganeh Abdollahinejad, Darren R. Flower, Amit K. Chattopadhyay

