"""
Generate figures for the Biological Case Study section.
Figure 1: GC Content and LLE comparison between TATA and CpG promoters
Figure 2: Propensity profiles for representative sequences
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Data from analysis
data = {
    'TATA': {
        'names': ['HBB_1', 'TBP_2', 'MYC_1', 'FOS_1', 'PCNA_1'],
        'gc': [45.33, 42.50, 61.00, 65.17, 58.50],
        'twist': [0.0618, 0.0124, 0.0169, 0.0271, 0.0154],
        'gc_scale': [0.0625, 0.0093, 0.0166, 0.0300, 0.0233]
    },
    'CpG': {
        'names': ['GAPDH_1', 'ACTB_1', 'EEF1A1_2', 'RPS6_1', 'UBB_2'],
        'gc': [63.33, 69.83, 53.50, 56.83, 57.17],
        'twist': [0.0287, 0.0388, 0.0252, np.nan, 0.0380],
        'gc_scale': [0.0243, 0.0423, 0.0266, 0.0428, 0.0362]
    }
}

DNA_SCALES = {
    "twist": {"A": 34.4, "T": 34.4, "G": 33.6, "C": 33.6},
    "gc_content": {"A": 0, "T": 0, "G": 1, "C": 1},
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
    """Convert DNA sequence to numerical profile."""
    scale = DNA_SCALES[scale_name]
    return np.array([scale.get(base, 0) for base in sequence.upper()])

def moving_average(profile, window_size=23):
    """Apply moving average smoothing."""
    result = []
    for i in range(len(profile) - window_size + 1):
        result.append(np.mean(profile[i:i + window_size]))
    return np.array(result)


def generate_figure1():
    """Generate comparison bar chart for GC content and LLE."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    colors = {'TATA': '#E74C3C', 'CpG': '#3498DB'}
    
    # Panel A: GC Content
    ax1 = axes[0]
    x_tata = np.arange(len(data['TATA']['names']))
    x_cpg = np.arange(len(data['CpG']['names'])) + len(data['TATA']['names']) + 1
    
    bars1 = ax1.bar(x_tata, data['TATA']['gc'], color=colors['TATA'], alpha=0.8, edgecolor='black', linewidth=0.5)
    bars2 = ax1.bar(x_cpg, data['CpG']['gc'], color=colors['CpG'], alpha=0.8, edgecolor='black', linewidth=0.5)
    
    ax1.axhline(y=np.mean(data['TATA']['gc']), color=colors['TATA'], linestyle='--', linewidth=2, label=f"TATA mean: {np.mean(data['TATA']['gc']):.1f}%")
    ax1.axhline(y=np.mean(data['CpG']['gc']), color=colors['CpG'], linestyle='--', linewidth=2, label=f"CpG mean: {np.mean(data['CpG']['gc']):.1f}%")
    
    ax1.set_xticks(list(x_tata) + list(x_cpg))
    ax1.set_xticklabels(data['TATA']['names'] + data['CpG']['names'], rotation=45, ha='right', fontsize=9)
    ax1.set_ylabel('GC Content (%)', fontsize=11)
    ax1.set_title('(A) GC Content by Promoter Class', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, 80)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Panel B: LLE (twist scale)
    ax2 = axes[1]
    
    tata_twist = [v if not np.isnan(v) else 0 for v in data['TATA']['twist']]
    cpg_twist = [v if not np.isnan(v) else 0 for v in data['CpG']['twist']]
    
    bars3 = ax2.bar(x_tata, tata_twist, color=colors['TATA'], alpha=0.8, edgecolor='black', linewidth=0.5)
    bars4 = ax2.bar(x_cpg, cpg_twist, color=colors['CpG'], alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Mark N/A value
    for i, v in enumerate(data['CpG']['twist']):
        if np.isnan(v):
            ax2.annotate('N/A', (x_cpg[i], 0.002), ha='center', fontsize=8, color='gray')
    
    tata_mean = np.nanmean(data['TATA']['twist'])
    cpg_mean = np.nanmean(data['CpG']['twist'])
    ax2.axhline(y=tata_mean, color=colors['TATA'], linestyle='--', linewidth=2, label=f"TATA mean: {tata_mean:.4f}")
    ax2.axhline(y=cpg_mean, color=colors['CpG'], linestyle='--', linewidth=2, label=f"CpG mean: {cpg_mean:.4f}")
    
    ax2.set_xticks(list(x_tata) + list(x_cpg))
    ax2.set_xticklabels(data['TATA']['names'] + data['CpG']['names'], rotation=45, ha='right', fontsize=9)
    ax2.set_ylabel('Lyapunov Exponent (LLE)', fontsize=11)
    ax2.set_title('(B) LLE Estimates (Twist Scale)', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 0.08)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('../figures/fig_case_study_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('../figures/fig_case_study_comparison.pdf', bbox_inches='tight')
    print("Saved: ../figures/fig_case_study_comparison.png/pdf")
    plt.close()


def generate_figure2():
    """Generate propensity profiles for representative sequences."""
    sequences = read_fasta('../fasta/all_promoters.fa')
    
    # Select representative sequences: TBP_2 (TATA) and ACTB_1 (CpG)
    selected = ['TBP_2', 'ACTB_1']
    colors = {'TBP_2': '#E74C3C', 'ACTB_1': '#3498DB'}
    labels = {'TBP_2': 'TBP (TATA-box)', 'ACTB_1': 'ACTB (CpG island)'}
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    
    tss_position = 500
    window_size = 23
    
    for row, scale_name in enumerate(['twist', 'gc_content']):
        scale_title = 'Twist Angle' if scale_name == 'twist' else 'GC Content'
        
        for col, seq_name in enumerate(selected):
            ax = axes[row, col]
            seq = sequences[seq_name]
            
            profile = sequence_to_profile(seq, scale_name)
            smoothed = moving_average(profile, window_size)
            
            x_raw = np.arange(len(profile))
            x_smooth = np.arange(window_size // 2, len(profile) - window_size // 2)
            
            ax.plot(x_raw, profile, alpha=0.3, color='gray', linewidth=0.5, label='Raw')
            ax.plot(x_smooth, smoothed, color=colors[seq_name], linewidth=1.5, label='Smoothed (w=23)')
            
            ax.axvline(x=tss_position, color='green', linestyle='--', linewidth=2, label='TSS')
            
            if scale_name == 'twist':
                ax.axvspan(tss_position - 40, tss_position - 20, alpha=0.2, color='orange', label='TATA region')
            
            ax.set_xlabel('Position (bp)', fontsize=10)
            ylabel = 'Twist (°)' if scale_name == 'twist' else 'GC (0/1)'
            ax.set_ylabel(ylabel, fontsize=10)
            
            panel_label = chr(65 + row * 2 + col)
            ax.set_title(f'({panel_label}) {labels[seq_name]} - {scale_title}', fontsize=11, fontweight='bold')
            
            if row == 0 and col == 0:
                ax.legend(loc='upper right', fontsize=8)
            
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_xlim(0, 600)
    
    plt.tight_layout()
    plt.savefig('../figures/fig_case_study_profiles.png', dpi=300, bbox_inches='tight')
    plt.savefig('../figures/fig_case_study_profiles.pdf', bbox_inches='tight')
    print("Saved: ../figures/fig_case_study_profiles.png/pdf")
    plt.close()


if __name__ == "__main__":
    print("Generating Case Study Figures...")
    print("=" * 50)
    generate_figure1()
    generate_figure2()
    print("=" * 50)
    print("Done!")
