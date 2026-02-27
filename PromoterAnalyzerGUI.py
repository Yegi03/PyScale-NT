#!/usr/bin/env python3
"""
PromoterAnalyzerGUI - A GUI tool for analyzing promoter sequences.
Analyzes TATA motifs, GC content, and Lyapunov complexity of DNA sequences.

Authors: Yegi Esfandiari, Joshua Thomas Oludare, Suresh Chattopadhyay
License: MIT
"""

import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
import numpy as np
import re
import os


class PromoterAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Promoter Sequence Analyzer")
        self.root.geometry("1000x800")
        
        self.sequences = {}
        self.results = []
        
        self.create_widgets()
    
    def create_widgets(self):
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill='both', expand=True, padx=10, pady=10)
        
        self.input_frame = ttk.Frame(notebook)
        notebook.add(self.input_frame, text="Input Sequences")
        self.create_input_tab()
        
        self.analysis_frame = ttk.Frame(notebook)
        notebook.add(self.analysis_frame, text="Analysis Settings")
        self.create_analysis_tab()
        
        self.results_frame = ttk.Frame(notebook)
        notebook.add(self.results_frame, text="Results")
        self.create_results_tab()
    
    def create_input_tab(self):
        load_frame = ttk.LabelFrame(self.input_frame, text="Load Sequences")
        load_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Button(load_frame, text="Load FASTA File", 
                   command=self.load_fasta).pack(side='left', padx=5, pady=5)
        ttk.Button(load_frame, text="Load Directory", 
                   command=self.load_directory).pack(side='left', padx=5, pady=5)
        ttk.Button(load_frame, text="Clear All", 
                   command=self.clear_sequences).pack(side='left', padx=5, pady=5)
        
        self.seq_count_label = ttk.Label(load_frame, text="Sequences loaded: 0")
        self.seq_count_label.pack(side='right', padx=10)
        
        paste_frame = ttk.LabelFrame(self.input_frame, text="Or Paste Sequence (FASTA format)")
        paste_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        self.paste_text = scrolledtext.ScrolledText(paste_frame, height=15)
        self.paste_text.pack(fill='both', expand=True, padx=5, pady=5)
        
        ttk.Button(paste_frame, text="Add Pasted Sequences", 
                   command=self.add_pasted).pack(pady=5)
        
        list_frame = ttk.LabelFrame(self.input_frame, text="Loaded Sequences")
        list_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        self.seq_listbox = tk.Listbox(list_frame, height=8)
        self.seq_listbox.pack(fill='both', expand=True, padx=5, pady=5)
    
    def create_analysis_tab(self):
        tss_frame = ttk.LabelFrame(self.analysis_frame, text="TSS Position")
        tss_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Label(tss_frame, text="TSS position in sequence (0-indexed):").pack(side='left', padx=5)
        self.tss_var = tk.IntVar(value=500)
        ttk.Entry(tss_frame, textvariable=self.tss_var, width=10).pack(side='left', padx=5)
        ttk.Label(tss_frame, text="(Default: 500 for -500 to +100 sequences)").pack(side='left', padx=5)
        
        motif_frame = ttk.LabelFrame(self.analysis_frame, text="TATA Motif Search")
        motif_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Label(motif_frame, text="Search region upstream of TSS:").pack(side='left', padx=5)
        self.upstream_start = tk.IntVar(value=40)
        self.upstream_end = tk.IntVar(value=20)
        ttk.Entry(motif_frame, textvariable=self.upstream_start, width=5).pack(side='left')
        ttk.Label(motif_frame, text="to").pack(side='left', padx=2)
        ttk.Entry(motif_frame, textvariable=self.upstream_end, width=5).pack(side='left')
        ttk.Label(motif_frame, text="bp upstream").pack(side='left', padx=5)
        
        lyap_frame = ttk.LabelFrame(self.analysis_frame, text="Lyapunov Estimation Parameters")
        lyap_frame.pack(fill='x', padx=10, pady=5)
        
        param_grid = ttk.Frame(lyap_frame)
        param_grid.pack(fill='x', padx=5, pady=5)
        
        ttk.Label(param_grid, text="Embedding dim (m):").grid(row=0, column=0, sticky='w', padx=5)
        self.lyap_m = tk.IntVar(value=2)
        ttk.Entry(param_grid, textvariable=self.lyap_m, width=8).grid(row=0, column=1, padx=5)
        
        ttk.Label(param_grid, text="Time delay (τ):").grid(row=0, column=2, sticky='w', padx=5)
        self.lyap_tau = tk.IntVar(value=1)
        ttk.Entry(param_grid, textvariable=self.lyap_tau, width=8).grid(row=0, column=3, padx=5)
        
        ttk.Label(param_grid, text="K neighbors:").grid(row=1, column=0, sticky='w', padx=5)
        self.lyap_k = tk.IntVar(value=15)
        ttk.Entry(param_grid, textvariable=self.lyap_k, width=8).grid(row=1, column=1, padx=5)
        
        ttk.Label(param_grid, text="Window size:").grid(row=1, column=2, sticky='w', padx=5)
        self.window_size = tk.IntVar(value=23)
        ttk.Entry(param_grid, textvariable=self.window_size, width=8).grid(row=1, column=3, padx=5)
        
        self.auto_fit_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(lyap_frame, text="Automatic fit-region selection", 
                        variable=self.auto_fit_var).pack(anchor='w', padx=5)
        
        scale_frame = ttk.LabelFrame(self.analysis_frame, text="Propensity Scales")
        scale_frame.pack(fill='x', padx=10, pady=5)
        
        self.scale_vars = {}
        scales = ['twist', 'gc_content', 'purine', 'rise', 'roll']
        for i, scale in enumerate(scales):
            var = tk.BooleanVar(value=(scale in ['twist', 'gc_content', 'purine']))
            self.scale_vars[scale] = var
            ttk.Checkbutton(scale_frame, text=scale, variable=var).grid(
                row=i//3, column=i%3, sticky='w', padx=10, pady=2)
        
        btn_frame = ttk.Frame(self.analysis_frame)
        btn_frame.pack(fill='x', padx=10, pady=20)
        
        ttk.Button(btn_frame, text="Run Analysis", 
                   command=self.run_analysis, style='Accent.TButton').pack(pady=10)
    
    def create_results_tab(self):
        self.results_text = scrolledtext.ScrolledText(self.results_frame, height=30, font=('Courier', 10))
        self.results_text.pack(fill='both', expand=True, padx=10, pady=10)
        
        btn_frame = ttk.Frame(self.results_frame)
        btn_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Button(btn_frame, text="Export Results", command=self.export_results).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Clear Results", command=self.clear_results).pack(side='left', padx=5)
    
    def load_fasta(self):
        filepath = filedialog.askopenfilename(
            filetypes=[("FASTA files", "*.fa *.fasta *.fna"), ("All files", "*.*")])
        if filepath:
            self.parse_fasta_file(filepath)
    
    def load_directory(self):
        dirpath = filedialog.askdirectory()
        if dirpath:
            for fname in os.listdir(dirpath):
                if fname.endswith(('.fa', '.fasta', '.fna')):
                    self.parse_fasta_file(os.path.join(dirpath, fname))
    
    def parse_fasta_file(self, filepath):
        try:
            current_name = None
            current_seq = []
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_name:
                            self.sequences[current_name] = ''.join(current_seq)
                        current_name = line[1:].split()[0]
                        current_seq = []
                    elif line:
                        current_seq.append(line.upper())
                if current_name:
                    self.sequences[current_name] = ''.join(current_seq)
            self.update_sequence_list()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {e}")
    
    def add_pasted(self):
        text = self.paste_text.get("1.0", tk.END).strip()
        if not text:
            return
        
        lines = text.split('\n')
        current_name = None
        current_seq = []
        
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    self.sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            elif line:
                current_seq.append(line.upper())
        
        if current_name:
            self.sequences[current_name] = ''.join(current_seq)
        
        self.update_sequence_list()
        self.paste_text.delete("1.0", tk.END)
    
    def clear_sequences(self):
        self.sequences = {}
        self.update_sequence_list()
    
    def update_sequence_list(self):
        self.seq_listbox.delete(0, tk.END)
        for name, seq in self.sequences.items():
            self.seq_listbox.insert(tk.END, f"{name} ({len(seq)} bp)")
        self.seq_count_label.config(text=f"Sequences loaded: {len(self.sequences)}")
    
    def run_analysis(self):
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded!")
            return
        
        self.results_text.delete("1.0", tk.END)
        self.results = []
        
        tss = self.tss_var.get()
        upstream_start = self.upstream_start.get()
        upstream_end = self.upstream_end.get()
        
        selected_scales = [s for s, v in self.scale_vars.items() if v.get()]
        
        self.results_text.insert(tk.END, "=" * 80 + "\n")
        self.results_text.insert(tk.END, "PROMOTER SEQUENCE ANALYSIS RESULTS\n")
        self.results_text.insert(tk.END, "=" * 80 + "\n\n")
        self.results_text.insert(tk.END, f"TSS position: {tss}\n")
        self.results_text.insert(tk.END, f"TATA search region: -{upstream_start} to -{upstream_end} bp\n")
        self.results_text.insert(tk.END, f"Scales: {', '.join(selected_scales)}\n\n")
        
        for name, seq in self.sequences.items():
            self.analyze_sequence(name, seq, tss, upstream_start, upstream_end, selected_scales)
        
        self.print_summary()
    
    def analyze_sequence(self, name, seq, tss, up_start, up_end, scales):
        self.results_text.insert(tk.END, "-" * 60 + "\n")
        self.results_text.insert(tk.END, f"Sequence: {name}\n")
        self.results_text.insert(tk.END, f"Length: {len(seq)} bp\n")
        
        gc = self.calculate_gc(seq)
        self.results_text.insert(tk.END, f"GC Content: {gc:.1f}%\n")
        
        if tss <= len(seq):
            region_start = max(0, tss - up_start)
            region_end = tss - up_end
            if region_end > region_start:
                tata_region = seq[region_start:region_end]
                tata_gc = self.calculate_gc(tata_region)
                self.results_text.insert(tk.END, f"TATA region ({-up_start} to {-up_end}): {tata_region}\n")
                self.results_text.insert(tk.END, f"TATA region GC: {tata_gc:.1f}%\n")
                
                motifs = self.find_tata_motifs(seq, region_start, region_end + 10, tss)
                if motifs:
                    self.results_text.insert(tk.END, "TATA motifs found:\n")
                    for m in motifs[:3]:
                        self.results_text.insert(tk.END, 
                            f"  {m['match']} at {m['upstream']}bp upstream\n")
                else:
                    self.results_text.insert(tk.END, "No TATA motif found\n")
        
        lyap_results = {}
        for scale in scales:
            profile = self.sequence_to_profile(seq, scale)
            smoothed = self.moving_average(profile, self.window_size.get())
            lyap_val, _ = self.return_lyapunov(smoothed)
            lyap_results[scale] = lyap_val
            lyap_str = f"{lyap_val:.4f}" if not np.isnan(lyap_val) else "N/A"
            self.results_text.insert(tk.END, f"LLE ({scale}): {lyap_str}\n")
        
        self.results.append({
            'name': name,
            'length': len(seq),
            'gc': gc,
            'has_tata': len(motifs) > 0 if 'motifs' in dir() else False,
            'lyapunov': lyap_results
        })
        
        self.results_text.insert(tk.END, "\n")
    
    def calculate_gc(self, seq):
        if len(seq) == 0:
            return 0
        gc = seq.count('G') + seq.count('C')
        return (gc / len(seq)) * 100
    
    def find_tata_motifs(self, seq, start, end, tss):
        region = seq[start:end]
        patterns = [
            (r'TATAAA', 'canonical'),
            (r'TATAA[AT]', 'variant'),
            (r'[TC]ATA[AT]A', 'variant'),
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
    
    def sequence_to_profile(self, seq, scale_name):
        scales = {
            "twist": {"A": 34.4, "T": 34.4, "G": 33.6, "C": 33.6},
            "rise": {"A": 3.32, "T": 3.32, "G": 3.38, "C": 3.38},
            "roll": {"A": 0.0, "T": 0.0, "G": 4.0, "C": 4.0},
            "gc_content": {"A": 0, "T": 0, "G": 1, "C": 1},
            "purine": {"A": 1, "T": 0, "G": 1, "C": 0},
        }
        scale = scales.get(scale_name, scales["gc_content"])
        return np.array([scale.get(b, 0) for b in seq.upper()])
    
    def moving_average(self, profile, window_size):
        if window_size >= len(profile):
            return np.array([np.mean(profile)])
        result = []
        for i in range(len(profile) - window_size + 1):
            result.append(np.mean(profile[i:i + window_size]))
        return np.array(result)
    
    def return_lyapunov(self, array, m=2, tau=1, k_max=30, k_neighbors=15):
        x = np.asarray(array, dtype=float)
        N = len(x)
        
        if N < 50:
            return (float('nan'), "series too short")
        
        M = N - (m - 1) * tau
        if M <= k_max + 5:
            k_max = max(10, M - 5)
        
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
    
    def print_summary(self):
        self.results_text.insert(tk.END, "=" * 80 + "\n")
        self.results_text.insert(tk.END, "SUMMARY\n")
        self.results_text.insert(tk.END, "=" * 80 + "\n\n")
        
        if not self.results:
            return
        
        gc_values = [r['gc'] for r in self.results]
        self.results_text.insert(tk.END, f"GC Content: mean={np.mean(gc_values):.1f}%, "
                                          f"range={min(gc_values):.1f}-{max(gc_values):.1f}%\n")
        
        tata_count = sum(1 for r in self.results if r.get('has_tata', False))
        self.results_text.insert(tk.END, f"Sequences with TATA motif: {tata_count}/{len(self.results)}\n")
        
        for scale in self.results[0].get('lyapunov', {}).keys():
            lyap_vals = [r['lyapunov'].get(scale, np.nan) for r in self.results]
            valid = [v for v in lyap_vals if not np.isnan(v)]
            if valid:
                self.results_text.insert(tk.END, 
                    f"LLE ({scale}): mean={np.mean(valid):.4f}, std={np.std(valid):.4f}\n")
    
    def export_results(self):
        filepath = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if filepath:
            with open(filepath, 'w') as f:
                f.write(self.results_text.get("1.0", tk.END))
            messagebox.showinfo("Success", f"Results saved to {filepath}")
    
    def clear_results(self):
        self.results_text.delete("1.0", tk.END)
        self.results = []


def main():
    root = tk.Tk()
    app = PromoterAnalyzerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
