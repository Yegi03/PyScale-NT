#!/usr/bin/env python3
"""
PyScale - Lyapunov Exponent Validation Script
Validates the Rosenstein algorithm implementation using the logistic map.
Expected result for r=4.0: λ = ln(2) ≈ 0.693

Authors: Joshua Davies, Yeganeh Abdollahinejad, Darren R. Flower, Amit K. Chattopadhyay
License: MIT
"""

import numpy as np
from math import log


def return_lyapunov(array, periodicity, width,
                    m=2, tau=1, k_max=30, theiler=0, eps=0.0,
                    fit_start=1, fit_end=10,
                    k_neighbors=15, auto_fit=True, r2_min=0.98, min_fit_len=6):

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
        meta = f"method=rosenstein_knn m={m} tau={tau} theiler={theiler} k={k_neighbors} k_max={k_max} eps={eps_cutoff} fit=[{fit_start},{fit_end}] r2={r2:.3f}"
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

        meta = f"method=rosenstein_knn m={m} tau={tau} theiler={theiler} k={k_neighbors} k_max={k_max} eps={eps_cutoff} fit=[{fit_start},{fit_end}] r2={r2:.3f}"
        return (float(slope), meta)


def logistic_map(r=4.0, x0=0.2, n=3000, discard=200):
    x = np.empty(n)
    x[0] = x0
    for i in range(1, n):
        x[i] = r * x[i-1] * (1 - x[i-1])
    return x[discard:]


def sine_wave(n=3000, freq=0.1):
    t = np.arange(n)
    return np.sin(2 * np.pi * freq * t)


if __name__ == "__main__":
    print("=" * 60)
    print("Lyapunov Exponent Verification")
    print("=" * 60)
    
    # Test 1: Chaotic logistic map
    print("\nTest 1: Logistic Map (r=4.0, chaotic)")
    print("-" * 60)
    
    x_chaotic = logistic_map(r=4.0, x0=0.2, n=3000, discard=200)
    expected = log(2)
    print(f"Series: {len(x_chaotic)} points, expected: {expected:.4f}")
    
    val1, meta1 = return_lyapunov(x_chaotic, periodicity=True, width=7, 
                                   k_neighbors=1, auto_fit=False)
    print(f"K=1:  {val1:.4f}  {meta1}")
    
    val15, meta15 = return_lyapunov(x_chaotic, periodicity=True, width=7, 
                                     k_neighbors=15, auto_fit=False)
    print(f"K=15: {val15:.4f}  {meta15}")
    
    val_auto, meta_auto = return_lyapunov(x_chaotic, periodicity=True, width=7, 
                                           k_neighbors=15, auto_fit=True)
    print(f"Auto: {val_auto:.4f}  {meta_auto}")
    
    error = abs(val_auto - expected)
    status = "PASS" if error < expected * 0.15 else "CHECK"
    print(f"Result: {status} (error: {error/expected*100:.1f}%)")
    
    # Test 2: Periodic logistic map
    print("\nTest 2: Logistic Map (r=3.2, periodic)")
    print("-" * 60)
    
    x_periodic = logistic_map(r=3.2, x0=0.2, n=3000, discard=500)
    val_p, meta_p = return_lyapunov(x_periodic, periodicity=True, width=7,
                                     k_neighbors=15, auto_fit=True, r2_min=0.95)
    print(f"Result: {val_p}  {meta_p}")
    status = "PASS" if (np.isnan(val_p) or val_p <= 0.1) else "CHECK"
    print(f"Status: {status}")
    
    # Test 3: Sine wave
    print("\nTest 3: Sine Wave")
    print("-" * 60)
    
    x_sine = sine_wave(n=3000, freq=0.05)
    val_s, meta_s = return_lyapunov(x_sine, periodicity=True, width=7,
                                     k_neighbors=15, auto_fit=True, r2_min=0.95)
    print(f"Result: {val_s}  {meta_s}")
    status = "PASS" if (np.isnan(val_s) or val_s <= 0.05) else "CHECK"
    print(f"Status: {status}")
    
    # Test 4: Stability
    print("\nTest 4: Stability (5 runs)")
    print("-" * 60)
    
    results = []
    for seed in range(5):
        x0 = 0.1 + 0.05 * seed
        x_test = logistic_map(r=4.0, x0=x0, n=3000, discard=200)
        val, _ = return_lyapunov(x_test, periodicity=True, width=7, k_neighbors=15, auto_fit=True)
        results.append(val)
        print(f"  x0={x0:.2f}: {val:.4f}")
    
    print(f"Mean: {np.nanmean(results):.4f} ± {np.nanstd(results):.4f}")
    
    # Test 5: K=1 vs K=15
    print("\nTest 5: K=1 vs K=15")
    print("-" * 60)
    
    results_k1, results_k15 = [], []
    for seed in range(5):
        x0 = 0.1 + 0.05 * seed
        x_test = logistic_map(r=4.0, x0=x0, n=2000, discard=200)
        val1, _ = return_lyapunov(x_test, periodicity=True, width=7, k_neighbors=1, auto_fit=True, r2_min=0.95)
        val15, _ = return_lyapunov(x_test, periodicity=True, width=7, k_neighbors=15, auto_fit=True, r2_min=0.95)
        results_k1.append(val1)
        results_k15.append(val15)
    
    print(f"K=1  std: {np.nanstd(results_k1):.4f}")
    print(f"K=15 std: {np.nanstd(results_k15):.4f}")
    
    print("\n" + "=" * 60)
    print("Done")
    print("=" * 60)
