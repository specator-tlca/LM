#!/usr/bin/env python
"""
Test Delta Detector with extended LMFDB data
After running extract_lmfdb_data.py in sage

This will:
1. Compare with more primes (thousands instead of hundreds)
2. Test higher powers λ(p^2), λ(p^3) directly
3. Search for more Lehmer candidates
"""

import os
import sys
sys.path.insert(0, '..')

from delta_detector.core import DeltaDetector
from delta_detector.data import load_lambda_from_csv
import matplotlib.pyplot as plt
import numpy as np
import csv

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH_ORIGINAL = r"E:\engGit\Gem\WN\LM\data\ramanujan_lambda.csv"
DATA_PATH_FULL = r"E:\engGit\Gem\WN\LM\data\ramanujan_lambda_full.csv"
TAU_PATH_FULL = r"E:\engGit\Gem\WN\LM\data\ramanujan_tau_full.csv"


def load_tau_values(filepath):
    """Load tau(n) values from CSV"""
    tau_dict = {}
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            tau_dict[int(row['n'])] = int(row['tau'])
    return tau_dict


def test_higher_powers():
    """Test detection of λ(p^a) for a > 1"""
    
    print("=== Testing Higher Powers λ(p^a) ===\n")
    
    # Load data
    if not os.path.exists(DATA_PATH_FULL):
        print(f"Extended data not found at {DATA_PATH_FULL}")
        print("Please run 'sage extract_lmfdb_data.py' first")
        return
    
    lambdas = load_lambda_from_csv(DATA_PATH_FULL)
    tau_values = load_tau_values(TAU_PATH_FULL)
    
    # Test cases: prime and power
    test_cases = [
        (2, 2),   # λ(4)
        (2, 3),   # λ(8)
        (3, 2),   # λ(9)
        (5, 2),   # λ(25)
        (7, 2),   # λ(49)
        (11, 2),  # λ(121)
        (13, 2),  # λ(169)
    ]
    
    results = []
    
    for p, a in test_cases:
        print(f"\nTesting λ({p}^{a}) = λ({p**a}):")
        
        # True value from tau
        n = p**a
        if n in tau_values:
            true_lambda = tau_values[n] / (n ** 5.5)
            print(f"  True value: {true_lambda:.12f}")
        else:
            true_lambda = None
            print(f"  True value not available")
            continue
        
        # Detect with our method
        detector = DeltaDetector(p0=p, a=a, sigma=1.3, T_mul=50.0, k=2)
        detector.set_lambda_values(lambdas)
        result = detector.estimate_lambda()
        
        print(f"  Detected:   {result.lambda_estimate:.12f}")
        print(f"  Error:      {abs(result.lambda_estimate - true_lambda):.6e}")
        
        # Verify Chebyshev relation
        if p in lambdas:
            # Compute via recursion
            lambda_p = lambdas[p]
            if a == 2:
                lambda_pa_theory = lambda_p**2 - 1
            elif a == 3:
                lambda_pa_theory = lambda_p**3 - 2*lambda_p
            else:
                # General Chebyshev U_a
                from delta_detector.core.chebyshev import lambda_p_power
                lambda_pa_theory = lambda_p_power(lambda_p, a)
            
            print(f"  Via Chebyshev: {lambda_pa_theory:.12f}")
            print(f"  Chebyshev error: {abs(lambda_pa_theory - true_lambda):.6e}")
        
        results.append({
            'p': p,
            'a': a,
            'true': true_lambda,
            'detected': result.lambda_estimate,
            'error': abs(result.lambda_estimate - true_lambda)
        })
    
    return results


def search_extended_lehmer():
    """Search for Lehmer candidates in extended data"""
    
    print("\n\n=== Extended Lehmer Candidate Search ===\n")
    
    if not os.path.exists(DATA_PATH_FULL):
        print("Extended data not found")
        return
    
    lambdas = load_lambda_from_csv(DATA_PATH_FULL)
    
    # Find smallest |λ(p)|
    prime_lambdas = [(p, lam) for p, lam in lambdas.items()]
    prime_lambdas.sort(key=lambda x: abs(x[1]))
    
    print(f"Found {len(prime_lambdas)} primes total")
    print("\nTop 30 smallest |λ(p)| values:")
    
    lehmer_candidates = []
    
    for i, (p, lam) in enumerate(prime_lambdas[:30]):
        print(f"{i+1:3d}. p={p:5d}: λ = {lam:15.12f}, |λ| = {abs(lam):12.6e}")
        
        # Test with detector
        if i < 10:  # Test first 10 in detail
            detector = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=100.0, k=2)
            detector.set_lambda_values(lambdas)
            result = detector.estimate_lambda_adaptive(verbose=False)
            
            error = abs(result.lambda_estimate - lam)
            print(f"     Detected: {result.lambda_estimate:15.12f}, error: {error:.6e}")
            
            lehmer_candidates.append({
                'p': p,
                'true_lambda': lam,
                'detected_lambda': result.lambda_estimate,
                'error': error,
                'T_used': result.T
            })
    
    # Statistics
    abs_lambdas = [abs(lam) for p, lam in prime_lambdas]
    print(f"\nStatistics of |λ(p)|:")
    print(f"  Minimum: {min(abs_lambdas):.6e}")
    print(f"  1st percentile: {np.percentile(abs_lambdas, 1):.6e}")
    print(f"  5th percentile: {np.percentile(abs_lambdas, 5):.6e}")
    print(f"  Median: {np.median(abs_lambdas):.6f}")
    
    return lehmer_candidates


def plot_extended_distribution():
    """Plot distribution of λ(p) with extended data"""
    
    print("\n\n=== Plotting Extended Distribution ===\n")
    
    if not os.path.exists(DATA_PATH_FULL):
        print("Extended data not found")
        return
    
    lambdas_orig = load_lambda_from_csv(DATA_PATH_ORIGINAL)
    lambdas_full = load_lambda_from_csv(DATA_PATH_FULL)
    
    print(f"Original data: {len(lambdas_orig)} primes")
    print(f"Extended data: {len(lambdas_full)} primes")
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Histogram comparison
    ax = axes[0, 0]
    vals_orig = list(lambdas_orig.values())
    vals_full = list(lambdas_full.values())
    
    bins = np.linspace(-2, 2, 50)
    ax.hist(vals_orig, bins=bins, alpha=0.5, label=f'Original ({len(vals_orig)})', density=True)
    ax.hist(vals_full, bins=bins, alpha=0.5, label=f'Extended ({len(vals_full)})', density=True)
    
    # Sato-Tate density
    x = np.linspace(-2, 2, 1000)
    sato_tate = np.sqrt(1 - (x/2)**2) * 2/np.pi
    ax.plot(x, sato_tate, 'r-', label='Sato-Tate', linewidth=2)
    
    ax.set_xlabel('λ(p)')
    ax.set_ylabel('Density')
    ax.set_title('Distribution of λ(p)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. |λ(p)| vs p
    ax = axes[0, 1]
    primes_full = sorted(lambdas_full.keys())[:5000]  # First 5000
    abs_vals = [abs(lambdas_full[p]) for p in primes_full]
    
    ax.semilogy(primes_full, abs_vals, ',', alpha=0.5)
    ax.set_xlabel('Prime p')
    ax.set_ylabel('|λ(p)|')
    ax.set_title('Magnitude of λ(p)')
    ax.grid(True, alpha=0.3)
    
    # Mark Lehmer candidates
    lehmer_threshold = 0.01
    lehmer_primes = [p for p in primes_full if abs(lambdas_full[p]) < lehmer_threshold]
    lehmer_vals = [abs(lambdas_full[p]) for p in lehmer_primes]
    ax.plot(lehmer_primes, lehmer_vals, 'ro', markersize=5, label=f'|λ| < {lehmer_threshold}')
    ax.legend()
    
    # 3. QQ plot vs Sato-Tate
    ax = axes[1, 0]
    from scipy import stats
    
    # Generate Sato-Tate samples
    def sample_sato_tate(n):
        # Use inverse transform sampling
        u = np.random.uniform(0, 1, n)
        return 2 * np.sin(np.pi * u)
    
    st_samples = sample_sato_tate(len(vals_full))
    
    # QQ plot
    sorted_data = np.sort(vals_full)
    sorted_st = np.sort(st_samples)
    ax.plot(sorted_st, sorted_data, 'b.', alpha=0.5, markersize=1)
    ax.plot([-2, 2], [-2, 2], 'r--', linewidth=2)
    ax.set_xlabel('Sato-Tate quantiles')
    ax.set_ylabel('Data quantiles')
    ax.set_title('QQ Plot vs Sato-Tate')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.2, 2.2)
    
    # 4. Small |λ(p)| zoom
    ax = axes[1, 1]
    small_lambda = [(p, lam) for p, lam in lambdas_full.items() if abs(lam) < 0.1]
    small_lambda.sort(key=lambda x: x[0])
    
    if small_lambda:
        primes_small = [p for p, _ in small_lambda]
        vals_small = [lam for _, lam in small_lambda]
        
        ax.plot(primes_small, vals_small, 'o', markersize=4)
        ax.axhline(0, color='k', linestyle='--', alpha=0.5)
        ax.set_xlabel('Prime p')
        ax.set_ylabel('λ(p)')
        ax.set_title(f'Small λ(p) values (|λ| < 0.1, n={len(small_lambda)})')
        ax.grid(True, alpha=0.3)
        
        # Annotate smallest
        smallest_idx = np.argmin([abs(v) for v in vals_small])
        p_min = primes_small[smallest_idx]
        lam_min = vals_small[smallest_idx]
        ax.annotate(f'p={p_min}\nλ={lam_min:.6f}', 
                   xy=(p_min, lam_min), 
                   xytext=(p_min+1000, lam_min+0.05),
                   arrowprops=dict(arrowstyle='->', color='red'))
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(SCRIPT_DIR, '..', 'results', 'extended_distribution.png')
    plt.savefig(output_path, dpi=150)
    print(f"Plot saved to {output_path}")


def test_extreme_accuracy():
    """Test detector accuracy with B-splines on Lehmer candidates"""
    
    print("\n\n=== Extreme Accuracy Test ===\n")
    
    if not os.path.exists(DATA_PATH_FULL):
        lambdas = load_lambda_from_csv(DATA_PATH_ORIGINAL)
        test_prime = 907  # Known smallest from original data
    else:
        lambdas = load_lambda_from_csv(DATA_PATH_FULL)
        # Find smallest |λ(p)|
        min_p = min(lambdas.keys(), key=lambda p: abs(lambdas[p]))
        test_prime = min_p
    
    true_lambda = lambdas[test_prime]
    
    print(f"Testing p={test_prime} with λ={true_lambda:.15f}")
    print("Using different B-spline orders:\n")
    
    for k in [1, 2, 3, 4]:
        print(f"B-spline order k={k}:")
        
        # Use very large T for maximum accuracy
        detector = DeltaDetector(
            p0=test_prime, 
            a=1, 
            sigma=1.3, 
            T_mul=200.0,  # Very large window
            P=2_000_000,   # Extended prime range
            k=k
        )
        detector.set_lambda_values(lambdas)
        
        # Estimate
        result = detector.estimate_lambda()
        error = abs(result.lambda_estimate - true_lambda)
        
        print(f"  Estimate: {result.lambda_estimate:.15f}")
        print(f"  Error:    {error:.6e}")
        
        # Check if we hit machine precision
        if error < 1e-14:
            print(f"  → Machine precision reached!")
        
        print(f"  T used: {result.T:.1f}")
        print(f"  c(T): {result.c_T:.6f}")
        print()


def main():
    # Check if extended data exists
    if os.path.exists(DATA_PATH_FULL):
        print("Extended LMFDB data found!")
    else:
        print("Extended data not found. Using original data.")
        print(f"To get extended data, run in WSL with sage:")
        print(f"  cd {os.path.dirname(os.path.abspath(__file__))}")
        print(f"  sage extract_lmfdb_data.py")
        print()
    
    # Run tests
    test_higher_powers()
    search_extended_lehmer()
    plot_extended_distribution()
    test_extreme_accuracy()
    
    print("\n\n=== Extended Data Test Complete ===")


if __name__ == "__main__":
    main()
