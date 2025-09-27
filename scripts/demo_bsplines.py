#!/usr/bin/env python
"""
Demo of B-spline kernels of different orders
Shows how higher order B-splines provide better frequency isolation
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from delta_detector.core import DeltaDetector
from delta_detector.data import load_lambda_from_csv
from delta_detector.core.weights import W_weight_vectorized, Phi_k

# Get the correct path to data file
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(SCRIPT_DIR, '..', 'data', 'ramanujan_lambda.csv')


def visualize_bspline_kernels():
    """Visualize B-spline kernels in frequency domain"""
    
    print("=== B-spline Kernel Comparison ===\n")
    
    # Create frequency range
    xi = np.linspace(-10, 10, 1000)
    
    # Plot different orders
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Linear scale
    for k in [1, 2, 3]:
        phi_vals = [Phi_k(x, k) for x in xi]
        ax1.plot(xi, phi_vals, label=f'k={k}', linewidth=2)
    
    ax1.set_xlabel('ξ')
    ax1.set_ylabel('Φ_k(ξ)')
    ax1.set_title('B-spline Kernels in Frequency Domain (Linear Scale)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xlim(-10, 10)
    
    # Log scale for sidelobes
    xi_pos = xi[xi > 0.1]
    for k in [1, 2, 3]:
        phi_vals = [Phi_k(x, k) for x in xi_pos]
        ax2.loglog(xi_pos, phi_vals, label=f'k={k}', linewidth=2)
    
    # Add reference lines for decay rates
    x_ref = np.logspace(0.5, 1.5, 100)
    ax2.loglog(x_ref, 0.1 * x_ref**(-2), '--', color='gray', alpha=0.5, label='~1/ξ²')
    ax2.loglog(x_ref, 0.01 * x_ref**(-4), '-.', color='gray', alpha=0.5, label='~1/ξ⁴')
    ax2.loglog(x_ref, 0.001 * x_ref**(-6), ':', color='gray', alpha=0.5, label='~1/ξ⁶')
    
    ax2.set_xlabel('ξ')
    ax2.set_ylabel('Φ_k(ξ)')
    ax2.set_title('Sidelobe Decay (Log Scale)')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend()
    ax2.set_xlim(1, 100)
    ax2.set_ylim(1e-10, 1)
    
    plt.tight_layout()
    
    # Save plot
    results_dir = os.path.join(SCRIPT_DIR, '..', 'results')
    os.makedirs(results_dir, exist_ok=True)
    output_path = os.path.join(results_dir, 'bspline_kernels.png')
    plt.savefig(output_path, dpi=150)
    print(f"Kernel comparison saved to {output_path}")
    
    # Print decay rates
    print("\nSidelobe decay rates:")
    print("k=1 (Fejér): ~1/ξ² (-40 dB/decade)")
    print("k=2: ~1/ξ⁴ (-80 dB/decade)")
    print("k=3: ~1/ξ⁶ (-120 dB/decade)")


def test_frequency_isolation():
    """Test how different kernels isolate nearby frequencies"""
    
    print("\n\n=== Frequency Isolation Test ===\n")
    
    # Two close frequencies
    freq1 = 1.0
    freq2 = 1.1  # 10% separation
    T = 100
    
    # Compute weights at both frequencies for different k
    x = np.linspace(0, 2, 1000)
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    for idx, k in enumerate([1, 2, 3]):
        ax = axes[idx]
        
        # Weight centered at freq1
        w1 = W_weight_vectorized(x, freq1, T, k)
        # Weight centered at freq2
        w2 = W_weight_vectorized(x, freq2, T, k)
        
        ax.plot(x, w1, 'b-', label=f'Centered at {freq1}', linewidth=2)
        ax.plot(x, w2, 'r--', label=f'Centered at {freq2}', linewidth=2)
        ax.axvline(freq1, color='b', linestyle=':', alpha=0.5)
        ax.axvline(freq2, color='r', linestyle=':', alpha=0.5)
        
        # Compute cross-talk
        crosstalk1 = W_weight_vectorized(np.array([freq2]), freq1, T, k)[0]
        crosstalk2 = W_weight_vectorized(np.array([freq1]), freq2, T, k)[0]
        
        ax.set_xlabel('Frequency')
        ax.set_ylabel('Weight')
        ax.set_title(f'B-spline k={k}\nCross-talk: {crosstalk1:.2e}')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xlim(0.5, 1.5)
    
    plt.tight_layout()
    
    # Save plot
    results_dir = os.path.join(SCRIPT_DIR, '..', 'results')
    output_path = os.path.join(results_dir, 'frequency_isolation.png')
    plt.savefig(output_path, dpi=150)
    print(f"Frequency isolation plot saved to {output_path}")
    
    # Print cross-talk values
    print("\nCross-talk between frequencies 1.0 and 1.1 (T=100):")
    for k in [1, 2, 3]:
        ct = W_weight_vectorized(np.array([freq2]), freq1, T, k)[0]
        print(f"k={k}: {ct:.6e} ({20*np.log10(ct):.1f} dB)")


def test_detection_with_different_k():
    """Compare detection accuracy with different B-spline orders"""
    
    print("\n\n=== Detection with Different B-spline Orders ===\n")
    
    # Load lambda values
    lambdas = load_lambda_from_csv(DATA_PATH)
    
    # Test primes with different properties
    test_cases = [
        (29, "Large λ"),
        (43, "Small λ"),
        (907, "Smallest λ (Lehmer)"),
    ]
    
    # Test each prime with different k values
    results = {p: {} for p, _ in test_cases}
    
    for p, desc in test_cases:
        print(f"\nTesting p={p} ({desc}, true λ={lambdas[p]:.6f}):")
        
        for k in [1, 2, 3]:
            # Create detector with specific k
            detector = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=20.0, k=k)
            detector.set_lambda_values(lambdas)
            
            # Estimate lambda
            result = detector.estimate_lambda()
            error = abs(result.lambda_estimate - lambdas[p])
            
            results[p][k] = {
                'estimate': result.lambda_estimate,
                'error': error,
                'J_T': result.J_T,
                'c_T': result.c_T
            }
            
            print(f"  k={k}: λ̂={result.lambda_estimate:.6f}, "
                  f"error={error:.2e}, c_T={result.c_T:.3f}")
    
    # Compare results
    print("\n\nSummary of errors by B-spline order:")
    print("Prime    k=1 error    k=2 error    k=3 error")
    print("-" * 45)
    for p, desc in test_cases:
        print(f"p={p:4}: ", end="")
        for k in [1, 2, 3]:
            print(f"{results[p][k]['error']:10.2e} ", end="")
        print()
    
    # Recommendation
    print("\nRecommendations:")
    print("- k=1: Best for general use (narrower main lobe)")
    print("- k=2: Good when dealing with closely-spaced frequencies")
    print("- k=3: Maximum isolation, but wider main lobe reduces resolution")


def test_computational_cost():
    """Measure computational cost of different orders"""
    
    print("\n\n=== Computational Cost Test ===\n")
    
    import time
    
    # Load lambda values
    lambdas = load_lambda_from_csv(DATA_PATH)
    
    p = 29
    n_runs = 5
    
    for k in [1, 2, 3]:
        detector = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=20.0, k=k)
        detector.set_lambda_values(lambdas)
        
        # Time multiple runs
        start = time.time()
        for _ in range(n_runs):
            _ = detector.estimate_lambda()
        elapsed = time.time() - start
        
        avg_time = elapsed / n_runs
        print(f"k={k}: {avg_time*1000:.2f} ms per estimation")
    
    print("\nNote: Higher k requires computing higher powers of sinc,")
    print("but the overhead is usually negligible compared to the sum over primes.")


if __name__ == "__main__":
    # Visualize kernels
    visualize_bspline_kernels()
    
    # Test frequency isolation
    test_frequency_isolation()
    
    # Test detection accuracy
    test_detection_with_different_k()
    
    # Test computational cost
    test_computational_cost()
    
    print("\n\nB-spline demonstration complete!")
    print("Check the results/ directory for visualizations.")
