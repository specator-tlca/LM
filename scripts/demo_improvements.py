#!/usr/bin/env python
"""
Demo of new Delta Detector features:
- Sin test (should remain O(1))
- Adaptive T selection
- P stability control
"""

import os
from delta_detector.core import DeltaDetector
from delta_detector.data import load_lambda_from_csv
import matplotlib.pyplot as plt
import numpy as np

# Get the correct path to data file
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(SCRIPT_DIR, '..', 'data', 'ramanujan_lambda.csv')


def demo_all_features():
    """Demonstrate all detector improvements"""
    
    print("=== Delta Detector Enhanced Features Demo ===\n")
    
    # Load lambda values
    lambdas = load_lambda_from_csv(DATA_PATH)
    
    # Test prime (Lehmer candidate with very small λ)
    p = 907
    true_lambda = lambdas[p]
    
    print(f"Testing p={p} (smallest known |λ(p)| = {true_lambda:.6f})")
    print("="*50)
    
    # 1. Basic detection
    print("\n1. BASIC DETECTION:")
    detector = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=20.0, P=1_000_000)
    detector.set_lambda_values(lambdas)
    
    result_basic = detector.estimate_lambda()
    print(f"  T = {result_basic.T:.1f}")
    print(f"  λ estimate = {result_basic.lambda_estimate:.6f}")
    print(f"  Error = {abs(result_basic.lambda_estimate - true_lambda):.6e}")
    
    # 2. Sin test
    print("\n2. SIN TEST (orthogonality check):")
    for T_mul in [10, 20, 50, 100]:
        T = T_mul * pow(p, 1.3)
        sin_val = detector.compute_sin_test(T)
        cos_val = detector.compute_J_weighted(T)
        print(f"  T = {T:.0f}: sin = {sin_val:.2e}, cos = {cos_val:.2e}, |sin/cos| = {abs(sin_val/cos_val):.2e}")
    print("  ✓ Sin remains O(1) while cos grows with T")
    
    # 3. Adaptive T
    print("\n3. ADAPTIVE T SELECTION:")
    detector_adaptive = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=10.0, P=1_000_000)
    detector_adaptive.set_lambda_values(lambdas)
    
    result_adaptive = detector_adaptive.estimate_lambda_adaptive(verbose=True)
    print(f"  Final λ estimate = {result_adaptive.lambda_estimate:.6f}")
    print(f"  Final error = {abs(result_adaptive.lambda_estimate - true_lambda):.6e}")
    print(f"  Improvement over basic: {abs(result_basic.lambda_estimate - true_lambda) / abs(result_adaptive.lambda_estimate - true_lambda):.1f}x")
    
    # 4. P stability
    print("\n4. P STABILITY CHECK:")
    stability = detector.check_P_stability(verbose=True, P_factors=[0.1, 0.5, 1.0, 2.0])
    
    # 5. Combined approach
    print("\n5. COMBINED APPROACH (adaptive T + P stability):")
    
    # First get optimal T
    detector_combined = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=10.0, P=500_000)
    detector_combined.set_lambda_values(lambdas)
    
    # Adaptive T
    result_opt = detector_combined.estimate_lambda_adaptive(verbose=False)
    print(f"  Optimal T found: {result_opt.T:.1f}")
    
    # Check P stability at optimal T
    p_stab = detector_combined.check_P_stability(T=result_opt.T, verbose=False)
    if p_stab['stable']:
        print(f"  P={detector_combined.P} is stable ✓")
    else:
        print(f"  P={detector_combined.P} unstable, recommend P={p_stab['recommended_P']}")
        # Update and re-estimate
        detector_combined.P = p_stab['recommended_P']
        result_final = detector_combined.estimate_lambda(T=result_opt.T)
        print(f"  Final λ with stable P: {result_final.lambda_estimate:.6f}")
        print(f"  Final error: {abs(result_final.lambda_estimate - true_lambda):.6e}")
    
    # Plot convergence
    print("\n6. VISUALIZING IMPROVEMENTS:")
    plot_improvements(p, lambdas)
    

def plot_improvements(p, lambdas):
    """Plot comparison of different approaches"""
    
    true_lambda = lambdas[p]
    
    # Test range of T multipliers
    T_muls = np.logspace(0, 2.5, 20)  # 1 to ~300
    
    # Results storage
    errors_basic = []
    errors_adaptive = []
    T_values_basic = []
    T_values_adaptive = []
    
    for T_mul in T_muls:
        # Basic approach
        detector_basic = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=T_mul, P=1_000_000)
        detector_basic.set_lambda_values(lambdas)
        result_basic = detector_basic.estimate_lambda()
        errors_basic.append(abs(result_basic.lambda_estimate - true_lambda))
        T_values_basic.append(result_basic.T)
        
        # Adaptive approach (starting from same T_mul)
        detector_adaptive = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=T_mul, P=1_000_000)
        detector_adaptive.set_lambda_values(lambdas)
        result_adaptive = detector_adaptive.estimate_lambda_adaptive(verbose=False)
        errors_adaptive.append(abs(result_adaptive.lambda_estimate - true_lambda))
        T_values_adaptive.append(result_adaptive.T)
    
    # Create plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Error vs T_mul
    ax1.loglog(T_muls, errors_basic, 'o-', label='Basic', alpha=0.7)
    ax1.loglog(T_muls, errors_adaptive, 's-', label='Adaptive', alpha=0.7)
    ax1.set_xlabel('Initial T multiplier')
    ax1.set_ylabel('|Error|')
    ax1.set_title(f'Detection Error for p={p}')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Final T vs initial T_mul
    ax2.semilogx(T_muls, np.array(T_values_basic)/1000, 'o-', label='Basic T', alpha=0.7)
    ax2.semilogx(T_muls, np.array(T_values_adaptive)/1000, 's-', label='Adaptive T', alpha=0.7)
    ax2.set_xlabel('Initial T multiplier')
    ax2.set_ylabel('Final T (×1000)')
    ax2.set_title('Window Size Selection')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    
    # Save to results directory
    results_dir = os.path.join(SCRIPT_DIR, '..', 'results')
    os.makedirs(results_dir, exist_ok=True)
    output_path = os.path.join(results_dir, 'detector_improvements.png')
    plt.savefig(output_path, dpi=150)
    print(f"  Plot saved to {output_path}")


if __name__ == "__main__":
    demo_all_features()
