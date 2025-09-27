#!/usr/bin/env python
"""
Test adaptive T selection
"""

import os
import sys
sys.path.insert(0, '..')

from delta_detector.core import DeltaDetector
from delta_detector.data import load_lambda_from_csv

# Get the correct path to data file
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(SCRIPT_DIR, '..', 'data', 'ramanujan_lambda.csv')


def test_adaptive_on_lehmer():
    """Test adaptive T on Lehmer candidates (very small λ)"""
    
    print("=== Testing adaptive T on Lehmer candidates ===\n")
    
    # Load lambda values
    lambdas = load_lambda_from_csv(DATA_PATH)
    
    # Test primes with different λ magnitudes
    test_cases = [
        (29, 1.162514717),    # Large λ
        (43, -0.017765279),   # Small λ  
        (907, -0.004564446),  # Smallest λ (Lehmer)
        (617, 0.010957076),   # Very small λ
    ]
    
    for p, true_lambda in test_cases:
        print(f"\nTesting p={p} (λ={true_lambda:.6f}):")
        
        # Create detector
        detector = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=10.0)
        detector.set_lambda_values(lambdas)
        
        # Test with fixed T first
        print("\nFixed T:")
        result_fixed = detector.estimate_lambda()
        print(f"  T = {result_fixed.T:.1f}")
        print(f"  |J(T)| = {abs(result_fixed.J_T):.6e}")
        print(f"  λ estimate = {result_fixed.lambda_estimate:.6f}")
        print(f"  Error = {abs(result_fixed.lambda_estimate - true_lambda):.6e}")
        
        # Test with adaptive T
        print("\nAdaptive T:")
        result_adaptive = detector.estimate_lambda_adaptive(
            threshold_ratio=0.01,
            max_iterations=5,
            verbose=True
        )
        print(f"  Final λ estimate = {result_adaptive.lambda_estimate:.6f}")
        print(f"  Final error = {abs(result_adaptive.lambda_estimate - true_lambda):.6e}")
        print(f"  T increased by factor: {result_adaptive.T / result_fixed.T:.1f}")


def test_adaptive_stability():
    """Test that adaptive T gives stable results"""
    
    print("\n\n=== Testing adaptive stability ===\n")
    
    lambdas = load_lambda_from_csv(DATA_PATH)
    
    p = 907  # Lehmer candidate
    
    # Test with different threshold ratios
    thresholds = [0.001, 0.01, 0.1]
    
    print(f"Testing p={p} with different thresholds:")
    results = []
    
    for thresh in thresholds:
        detector = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=10.0)
        detector.set_lambda_values(lambdas)
        
        result = detector.estimate_lambda_adaptive(
            threshold_ratio=thresh,
            verbose=False
        )
        
        results.append(result)
        print(f"  threshold={thresh:.3f}: T={result.T:.1f}, λ={result.lambda_estimate:.6f}")
    
    # Check consistency
    estimates = [r.lambda_estimate for r in results]
    mean_est = sum(estimates) / len(estimates)
    max_dev = max(abs(e - mean_est) for e in estimates)
    
    print(f"\nConsistency check:")
    print(f"  Mean estimate: {mean_est:.6f}")
    print(f"  Max deviation: {max_dev:.6e}")
    
    if max_dev < 1e-4:
        print("  ✓ Estimates are consistent across thresholds")
    else:
        print("  ⚠ Estimates vary significantly with threshold")


def test_adaptive_vs_manual():
    """Compare adaptive with manual T selection"""
    
    print("\n\n=== Comparing adaptive vs manual T selection ===\n")
    
    lambdas = load_lambda_from_csv(DATA_PATH)
    
    # Small λ cases
    small_lambda_primes = [43, 617, 907]
    
    for p in small_lambda_primes:
        true_lambda = lambdas[p]
        print(f"\np={p} (true λ={true_lambda:.6f}):")
        
        # Adaptive
        detector_adaptive = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=10.0)
        detector_adaptive.set_lambda_values(lambdas)
        result_adaptive = detector_adaptive.estimate_lambda_adaptive()
        
        # Manual with large T_mul
        detector_manual = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=100.0)
        detector_manual.set_lambda_values(lambdas)
        result_manual = detector_manual.estimate_lambda()
        
        print(f"  Adaptive: T={result_adaptive.T:.1f}, λ={result_adaptive.lambda_estimate:.6f}")
        print(f"  Manual:   T={result_manual.T:.1f}, λ={result_manual.lambda_estimate:.6f}")
        print(f"  Adaptive error: {abs(result_adaptive.lambda_estimate - true_lambda):.6e}")
        print(f"  Manual error:   {abs(result_manual.lambda_estimate - true_lambda):.6e}")


if __name__ == "__main__":
    test_adaptive_on_lehmer()
    test_adaptive_stability()
    test_adaptive_vs_manual()
    
    print("\n\nAdaptive T test complete!")
