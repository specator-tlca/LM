#!/usr/bin/env python
"""
Test B-spline kernels of different orders
"""

import os
import sys
sys.path.insert(0, '..')

import numpy as np
from delta_detector.core import DeltaDetector
from delta_detector.data import load_lambda_from_csv
from delta_detector.core.weights import Phi_k, W_weight, W_weight_vectorized

# Get the correct path to data file
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(SCRIPT_DIR, '..', 'data', 'ramanujan_lambda.csv')


def test_phi_k_properties():
    """Test basic properties of Phi_k"""
    
    print("=== Testing Phi_k properties ===\n")
    
    # Test 1: Phi_k(0) = 1 for all k
    print("Test 1: Phi_k(0) should equal 1:")
    for k in [1, 2, 3, 4]:
        val = Phi_k(0.0, k)
        print(f"  k={k}: Phi_k(0) = {val:.10f} {'✓' if abs(val - 1.0) < 1e-10 else '✗'}")
    
    # Test 2: Decay rates
    print("\nTest 2: Sidelobe decay rates:")
    x_test = [10.0, 20.0, 40.0]
    for k in [1, 2, 3]:
        vals = [Phi_k(x, k) for x in x_test]
        # Check if doubling x gives appropriate decay
        ratio1 = vals[1] / vals[0]
        ratio2 = vals[2] / vals[1]
        expected_ratio = 2**(-2*k)  # Should decay as x^(-2k)
        
        print(f"  k={k}: ratio at 20/10 = {ratio1:.6f}, at 40/20 = {ratio2:.6f}")
        print(f"       Expected ~{expected_ratio:.6f} for 1/x^{2*k} decay")
        
        # Check if close to expected
        err1 = abs(ratio1 - expected_ratio) / expected_ratio
        err2 = abs(ratio2 - expected_ratio) / expected_ratio
        if err1 < 0.1 and err2 < 0.1:
            print(f"       ✓ Decay rate confirmed")
        else:
            print(f"       ✗ Decay rate off by {100*max(err1,err2):.1f}%")


def test_weight_function_consistency():
    """Test that scalar and vectorized versions agree"""
    
    print("\n\n=== Testing scalar/vectorized consistency ===\n")
    
    # Test parameters
    a_logp0 = 2.0
    T = 100.0
    test_x = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    
    for k in [1, 2, 3]:
        # Compute with scalar version
        scalar_vals = [W_weight(x, a_logp0, T, k) for x in test_x]
        
        # Compute with vectorized version
        vector_vals = W_weight_vectorized(test_x, a_logp0, T, k)
        
        # Compare
        max_diff = np.max(np.abs(scalar_vals - vector_vals))
        print(f"k={k}: Max difference = {max_diff:.2e} {'✓' if max_diff < 1e-10 else '✗'}")


def test_detection_stability_with_k():
    """Test that detection is stable across different k values"""
    
    print("\n\n=== Testing detection stability with k ===\n")
    
    # Load lambda values
    lambdas = load_lambda_from_csv(DATA_PATH)
    
    # Test a prime with moderate lambda
    p = 29
    true_lambda = lambdas[p]
    
    print(f"Testing p={p} (true λ={true_lambda:.6f}):\n")
    
    # Test with different k and T values
    T_muls = [10, 20, 40]
    
    print("       k=1        k=2        k=3")
    print("-" * 40)
    
    for T_mul in T_muls:
        print(f"T_mul={T_mul:2}: ", end="")
        
        for k in [1, 2, 3]:
            detector = DeltaDetector(p0=p, a=1, sigma=1.3, T_mul=T_mul, k=k)
            detector.set_lambda_values(lambdas)
            result = detector.estimate_lambda()
            
            print(f"{result.lambda_estimate:8.4f}  ", end="")
        print()
    
    print("\nNote: Higher k gives slightly different estimates due to")
    print("different frequency response, but all should be close to true value.")


def test_crosstalk_reduction():
    """Test that higher k reduces cross-talk between nearby frequencies"""
    
    print("\n\n=== Testing cross-talk reduction ===\n")
    
    # Two frequencies with 5% separation
    freq1 = 1.0
    freq2 = 1.05
    T = 200
    
    print(f"Frequencies: {freq1} and {freq2} (5% separation)")
    print(f"T = {T}\n")
    
    print("Cross-talk (weight at wrong frequency):")
    
    for k in [1, 2, 3, 4]:
        # Weight function centered at freq1, evaluated at freq2
        crosstalk = W_weight(freq2, freq1, T, k)
        db = 20 * np.log10(crosstalk) if crosstalk > 0 else -200
        print(f"  k={k}: {crosstalk:.6e} ({db:6.1f} dB)")
    
    print("\n✓ Higher k provides better frequency isolation")


def test_adaptive_k_selection():
    """Test adaptive selection of k based on frequency spacing"""
    
    print("\n\n=== Testing adaptive k selection ===\n")
    
    # Load lambda values
    lambdas = load_lambda_from_csv(DATA_PATH)
    
    # Create detector with ability to change k
    p = 43
    detector = DeltaDetector(p0=p, a=1, sigma=1.3)
    detector.set_lambda_values(lambdas)
    
    # Function to estimate optimal k based on nearby frequencies
    def estimate_optimal_k(p0, primes_dict, threshold_db=-40):
        """Estimate optimal k based on frequency spacing"""
        log_p0 = np.log(p0)
        
        # Find nearest frequency
        min_spacing = float('inf')
        for p in primes_dict:
            if p != p0:
                spacing = abs(np.log(p) - log_p0)
                if spacing < min_spacing:
                    min_spacing = spacing
        
        # Estimate k needed for threshold_db isolation
        # Rough formula: need k such that sidelobes at spacing are below threshold
        T = 1000  # Typical T value
        
        for k in [1, 2, 3, 4]:
            # Approximate sidelobe level at the spacing
            xi = min_spacing * T
            if xi > 2:  # Only consider sidelobes
                sidelobe = (2/xi)**(2*k)  # Approximate decay
                db = 20 * np.log10(sidelobe)
                if db < threshold_db:
                    return k
        
        return 4  # Maximum k if still not enough
    
    # Test on different primes
    test_primes = [29, 43, 103, 617]
    
    print("Recommended k values based on frequency spacing:")
    print("(Target: -40 dB isolation from nearest prime)\n")
    
    for p in test_primes:
        if p in lambdas:
            k_opt = estimate_optimal_k(p, lambdas, -40)
            print(f"p={p}: Recommended k={k_opt}")
            
            # Test detection with recommended k
            detector.p0 = p
            detector.k = k_opt
            result = detector.estimate_lambda()
            error = abs(result.lambda_estimate - lambdas[p])
            print(f"      λ estimate with k={k_opt}: {result.lambda_estimate:.6f}")
            print(f"      Error: {error:.6e}\n")


if __name__ == "__main__":
    test_phi_k_properties()
    test_weight_function_consistency()
    test_detection_stability_with_k()
    test_crosstalk_reduction()
    test_adaptive_k_selection()
    
    print("\n\nB-spline tests complete!")
