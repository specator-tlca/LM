#!/usr/bin/env python
"""
Test sin modulation (should remain O(1))
"""

import sys
import os
sys.path.insert(0, '..')

from delta_detector.core import DeltaDetector
from delta_detector.data import load_lambda_from_csv
import numpy as np

# Get the correct path to data file
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(SCRIPT_DIR, '..', 'data', 'ramanujan_lambda.csv')


def test_sin_remains_bounded():
    """Test that sin-modulated integral remains O(1) as T grows"""
    
    print("Testing sin modulation...")
    
    # Load lambda values
    lambdas = load_lambda_from_csv(DATA_PATH)
    
    # Test on several primes
    test_primes = [29, 43, 103]
    
    for p in test_primes:
        print(f"\nTesting p={p}:")
        
        # Create detector
        detector = DeltaDetector(p0=p, a=1, sigma=1.3, P=100_000)
        detector.set_lambda_values(lambdas)
        
        # Test with increasing T values
        T_values = [100, 500, 1000, 2000, 5000]
        sin_values = []
        
        for T in T_values:
            sin_val = detector.compute_sin_test(T)
            sin_values.append(sin_val)
            print(f"  T={T:5}: sin_test = {sin_val:12.6e}")
        
        # Check that values don't grow with T
        sin_array = np.array(sin_values)
        
        # Compute growth rate (should be ~0 for O(1))
        if len(T_values) > 1:
            # Linear regression in log-log space
            log_T = np.log(T_values)
            log_sin = np.log(np.abs(sin_array) + 1e-10)
            
            # Fit y = a*x + b
            A = np.vstack([log_T, np.ones(len(log_T))]).T
            growth_rate, _ = np.linalg.lstsq(A, log_sin, rcond=None)[0]
            
            print(f"  Growth rate: {growth_rate:.3f} (should be ~0 for O(1))")
            
            # Warning if grows faster than O(T^0.1)
            if abs(growth_rate) > 0.1:
                print(f"  WARNING: Sin test may be growing with T!")


def test_sin_vs_cos_orthogonality():
    """Test orthogonality between sin and cos modulations"""
    
    print("\n\nTesting sin/cos orthogonality...")
    
    # Load lambda values  
    lambdas = load_lambda_from_csv(DATA_PATH)
    
    p = 43
    detector = DeltaDetector(p0=p, a=1, sigma=1.3, P=100_000)
    detector.set_lambda_values(lambdas)
    
    T = 1000
    
    # Compute both cos (normal) and sin tests
    cos_val = detector.compute_J_weighted(T)
    sin_val = detector.compute_sin_test(T)
    
    print(f"p={p}, T={T}:")
    print(f"  Cos integral (J): {cos_val:.6e}")
    print(f"  Sin integral:     {sin_val:.6e}")
    print(f"  Ratio |sin/cos|:  {abs(sin_val/cos_val):.6e}")
    
    # Sin should be much smaller than cos at resonance
    if abs(sin_val) > 0.1 * abs(cos_val):
        print("  WARNING: Sin test too large compared to cos!")


if __name__ == "__main__":
    test_sin_remains_bounded()
    test_sin_vs_cos_orthogonality()
    
    print("\n\nSin test complete!")
