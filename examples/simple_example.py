#!/usr/bin/env python
"""
Simple example of using Delta Detector

This example shows how to:
1. Detect λ(29) 
2. Verify stability
3. Create a sign map
"""

from delta_detector.core import DeltaDetector, get_primes_cached
from delta_detector.data import sample_lambda_ST
import matplotlib.pyplot as plt


def main():
    # Example 1: Detect λ(29)
    print("=== Example 1: Detecting λ(29) ===")
    
    # Create detector
    detector = DeltaDetector(p0=29, sigma=1.3)
    
    # Generate sample lambda values (Sato-Tate)
    # In practice, you would load real values from a file
    primes, _ = get_primes_cached(1_000_000)
    lambda_values = sample_lambda_ST(primes[:10000])  # Use subset for speed
    detector.set_lambda_values(lambda_values)
    
    # Estimate λ(29)
    result = detector.estimate_lambda()
    
    print(f"λ(29) ≈ {result.lambda_estimate:.6f}")
    print(f"Sign: {result.sign}")
    print(f"Confidence: {result.confidence:.1%}")
    print(f"Window size T: {result.T:.1f}")
    print(f"J(T): {result.J_T:.6f}")
    
    # Example 2: Stability analysis
    print("\n=== Example 2: Stability Analysis ===")
    
    stability = detector.analyze_stability()
    print(f"Mean estimate: {stability['mean']:.6f}")
    print(f"Std deviation: {stability['std']:.6f}")
    print(f"Coefficient of variation: {stability['cv']:.1%}")
    print(f"Sign stable across T values: {stability['sign_stable']}")
    
    # Example 3: Quick sign detection for small primes
    print("\n=== Example 3: Signs of λ(p) for p < 100 ===")
    
    small_primes = [p for p in primes if p < 100]
    signs = {}
    
    for p in small_primes[:10]:  # Just first 10 for demo
        det = DeltaDetector(p0=p, sigma=1.3, T_mul=10)  # Smaller T for speed
        det.set_lambda_values(lambda_values)
        res = det.estimate_lambda()
        signs[p] = res.sign
        print(f"λ({p}): {res.sign} (estimate: {res.lambda_estimate:.3f})")
    
    # Example 4: Visualize convergence
    print("\n=== Example 4: Plotting convergence ===")
    
    from delta_detector.visualize import plot_convergence
    plot_convergence(detector, save_to='convergence_p29.png')
    print("Plot saved to convergence_p29.png")
    
    print("\nDone! This was a basic demo of Delta Detector.")
    print("For more examples, see the examples/ directory.")


if __name__ == "__main__":
    main()
