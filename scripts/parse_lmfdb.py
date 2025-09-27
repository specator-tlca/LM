#!/usr/bin/env python
"""
Parse LMFDB Ramanujan tau data from full traces file
"""

import json
import csv
from pathlib import Path


def parse_traces_file(filepath):
    """Parse traces.txt file from LMFDB"""
    with open(filepath, 'r') as f:
        # Read all lines
        lines = f.readlines()
    
    # Find the line with the data (skip comments)
    data_line = None
    for line in lines:
        if not line.startswith('#') and line.strip():
            data_line = line.strip()
            break
    
    if not data_line:
        raise ValueError("No data found in file")
    
    # Use json to parse the list
    traces = json.loads(data_line)
    
    # Create dict: index -> value
    traces_dict = {n: tau for n, tau in enumerate(traces)}
    
    return traces_dict


def is_prime(n):
    """Simple primality test"""
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0:
            return False
    return True


def extract_lambda_values(traces_dict, weight=12):
    """Extract normalized λ(p) values from τ(n) traces"""
    lambda_dict = {}
    
    for n, tau_n in traces_dict.items():
        if n > 1 and is_prime(n):
            # Normalize: λ(p) = τ(p) / p^((k-1)/2)
            lambda_p = tau_n / (n ** ((weight - 1) / 2))
            lambda_dict[n] = lambda_p
    
    return lambda_dict


def main():
    # Parse traces
    print("Parsing 1.12.a.a.traces.txt...")
    traces = parse_traces_file("1.12.a.a.traces.txt")
    print(f"Loaded {len(traces)} τ values")
    
    # Extract lambdas
    lambdas = extract_lambda_values(traces)
    print(f"Found {len(lambdas)} prime λ values")
    
    # Show some values
    print("\nFirst few τ values:")
    for n in range(1, 11):
        print(f"  τ({n}) = {traces.get(n, 0)}")
    
    print("\nFirst few λ values:")
    primes = [p for p in sorted(lambdas.keys()) if p < 30]
    for p in primes:
        print(f"  λ({p}) = {lambdas[p]:.10f}")
    
    # Check Ramanujan-Petersson bound
    print("\nChecking |λ(p)| ≤ 2...")
    violations = 0
    max_lambda = 0
    max_p = 0
    
    for p, lam in lambdas.items():
        if abs(lam) > 2:
            violations += 1
        if abs(lam) > max_lambda:
            max_lambda = abs(lam)
            max_p = p
    
    print(f"  Max |λ(p)|: {max_lambda:.6f} at p={max_p}")
    print(f"  Violations: {violations}")
    
    # Save to CSV
    output_file = "ramanujan_lambda.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['p', 'lambda'])
        
        for p in sorted(lambdas.keys()):
            writer.writerow([p, lambdas[p]])
    
    print(f"\nSaved λ values to {output_file}")
    
    # Also save in the format expected by delta_detector
    output_file2 = "data/ramanujan_lambda.csv"
    Path("data").mkdir(exist_ok=True)
    with open(output_file2, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['p', 'lambda'])
        
        for p in sorted(lambdas.keys()):
            writer.writerow([p, lambdas[p]])
    
    print(f"Also saved to {output_file2}")
    
    # Check for small lambda (Lehmer test)
    print("\nSmallest |λ(p)| values (Lehmer test):")
    sorted_by_abs = sorted(lambdas.items(), key=lambda x: abs(x[1]))
    for p, lam in sorted_by_abs[:10]:
        tau_p = traces.get(p, 0)
        print(f"  p={p}: λ(p)={lam:.12f}, τ(p)={tau_p}")


if __name__ == "__main__":
    main()
