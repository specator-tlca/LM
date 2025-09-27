#!/usr/bin/env sage
"""Generate Hecke eigenvalues using SageMath"""

# This script requires SageMath
# Install: https://www.sagemath.org/

import csv

def generate_ramanujan_tau(limit=10000):
    """Generate Ramanujan tau using SageMath"""
    # Delta function - weight 12 cusp form
    Delta = CuspForms(1, 12).0
    
    # Get coefficients
    coeffs = Delta.qexp(limit).coefficients()
    
    # Extract for primes
    tau_dict = {}
    for p in primes(2, limit):
        tau_dict[p] = coeffs[p]
    
    return tau_dict


def generate_eigenvalues_level_N(N=11, k=2, limit=1000):
    """Generate eigenvalues for newforms of level N, weight k"""
    # Get newforms
    S = CuspForms(N, k)
    newforms = S.newforms()
    
    if not newforms:
        print(f"No newforms for level {N}, weight {k}")
        return None
    
    # Take first newform
    f = newforms[0]
    
    # Get Hecke eigenvalues
    eigenvalues = {}
    for p in primes(2, limit):
        if p.divides(N):
            # Bad prime
            ap = f.hecke_eigenvalue(p)
        else:
            # Good prime
            ap = f.hecke_eigenvalue(p)
        
        # Normalize to get lambda
        lambda_p = ap / p^((k-1)/2)
        eigenvalues[p] = float(lambda_p)
    
    return eigenvalues


def save_to_csv(data, filename, form_type="lambda"):
    """Save eigenvalues to CSV"""
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['p', form_type])
        
        for p in sorted(data.keys()):
            writer.writerow([p, data[p]])
    
    print(f"Saved to {filename}")


# Example usage in SageMath:
if __name__ == "__main__":
    # Ramanujan tau
    print("Generating Ramanujan tau...")
    tau = generate_ramanujan_tau(5000)
    
    # Convert to lambda
    lambda_ramanujan = {p: tau[p]/p^(11/2) for p in tau}
    save_to_csv(lambda_ramanujan, "ramanujan_lambda.csv")
    
    # Level 11, weight 2 form
    print("\nGenerating level 11, weight 2 eigenvalues...")
    lambda_11 = generate_eigenvalues_level_N(11, 2, 5000)
    if lambda_11:
        save_to_csv(lambda_11, "level11_weight2_lambda.csv")
    
    # Some examples
    print("\nFirst few values:")
    for p in primes(2, 30):
        print(f"λ_{{{p}}} = {lambda_ramanujan[p]:.6f}")
