#!/usr/bin/env python
"""Download Hecke eigenvalues from LMFDB"""

import requests
import csv
from math import pow

def download_ramanujan_tau(limit=10000):
    """
    Download Ramanujan tau function coefficients from LMFDB
    Delta form: weight 12, level 1
    """
    url = f"https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/1/12/a/a/download/qexp/{limit}"
    
    print(f"Downloading tau(n) for n up to {limit}...")
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    
    try:
        response = requests.get(url, headers=headers, timeout=30)
        
        if response.status_code != 200:
            print(f"Failed to download (status {response.status_code}). Try manually:")
            print(url)
            return None
    except Exception as e:
        print(f"Download error: {e}")
        print("Try manually or use create_ramanujan_data.py instead")
        return None
    
    # Parse the data
    # Format: q-expansion coefficients, one per line
    lines = response.text.strip().split('\n')
    
    # First line might be header or 0
    tau_values = {}
    for n, line in enumerate(lines):
        if n == 0:
            continue  # Skip n=0
        try:
            tau_values[n] = int(line.strip())
        except:
            continue
    
    return tau_values


def is_prime(n):
    """Simple primality test"""
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True


def tau_to_lambda(tau_values, output_file='ramanujan_lambda.csv'):
    """
    Convert tau(p) to lambda(p) = tau(p)/p^(11/2)
    Save only prime values
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['p', 'tau', 'lambda'])
        
        count = 0
        for p, tau in tau_values.items():
            if is_prime(p):
                lambda_p = tau / pow(p, 11.0/2.0)
                writer.writerow([p, tau, lambda_p])
                count += 1
                
                if count < 10:  # Show first few
                    print(f"p={p}: τ(p)={tau}, λ(p)={lambda_p:.6f}")
        
        print(f"\nSaved {count} prime values to {output_file}")


def download_other_forms():
    """Show how to get other modular forms"""
    
    forms = [
        {
            'name': 'First non-CM form of weight 2',
            'url': 'https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/11/2/a/a/',
            'weight': 2,
            'level': 11
        },
        {
            'name': 'Weight 16 cusp form',
            'url': 'https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/1/16/a/a/',
            'weight': 16,
            'level': 1
        },
        {
            'name': 'Weight 18 cusp form',
            'url': 'https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/1/18/a/a/',
            'weight': 18,
            'level': 1
        }
    ]
    
    print("\nOther modular forms available:")
    for form in forms:
        print(f"\n{form['name']}:")
        print(f"  Weight: {form['weight']}, Level: {form['level']}")
        print(f"  URL: {form['url']}")
        print(f"  Download: {form['url']}download/qexp/10000")


if __name__ == "__main__":
    # Download Ramanujan tau
    tau_values = download_ramanujan_tau(10000)
    
    if tau_values:
        tau_to_lambda(tau_values)
    
    # Show other available forms
    download_other_forms()
