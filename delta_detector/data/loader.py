import csv
import json
from typing import Dict, List
from math import pow
import os


def load_lambda_from_csv(filepath: str) -> Dict[int, float]:
    """
    Load lambda values from CSV file
    
    Supported formats:
    1) p,lambda
    2) p,tau  (will convert τ to λ = τ/p^(11/2))
    3) p,value,type  (type = 'lambda' or 'tau')
    """
    lambda_dict: Dict[int, float] = {}
    
    with open(filepath, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        
        # Skip header if present
        first_row = next(reader, None)
        if first_row and not first_row[0].isdigit():
            # It's a header, continue
            pass
        elif first_row:
            # It's data, process it
            reader = [first_row] + list(reader)
        
        for row in reader:
            if not row or row[0].strip().startswith('#'):
                continue
            if len(row) < 2:
                continue
                
            try:
                p = int(row[0].strip())
                val = float(row[1].strip())
                
                # Determine if it's lambda or tau
                if len(row) >= 3:
                    value_type = row[2].strip().lower()
                else:
                    # Heuristic: if |val| <= 2, assume lambda; else tau
                    value_type = 'lambda' if abs(val) <= 2.1 else 'tau'
                
                if value_type == 'lambda':
                    lambda_dict[p] = val
                else:  # tau
                    lambda_dict[p] = val / pow(p, 11.0 / 2.0)
                    
            except (ValueError, IndexError):
                continue
    
    return lambda_dict


def load_lambda_from_json(filepath: str) -> Dict[int, float]:
    """Load lambda values from JSON file"""
    with open(filepath, 'r') as f:
        data = json.load(f)
    
    # Handle different JSON structures
    if isinstance(data, dict):
        # {"2": 0.5, "3": -0.8, ...}
        return {int(k): float(v) for k, v in data.items()}
    elif isinstance(data, list):
        # [{"p": 2, "lambda": 0.5}, ...]
        return {item['p']: item['lambda'] for item in data}
    else:
        raise ValueError("Unknown JSON format")


def save_lambda_to_csv(lambda_dict: Dict[int, float], filepath: str):
    """Save lambda values to CSV file"""
    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['p', 'lambda'])
        for p in sorted(lambda_dict.keys()):
            writer.writerow([p, lambda_dict[p]])


def restrict_to_primes(lambda_dict: Dict[int, float], primes: List[int]) -> Dict[int, float]:
    """Filter lambda dictionary to only include specified primes"""
    prime_set = set(primes)
    return {p: v for p, v in lambda_dict.items() if p in prime_set}
