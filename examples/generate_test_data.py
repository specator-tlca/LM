#!/usr/bin/env python
"""Generate test data for Delta Detector"""

import csv
from delta_detector.core import primes_up_to
from delta_detector.data import sample_lambda_ST

# Generate primes
primes = primes_up_to(10000)

# Generate Sato-Tate values
lambda_values = sample_lambda_ST(primes, seed=42)

# Save to CSV
with open('test_lambda.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['p', 'lambda'])
    for p in sorted(lambda_values.keys()):
        writer.writerow([p, lambda_values[p]])

print(f"Generated test_lambda.csv with {len(lambda_values)} values")
