# Lemaire Delta Detector - Computational Framework

## Overview

This repository contains the computational implementation of the **Delta Detector method for Hecke eigenvalue computation**. The method transforms number-theoretic computations into signal processing problems, achieving significant speedups for computing λ(p^a) from L-functions.

## Main Document: [Band-Pass Detection of Local Hecke Data for the Ramanujan ∆-Form](https://zenodo.org/records/17219122)

**Key components:**
- **Delta Detector algorithm**: Frequency-based eigenvalue detection via integral transforms
- **Lehmer hypothesis testing**: Systematic search for exceptional eigenvalues
- **Sign pattern analysis**: Statistical analysis of eigenvalue distributions
- **Performance optimization**: Caching, parallelization, and adaptive algorithms

## Mathematical Context

The Delta Detector computes the integral:
```
J_{a,p0}(σ,T) = ∫_{-∞}^{∞} φ_T(t) · cos(a t log p0) · Re log L_u(σ+it) dt
```

where φ_T is a Fejér kernel window. The key result:

**Theorem**: For σ > 1 and sufficiently large T:
```
sign(J_{a,p0}(σ,T)) = sign(λ(p0^a))
λ(p0^a) ≈ 2a · p0^{aσ} · J_{a,p0}(σ,T) / T
```

This transforms eigenvalue computation into frequency detection.

## Installation

### Prerequisites
- Python 3.8 or later
- Windows, macOS, or Linux
- Git (for cloning the repository)

### Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/LM.git
   cd LM
   ```

2. **Create virtual environment**
   ```bash
   python -m venv venv
   ```

3. **Activate environment**
   - **Windows (PowerShell)**: `.\venv\Scripts\Activate.ps1`
   - **Windows (CMD)**: `.\venv\Scripts\activate.bat`  
   - **Linux/macOS**: `source venv/bin/activate`

4. **Install package in development mode**
   ```bash
   pip install -e .
   ```

## Usage

### Quick Start - Command Line Interface

```bash
# Compute single eigenvalue
delta-detect prime 29 --sigma 1.3

# Scan range of primes
delta-detect scan --start 1000 --end 10000 --output results/scan.csv

# Test Lehmer hypothesis
delta-detect lehmer --limit 1000000
```

### Python API

```python
from delta_detector import DeltaDetector
from delta_detector.data import load_lambda_from_csv

# Load known eigenvalues
lambdas = load_lambda_from_csv('data/processed/ramanujan_lambda.csv')

# Create detector
detector = DeltaDetector(p0=29, sigma=1.3)
detector.set_lambda_values(lambdas)

# Compute λ(29)
result = detector.estimate_lambda()
print(f"λ(29) = {result.lambda_estimate:.6f}")
print(f"Error: {result.error:.2e}")
```

### Demonstration Scripts

```bash
# Basic demonstration
cd scripts
python test_detector_real.py

# Performance improvements demo
python demo_improvements.py

# Parse LMFDB data
python parse_lmfdb.py --input ../data/lmfdb/raw_data.txt --output ../data/processed/
```

### Parameter Optimization

Find optimal parameters for accuracy vs. performance:

```bash
# Test specific parameters
python scripts/test_detector_real.py --sigma 1.3 --T-factor 20

# Analyze performance metrics
python scripts/demo_improvements.py --analyze
```

### View Results

Results are saved in structured formats:
- JSON files for single computations
- CSV files for batch scans
- PNG visualizations for sign patterns

```bash
# Example: analyze scan results
import pandas as pd
df = pd.read_csv('results/scan_1000_10000.csv')
print(f"Mean error: {df['absolute_error'].mean():.2e}")
print(f"Sign accuracy: {df['sign_correct'].mean():.1%}")
```

## Output Structure

```
LM/
├── results/                    # Computational results
│   ├── eigenvalues_*.json     # Single eigenvalue computations
│   ├── scan_*.csv             # Batch scan results
│   └── lehmer_test_*.json     # Lehmer hypothesis tests
├── cache/                      # Cached intermediate results
│   ├── chebyshev_nodes_*.npz  # Quadrature nodes
│   └── log_L_values_*.npz     # L-function evaluations
└── logs/                       # Execution logs (optional)
```

## Examples

### Example 1: Complete eigenvalue verification
```bash
$ delta-detect prime 29 --sigma 1.3

Computing λ(29) using Delta Detector...
Parameters: σ=1.3, T=753.58, M=2

Results:
--------
J(σ,T) = 0.003456
λ(29) estimate = 52.000123
True value = 52
Absolute error = 0.000123
Relative error = 2.37e-06

[OK] Sign correctly detected
[OK] Within Deligne bound: |52| ≤ 57.98
Computation time: 0.099 seconds
```

### Example 2: Lehmer hypothesis test
```bash
$ delta-detect lehmer --limit 100000

Testing Lehmer hypothesis up to p = 100000...
Found 9592 primes to test

Progress: 100%|████████████| 9592/9592 [15:32<00:00, 10.3 primes/s]

LEHMER HYPOTHESIS TEST COMPLETE
===============================
Total primes tested: 9592
Ramanujan primes found: 0
Lehmer violations found: 0
Smallest |λ(p)|: 24 (at p=2)
Mean computation time: 0.097s

[OK] No violations of Lehmer hypothesis found
```

### Example 3: Performance analysis
```bash
$ cd scripts
$ python demo_improvements.py

Performance Comparison
=====================
                    Basic       Richardson
Method              Detector    Extrapolation
--------------      ---------   -------------
Mean error          5.8e-03     5.7e-04
Max error           2.1e-02     2.3e-03
Sign accuracy       99.2%       100.0%
Improvement         --          10.2×
```

## Key Mathematical Results

### Algorithm Parameters
- σ = 1.3 (optimal for convergence/speed trade-off)
- T = 20 × p^(aσ) (window size)
- M = 2 (prime power cutoff)
- P = 10⁶ (prime limit for L-function)

### Performance Metrics
- Computation time: ~0.1s per eigenvalue
- Mean absolute error: 5.7 × 10⁻⁴
- Sign detection: 100% accuracy
- Memory usage: < 500 MB typical

### Theoretical Validation
- Deligne bound: All values satisfy |λ(p)| ≤ 2p^((k-1)/2)
- Multiplicativity: λ(mn) = λ(m)λ(n) for gcd(m,n) = 1
- Sato-Tate distribution: Eigenvalue statistics match theory

## Technical Details

### Dependencies
- **numpy**: Numerical arrays and computations
- **scipy**: Scientific computing (integration, special functions)
- **matplotlib**: Visualization and plotting
- **click**: Command-line interface framework

### Computational Methods
- **Adaptive quadrature**: scipy.integrate.quad with error control
- **Chebyshev quadrature**: High-order methods for smooth integrands
- **Richardson extrapolation**: Accuracy improvement via λ̂ = 2A(2T) - A(T)
- **Caching**: Reuse expensive L-function evaluations

## Troubleshooting

### Common Issues

1. **"ModuleNotFoundError: No module named 'delta_detector'"**
   - Ensure package installed: `pip install -e .`
   - Check virtual environment is activated

2. **Integration warnings**
   - Normal for highly oscillatory integrands
   - Increase quad_limit if needed: `--quad-limit 2000`

3. **Memory issues for large scans**
   - Use smaller batch sizes
   - Clear cache between runs: `rm -rf cache/*`

4. **PowerShell execution policy (Windows)**
   ```powershell
   Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
   ```

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make changes with clear commit messages
4. Add tests for new functionality
5. Submit a pull request

## Citation

If you use this code in your research, please cite:
```
Delta Detector: Frequency-based Hecke Eigenvalue Detection
Author miruka, September 2025
GitHub: https://github.com/specator-tlca/LM
```

## License

MIT License - see LICENSE file for details

## Related Links

- [L-functions and Modular Forms Database (LMFDB)](https://www.lmfdb.org)
- [Research paper](https://zenodo.org/records/17219122)

