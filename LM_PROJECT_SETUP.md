# LM Project Setup Instructions

## Project Structure

```
E:\engGit\Gem\WN\LM\              # Lemaire Project (Delta Detector)
├── venv/                          # Python virtual environment
├── delta_detector/                # Main package
│   ├── core/                      # Core algorithms
│   │   ├── detector.py           # Main DeltaDetector class
│   │   ├── chebyshev.py         # Chebyshev computations
│   │   ├── weights.py           # Weight functions (Fejér kernel)
│   │   ├── primes.py            # Prime number utilities
│   │   └── params.py            # Default parameters
│   ├── analysis/                  # Analytical tools
│   │   ├── lehmer.py            # Lehmer hypothesis testing
│   │   └── sign_map.py          # Sign map analysis
│   ├── data/                      # Data loading
│   │   ├── loader.py            # Eigenvalue loader
│   │   └── sato_tate.py         # Sato-Tate distribution
│   ├── visualize/                 # Visualization
│   │   └── plots.py             # Plots and diagrams
│   └── cli.py                    # Command line interface
├── data/                          # Eigenvalue data
│   ├── lmfdb/                    # Raw LMFDB data
│   └── processed/                # Processed tables
├── scripts/                       # Demonstration scripts
│   ├── demo_bsplines.py         # B-splines demo
│   ├── demo_improvements.py      # Method improvements demo
│   ├── parse_lmfdb.py           # LMFDB data parser
│   ├── test_detector_real.py     # Tests on real data
│   └── test_extended_data.py     # Extended tests
├── results/                       # Computation results
├── cache/                         # Intermediate computation cache
├── tests/                         # Unit tests
├── setup.py                       # Setup script
├── requirements.txt               # Python dependencies
└── README.md                      # Documentation
```

## Environment Setup

### 1. Open Terminal/PowerShell
```bash
cd E:\engGit\Gem\WN\LM
```

### 2. Create Virtual Environment (if not already created)
```bash
python -m venv venv
```

### 3. Activate Environment
- **Windows (PowerShell)**: `.\venv\Scripts\Activate.ps1`
- **Windows (CMD)**: `.\venv\Scripts\activate.bat`
- **Linux/macOS**: `source venv/bin/activate`

### 4. Install Package in Development Mode
```bash
pip install -e .
```

This will install:
- Core dependencies from requirements.txt:
  - numpy (arrays and computations)
  - scipy (scientific computing) 
  - matplotlib (plotting)
  - click (CLI framework)
- delta_detector package in development mode

### 5. Install Additional Dependencies (optional)
```bash
pip install -r requirements-extras.txt
```

## Usage

### Command Line Interface (CLI)

After installation, the `delta-detect` command is available:

#### 1. Compute Single Eigenvalue
```bash
# Compute λ(p) for prime p
delta-detect prime 29 --sigma 1.3

# With additional parameters
delta-detect prime 29 --sigma 1.3 --T-factor 20 --M 2
```

#### 2. Scan Range of Primes
```bash
# Scan primes from 1000 to 10000
delta-detect scan --start 1000 --end 10000 --output results/scan.csv
```

#### 3. Test Lehmer Hypothesis
```bash
# Test hypothesis up to limit
delta-detect lehmer --limit 1000000 --output results/lehmer_test.csv
```

### Python API

#### Basic Example
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

#### 1. Main Demonstration
```bash
cd scripts
python test_detector_real.py
```

#### 2. Improvements Demo
```bash
python demo_improvements.py
```

#### 3. Parse LMFDB Data
```bash
python parse_lmfdb.py --input ../data/lmfdb/raw_data.txt --output ../data/processed/
```

## Method Parameters

### Core Parameters
- **σ (sigma)**: Real part of s = σ + it (default: 1.3)
- **T**: Window size, T = T_factor × p0^(aσ) (default: T_factor = 20)
- **M**: Maximum power in sum over p^m (default: 2)
- **P**: Prime cutoff (default: 10^6)

### Integration Parameters
- **quad_limit**: Limit of subintervals for quad (default: 1000)
- **quad_epsabs**: Absolute accuracy (default: 1e-12)
- **quad_epsrel**: Relative accuracy (default: 1e-10)

## Data Structure

### Input Data
CSV files with eigenvalues in format:
```csv
n,lambda_n
2,-24
3,252
5,-4830
...
```

### Output Data
- **results/**: Computation results
- **cache/**: Cached intermediate computations
- **logs/**: Execution logs (if enabled)

## Quick Start for New Session

```bash
# 1. Navigate to project directory
cd E:\engGit\Gem\WN\LM

# 2. Activate environment
.\venv\Scripts\Activate.ps1  # PowerShell
# or
.\venv\Scripts\activate.bat   # CMD

# 3. Quick test
delta-detect prime 29 --sigma 1.3

# 4. Run demo
cd scripts
python test_detector_real.py
```

## Installation Verification

```bash
# Check environment is activated
where python
# Should show: E:\engGit\Gem\WN\LM\venv\Scripts\python.exe

# Check installed packages
pip list | findstr "delta"
# Should show: delta-detector ... (editable)

# Check CLI
delta-detect --help
```

## Troubleshooting

1. **"delta-detect: command not found"**
   - Ensure package is installed: `pip install -e .`
   - Check venv activation

2. **"ModuleNotFoundError: No module named 'delta_detector'"**
   - Install package in development mode: `pip install -e .`

3. **numpy/scipy errors**
   - Update dependencies: `pip install -U numpy scipy`

4. **PowerShell execution policy**
   ```powershell
   Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
   ```

## Performance

- Single eigenvalue computation: ~0.1 sec
- Scan 1000 primes: ~2 min
- Memory usage: < 500 MB for standard tasks

## Related Projects

- **RH Project**: E:\engGit\Gem\WN\RH - Riemann Hypothesis computations
- **LMFDB**: https://www.lmfdb.org - Source of eigenvalue data