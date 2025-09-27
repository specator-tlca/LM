import numpy as np
from typing import Dict, List
from math import sin, cos, pi
import random


def sample_theta_ST(rng: random.Random) -> float:
    """
    Sample θ from Sato-Tate distribution
    Density: f(θ) = (2/π) sin²(θ) for θ ∈ [0, π]
    """
    while True:
        theta = rng.uniform(0.0, pi)
        u = rng.random()
        if u <= (sin(theta) ** 2):
            return theta


def sample_lambda_ST(primes: List[int], seed: int = 123456) -> Dict[int, float]:
    """
    Generate λ(p) values according to Sato-Tate distribution
    λ(p) = 2 cos(θ) where θ ~ ST distribution
    
    Args:
        primes: List of primes
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary mapping p -> λ(p)
    """
    rng = random.Random(seed)
    result: Dict[int, float] = {}
    
    for p in primes:
        theta = sample_theta_ST(rng)
        result[p] = 2.0 * cos(theta)
    
    return result


def sample_lambda_ST_vectorized(primes: np.ndarray, seed: int = 123456) -> np.ndarray:
    """
    Vectorized version using numpy
    More efficient for large prime lists
    """
    np.random.seed(seed)
    n = len(primes)
    
    # Use rejection sampling
    lambda_values = np.zeros(n)
    remaining = np.ones(n, dtype=bool)
    
    while np.any(remaining):
        # Sample theta uniformly
        theta = np.random.uniform(0, pi, np.sum(remaining))
        # Accept/reject based on sin²(θ)
        u = np.random.random(np.sum(remaining))
        accept = u <= np.sin(theta) ** 2
        
        # Store accepted values
        indices = np.where(remaining)[0]
        accepted_indices = indices[accept]
        lambda_values[accepted_indices] = 2.0 * np.cos(theta[accept])
        remaining[accepted_indices] = False
    
    return lambda_values


def verify_sato_tate_distribution(lambda_values: Dict[int, float], bins: int = 50) -> Dict:
    """
    Verify that lambda values follow Sato-Tate distribution
    Returns statistics for validation
    """
    values = np.array(list(lambda_values.values()))
    
    # Theoretical Sato-Tate density for λ
    # f(λ) = (1/π) * sqrt(1 - (λ/2)²) for |λ| ≤ 2
    
    hist, bin_edges = np.histogram(values, bins=bins, range=(-2, 2), density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Theoretical density
    theoretical = np.zeros_like(bin_centers)
    mask = np.abs(bin_centers) < 2
    theoretical[mask] = (1/pi) * np.sqrt(1 - (bin_centers[mask]/2)**2)
    
    # Compute chi-squared statistic
    expected_counts = theoretical * len(values) * (bin_edges[1] - bin_edges[0])
    observed_counts = hist * len(values) * (bin_edges[1] - bin_edges[0])
    
    # Avoid division by zero
    mask = expected_counts > 5  # Standard chi-squared requirement
    if np.any(mask):
        chi2 = np.sum((observed_counts[mask] - expected_counts[mask])**2 / expected_counts[mask])
    else:
        chi2 = np.inf
    
    return {
        'mean': np.mean(values),
        'std': np.std(values),
        'min': np.min(values),
        'max': np.max(values),
        'chi2': chi2,
        'bins_used': np.sum(mask),
        'theoretical_mean': 0.0,  # ST distribution has mean 0
        'theoretical_std': 1.0     # ST distribution has std 1
    }
