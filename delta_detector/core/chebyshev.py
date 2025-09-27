from typing import Dict, Optional
import numpy as np

def lambda_p_power(lambda_p: float, m: int) -> float:
    """
    Compute λ(p^m) using Hecke/Chebyshev recursion
    λ(p^0) = 1
    λ(p^1) = λ(p)
    λ(p^{m+1}) = λ(p) * λ(p^m) - λ(p^{m-1})
    """
    if m == 0:
        return 1.0
    if m == 1:
        return float(lambda_p)
    
    u_prev2 = 1.0
    u_prev1 = float(lambda_p)
    
    for _ in range(2, m + 1):
        u = lambda_p * u_prev1 - u_prev2
        u_prev2, u_prev1 = u_prev1, u
    
    return u_prev1


class LambdaPMCache:
    """Cache for λ(p^m) values to avoid recomputation"""
    
    def __init__(self):
        self.cache: Dict[int, Dict[int, float]] = {}
    
    def get(self, p: int, m: int, lambda_p: float) -> float:
        """Get λ(p^m) with caching"""
        d = self.cache.setdefault(p, {0: 1.0, 1: float(lambda_p)})
        
        if m in d:
            return d[m]
        
        # Fill cache up to m
        k = max(d.keys())
        u_prev2 = d[k - 1] if k >= 1 else 1.0
        u_prev1 = d[k]
        
        for r in range(k + 1, m + 1):
            u = lambda_p * u_prev1 - u_prev2
            d[r] = u
            u_prev2, u_prev1 = u_prev1, u
        
        return d[m]
    
    def clear(self):
        """Clear the cache"""
        self.cache.clear()


def lambda_p_power_vectorized(lambda_p: np.ndarray, m: int) -> np.ndarray:
    """Vectorized computation of λ(p^m) for multiple primes"""
    if m == 0:
        return np.ones_like(lambda_p)
    if m == 1:
        return lambda_p.copy()
    
    u_prev2 = np.ones_like(lambda_p)
    u_prev1 = lambda_p.copy()
    
    for _ in range(2, m + 1):
        u = lambda_p * u_prev1 - u_prev2
        u_prev2, u_prev1 = u_prev1, u
    
    return u_prev1
