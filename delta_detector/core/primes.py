import numpy as np
from math import log, isqrt
from typing import List, Dict, Tuple
import pickle
import os

CACHE_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "cache")


def primes_up_to(n: int) -> List[int]:
    """
    Simple sieve of Eratosthenes
    Fast enough up to ~10^7 on a modern laptop
    """
    if n < 2:
        return []
    
    # Create boolean array
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    
    # Sieve
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            step = p
            start = p * p
            sieve[start:n + 1:step] = b"\x00" * (((n - start) // step) + 1)
    
    # Extract primes
    return [i for i in range(2, n + 1) if sieve[i]]


def primes_up_to_numpy(n: int) -> np.ndarray:
    """Numpy version of prime sieve for better performance"""
    if n < 2:
        return np.array([], dtype=np.int32)
    
    sieve = np.ones(n + 1, dtype=bool)
    sieve[0:2] = False
    
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = False
    
    return np.where(sieve)[0].astype(np.int32)


def prime_logs(primes: List[int]) -> Dict[int, float]:
    """Compute logarithms of primes"""
    return {p: log(p) for p in primes}


def prime_logs_numpy(primes: np.ndarray) -> Tuple[np.ndarray, Dict[int, float]]:
    """Numpy version returning both array and dict"""
    logs = np.log(primes.astype(float))
    log_dict = {int(p): float(l) for p, l in zip(primes, logs)}
    return logs, log_dict


def primes_and_logs(P: int) -> Tuple[List[int], Dict[int, float]]:
    """Get primes up to P and their logarithms"""
    ps = primes_up_to(P)
    return ps, prime_logs(ps)


class PrimeCache:
    """Cache for primes and their logarithms"""
    
    def __init__(self):
        self.primes: List[int] = []
        self.logs: Dict[int, float] = {}
        self.max_computed = 0
        
    def get_primes_up_to(self, n: int) -> Tuple[List[int], Dict[int, float]]:
        """Get primes up to n, computing if necessary"""
        if n > self.max_computed:
            self.primes = primes_up_to(n)
            self.logs = prime_logs(self.primes)
            self.max_computed = n
        
        # Return only primes <= n
        filtered_primes = [p for p in self.primes if p <= n]
        filtered_logs = {p: self.logs[p] for p in filtered_primes}
        return filtered_primes, filtered_logs
    
    def save_to_file(self, filename: str):
        """Save cache to pickle file"""
        os.makedirs(CACHE_DIR, exist_ok=True)
        filepath = os.path.join(CACHE_DIR, filename)
        with open(filepath, 'wb') as f:
            pickle.dump({
                'primes': self.primes,
                'logs': self.logs,
                'max_computed': self.max_computed
            }, f)
    
    def load_from_file(self, filename: str) -> bool:
        """Load cache from pickle file, return True if successful"""
        filepath = os.path.join(CACHE_DIR, filename)
        if os.path.exists(filepath):
            try:
                with open(filepath, 'rb') as f:
                    data = pickle.load(f)
                self.primes = data['primes']
                self.logs = data['logs']
                self.max_computed = data['max_computed']
                return True
            except:
                return False
        return False


# Global cache instance
_prime_cache = PrimeCache()

def get_primes_cached(P: int) -> Tuple[List[int], Dict[int, float]]:
    """Get primes using global cache"""
    return _prime_cache.get_primes_up_to(P)
