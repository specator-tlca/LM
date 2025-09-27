from dataclasses import dataclass
from math import log
from typing import Optional

@dataclass
class Defaults:
    """Default parameters for Delta Detector"""
    sigma: float = 1.3          # Real part (1 + epsilon)
    a: int = 1                  # Power in λ(p0^a)
    T_mul: float = 20.0         # T ≈ T_mul * p0^(a*sigma)
    M: int = 2                  # Max power m in weighted sum
    P: int = 1_000_000         # Prime cutoff (budget)
    kernel: str = "fejer"       # Window kernel type
    k: int = 1                  # Order for B-spline kernel
    richardson: bool = True     # Use Richardson extrapolation
    checkP: bool = False        # Check stability with 2P
    mode: str = "weighted"      # Computation mode
    step_c: float = 0.5         # Step size factor for product mode


def default_T(p0: int, a: int, sigma: float, T_mul: float) -> float:
    """Compute default window size T"""
    return float(T_mul) * (p0 ** (a * sigma))


def default_h(P: int, step_c: float) -> float:
    """Compute default step size for integration"""
    return step_c / max(1.0, log(max(3, P)))
