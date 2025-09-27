"""Core computation modules for Delta Detector"""

from .detector import DeltaDetector, DetectorResult
from .params import Defaults, default_T, default_h
from .weights import W_weight, c_at_resonance, Phi, Phi_k
from .chebyshev import lambda_p_power, LambdaPMCache
from .primes import primes_up_to, get_primes_cached

__all__ = [
    "DeltaDetector",
    "DetectorResult", 
    "Defaults",
    "default_T",
    "W_weight",
    "c_at_resonance",
    "lambda_p_power",
    "primes_up_to",
    "get_primes_cached"
]
