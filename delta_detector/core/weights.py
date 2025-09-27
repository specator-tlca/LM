"""
Weight functions for frequency-domain filtering using B-spline kernels.

B-spline properties by order:
- k=1 (Fejér): Triangular window, sidelobe decay ~ 1/ω²
- k=2: Smoother transition, sidelobe decay ~ 1/ω⁴
- k=3: Even smoother, sidelobe decay ~ 1/ω⁶

Higher k gives:
+ Better sidelobe suppression (reduces cross-talk between frequencies)
+ Smoother frequency response
- Wider main lobe (less frequency resolution)
- More computation (higher powers)

Recommended values:
- k=1: Default, good balance
- k=2: When λ values are close in frequency space
- k=3: For maximum isolation, when computation time allows
"""

import numpy as np
from math import sin, pi

def sinc(x: float) -> float:
    """Safe sinc function: sin(x)/x"""
    if x == 0.0:
        return 1.0
    return sin(x) / x


def Phi(x: float) -> float:
    """Fejér kernel in frequency domain: (sin(ξ/2)/(ξ/2))^2"""
    return sinc(0.5 * x) ** 2


def Phi_k(x: float, k: int) -> float:
    """B-spline of order k in frequency domain: (sinc(ξ/2))^(2k)"""
    if k <= 1:
        return Phi(x)
    s = sinc(0.5 * x)
    return s ** (2 * k)


def W_weight(x: float, a_logp0: float, T: float, k: int = 1) -> float:
    """
    Analytic weight function for frequency x
    W(x) = 0.5 * (Φ((x - a*log(p0))*T) + Φ((x + a*log(p0))*T))
    """
    return 0.5 * (Phi_k((x - a_logp0) * T, k) + Phi_k((x + a_logp0) * T, k))


def c_at_resonance(a_logp0: float, T: float, k: int = 1) -> float:
    """
    Value of weight at resonance point x = a*log(p0)
    c_{a,p0}(T) = 0.5 * (1 + Φ(2*a*log(p0)*T))
    """
    return 0.5 * (1.0 + Phi_k(2.0 * a_logp0 * T, k))


def W_weight_vectorized(x: np.ndarray, a_logp0: float, T: float, k: int = 1) -> np.ndarray:
    """Vectorized version of W_weight for numpy arrays"""
    def sinc_vec(y):
        # numpy's sinc(x) = sin(pi*x)/(pi*x), but we need sin(x)/x
        # so we divide by pi to get the correct normalization
        return np.sinc(y / np.pi)
    
    def Phi_k_vec(y, order):
        s = sinc_vec(0.5 * y)
        return s ** (2 * order)
    
    return 0.5 * (Phi_k_vec((x - a_logp0) * T, k) + Phi_k_vec((x + a_logp0) * T, k))
