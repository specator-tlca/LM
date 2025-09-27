import numpy as np
from math import log, pow
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import warnings

from .params import Defaults, default_T
from .primes import get_primes_cached
from .weights import W_weight, c_at_resonance, W_weight_vectorized
from .chebyshev import LambdaPMCache


@dataclass
class DetectorResult:
    """Result of eigenvalue detection"""
    lambda_estimate: float
    sign: str
    confidence: float
    J_T: float
    J_2T: Optional[float]
    A_T: float
    A_2T: Optional[float]
    c_T: float
    T: float
    p0: int
    a: int
    sigma: float
    P: int


class DeltaDetector:
    """
    Frequency-based detector for Hecke eigenvalues λ(p0^a)
    
    Uses bandpass filtering in frequency domain to extract eigenvalues
    from the logarithmic derivative of L-functions.
    """
    
    def __init__(
        self,
        p0: int,
        a: int = 1,
        sigma: float = 1.3,
        T_mul: float = 20.0,
        M: int = 2,
        P: int = 1_000_000,
        kernel: str = "fejer",
        k: int = 1,
        richardson: bool = True
    ):
        """
        Initialize detector
        
        Args:
            p0: Target prime
            a: Power (for λ(p0^a))
            sigma: Real part of s = σ + it (must be > 1)
            T_mul: Multiplier for default T
            M: Maximum power m in sum over p^m
            P: Prime cutoff
            kernel: Window kernel type ('fejer' or 'bspline')
            k: Order for B-spline kernel
            richardson: Use Richardson extrapolation
        """
        self.p0 = p0
        self.a = a
        self.sigma = sigma
        self.T_mul = T_mul
        self.M = M
        self.P = P
        self.kernel = kernel
        self.k = k
        self.richardson = richardson
        
        # Precompute some values
        self.a_logp0 = a * log(p0)
        self.cache = LambdaPMCache()
        
        # Lambda values (to be loaded)
        self.lambda_p: Dict[int, float] = {}
        
    def set_lambda_values(self, lambda_p: Dict[int, float]):
        """Set the Hecke eigenvalues λ(p)"""
        self.lambda_p = lambda_p
        self.cache.clear()
        
    def compute_J_weighted(self, T: float, P: Optional[int] = None) -> float:
        """
        Compute J using weighted sum (fast method)
        
        J_{a,p0}(σ,T) = T * Σ_{p,m} λ(p^m)/(m*p^{mσ}) * W(m*log(p))
        """
        if P is None:
            P = self.P
            
        primes, logp = get_primes_cached(P)
        
        if not self.lambda_p:
            warnings.warn("No lambda values loaded, using zeros")
            
        S = 0.0
        
        for p in primes:
            lp = logp[p]
            lam1 = self.lambda_p.get(p, 0.0)
            
            for m in range(1, self.M + 1):
                lam_pm = self.cache.get(p, m, lam1)
                denom = m * pow(p, m * self.sigma)
                x = m * lp
                w = W_weight(x, self.a_logp0, T, self.k)
                S += (lam_pm / denom) * w
        
        return T * S
    
    def compute_J_weighted_vectorized(self, T: float, P: Optional[int] = None) -> float:
        """Vectorized version for better performance with numpy"""
        if P is None:
            P = self.P
            
        primes, logp = get_primes_cached(P)
        primes_array = np.array(primes)
        
        if not self.lambda_p:
            warnings.warn("No lambda values loaded, using zeros")
            return 0.0
        
        # Vectorize lambda values
        lambda_array = np.array([self.lambda_p.get(p, 0.0) for p in primes])
        log_array = np.array([logp[p] for p in primes])
        
        S = 0.0
        
        for m in range(1, self.M + 1):
            # Compute λ(p^m) for all primes at once
            if m == 1:
                lam_pm_array = lambda_array
            else:
                lam_pm_array = np.array([
                    self.cache.get(p, m, lam) 
                    for p, lam in zip(primes, lambda_array)
                ])
            
            # Compute denominators
            denom_array = m * np.power(primes_array, m * self.sigma)
            
            # Compute weights
            x_array = m * log_array
            w_array = W_weight_vectorized(x_array, self.a_logp0, T, self.k)
            
            # Sum contribution
            S += np.sum((lam_pm_array / denom_array) * w_array)
        
        return T * S
    
    def estimate_lambda(
        self,
        T: Optional[float] = None,
        use_vectorized: bool = True
    ) -> DetectorResult:
        """
        Estimate λ(p0^a) using the detector
        
        Args:
            T: Window size (None for default)
            use_vectorized: Use numpy vectorized computation
            
        Returns:
            DetectorResult with estimate and diagnostics
        """
        if T is None:
            T = default_T(self.p0, self.a, self.sigma, self.T_mul)
        
        # Choose computation method
        compute_J = self.compute_J_weighted_vectorized if use_vectorized else self.compute_J_weighted
        
        # Compute J(T)
        J_T = compute_J(T)
        A_T = J_T / T
        
        # Richardson extrapolation
        if self.richardson:
            J_2T = compute_J(2 * T)
            A_2T = J_2T / (2 * T)
            
            # λ̂ = 2*a*p0^{aσ} * (2*A(2T) - A(T))
            lambda_hat = 2 * self.a * pow(self.p0, self.a * self.sigma) * (2 * A_2T - A_T)
        else:
            J_2T = None
            A_2T = None
            # Simple estimate
            c_T = c_at_resonance(self.a_logp0, T, self.k)
            lambda_hat = self.a * pow(self.p0, self.a * self.sigma) * A_T / c_T
        
        # Determine sign
        if abs(J_T) < 1e-12:
            sign = "0"
        else:
            sign = "+" if J_T > 0 else "-"
        
        # Compute confidence (heuristic)
        c_T = c_at_resonance(self.a_logp0, T, self.k)
        if self.richardson and A_2T is not None:
            # Stability of Richardson
            stability = abs(A_2T - A_T) / max(abs(A_T), 1e-12)
            confidence = min(0.99, 1.0 / (1.0 + stability))
        else:
            # Based on resonance strength
            confidence = c_T
        
        return DetectorResult(
            lambda_estimate=lambda_hat,
            sign=sign,
            confidence=confidence,
            J_T=J_T,
            J_2T=J_2T,
            A_T=A_T,
            A_2T=A_2T,
            c_T=c_T,
            T=T,
            p0=self.p0,
            a=self.a,
            sigma=self.sigma,
            P=self.P
        )
    
    def compute_sin_test(self, T: float, P: Optional[int] = None) -> float:
        """
        Sin-modulated integral - should remain O(1)
        
        Uses sin instead of cos in the modulation. If this grows with T,
        it indicates problems with parameters or implementation.
        
        Returns:
            float: Sin test value (should be O(1))
        """
        if P is None:
            P = self.P
            
        primes, logp = get_primes_cached(P)
        
        if not self.lambda_p:
            warnings.warn("No lambda values loaded, using zeros")
            return 0.0
            
        S = 0.0
        
        for p in primes:
            lp = logp[p]
            lam1 = self.lambda_p.get(p, 0.0)
            
            for m in range(1, self.M + 1):
                lam_pm = self.cache.get(p, m, lam1)
                denom = m * pow(p, m * self.sigma)
                
                # For sin test, we need to modify the weight function
                # The key is that sin is orthogonal to cos in our integral
                # Use fixed number of points for numerical integration
                N_POINTS = 1001  # Fixed number of integration points
                t_vals = np.linspace(-T, T, N_POINTS)
                dt = 2.0 * T / (N_POINTS - 1)
                
                # Compute integral with sin modulation
                integral = 0.0
                for t in t_vals:
                    # Window function
                    window = max(0, 1 - abs(t)/T) if abs(t) <= T else 0
                    # Sin modulation (instead of cos)
                    modulation = np.sin(t * self.a_logp0)
                    # Term from Re log L
                    term = np.cos(t * m * lp)
                    integral += window * modulation * term * dt
                
                S += (lam_pm / denom) * integral
        
        return S
    
    def analyze_stability(
        self,
        T_range: Optional[List[float]] = None,
        P_factors: List[float] = [1.0, 2.0]
    ) -> Dict:
        """
        Analyze stability of estimates across different T and P values
        
        Returns dictionary with stability metrics
        """
        if T_range is None:
            T_base = default_T(self.p0, self.a, self.sigma, self.T_mul)
            T_range = [T_base * f for f in [0.5, 1.0, 2.0, 4.0]]
        
        results = {
            'T_values': T_range,
            'estimates': [],
            'signs': [],
            'P_stability': {}
        }
        
        # Vary T
        for T in T_range:
            res = self.estimate_lambda(T)
            results['estimates'].append(res.lambda_estimate)
            results['signs'].append(res.sign)
        
        # Check P stability
        for factor in P_factors:
            P_test = int(self.P * factor)
            res = self.estimate_lambda(use_vectorized=True)
            results['P_stability'][f'P*{factor}'] = res.lambda_estimate
        
        # Compute variation metrics
        estimates = np.array(results['estimates'])
        results['mean'] = np.mean(estimates)
        results['std'] = np.std(estimates)
        results['cv'] = results['std'] / abs(results['mean']) if results['mean'] != 0 else float('inf')
        results['sign_stable'] = len(set(results['signs'])) == 1
        
        return results
    
    def estimate_lambda_adaptive(
        self,
        threshold_ratio: float = 0.01,
        max_iterations: int = 5,
        T_factor: float = 2.0,
        verbose: bool = False
    ) -> DetectorResult:
        """
        Estimate λ(p0^a) with adaptive T selection
        
        For small |λ(p)|, automatically increases T until signal is strong enough.
        
        Args:
            threshold_ratio: Minimum |J(T)|/T relative to expected scale
            max_iterations: Maximum number of T doublings
            T_factor: Factor to increase T by each iteration
            verbose: Print diagnostic information
            
        Returns:
            DetectorResult with adaptively chosen T
        """
        # Start with default T
        T = default_T(self.p0, self.a, self.sigma, self.T_mul)
        
        # Expected scale of J(T) for |λ| ~ 1
        expected_scale = T / (self.a * pow(self.p0, self.a * self.sigma))
        threshold = threshold_ratio * expected_scale
        
        best_result = None
        iterations = 0
        
        while iterations < max_iterations:
            # Compute J(T)
            result = self.estimate_lambda(T)
            
            if verbose:
                print(f"  Iteration {iterations}: T = {T:.1f}, |J(T)| = {abs(result.J_T):.6e}")
            
            # Check if signal is strong enough
            if abs(result.J_T) >= threshold:
                if verbose:
                    print(f"  Signal strong enough, using T = {T:.1f}")
                best_result = result
                break
            
            # Save result in case we hit max iterations
            if best_result is None or abs(result.J_T) > abs(best_result.J_T):
                best_result = result
            
            # Increase T
            T *= T_factor
            iterations += 1
            
            if verbose and iterations < max_iterations:
                print(f"  Signal too weak, increasing T to {T:.1f}")
        
        if iterations == max_iterations and verbose:
            print(f"  Max iterations reached, using best T = {best_result.T:.1f}")
        
        return best_result
    
    def check_P_stability(
        self,
        T: Optional[float] = None,
        P_factors: List[float] = [1.0, 2.0, 4.0],
        tolerance: float = 0.05,
        verbose: bool = False
    ) -> Dict:
        """
        Check stability of estimate with respect to prime cutoff P
        
        Args:
            T: Window size (None for default)
            P_factors: Factors to multiply P by
            tolerance: Maximum acceptable relative change
            verbose: Print diagnostic information
            
        Returns:
            Dictionary with stability analysis:
            - 'stable': bool, whether estimates are stable
            - 'estimates': list of estimates for each P
            - 'P_values': list of P values used
            - 'relative_changes': relative changes between consecutive P values
            - 'recommended_P': suggested P value if current is unstable
        """
        if T is None:
            T = default_T(self.p0, self.a, self.sigma, self.T_mul)
        
        results = {
            'P_values': [],
            'estimates': [],
            'J_values': [],
            'relative_changes': []
        }
        
        # Test different P values
        for factor in P_factors:
            P_test = int(self.P * factor)
            J_test = self.compute_J_weighted_vectorized(T, P_test)
            
            # Estimate lambda
            A_test = J_test / T
            c_T = c_at_resonance(self.a_logp0, T, self.k)
            
            if self.richardson:
                J_2T_test = self.compute_J_weighted_vectorized(2*T, P_test)
                A_2T_test = J_2T_test / (2*T)
                lambda_test = 2 * self.a * pow(self.p0, self.a * self.sigma) * (2 * A_2T_test - A_test)
            else:
                lambda_test = self.a * pow(self.p0, self.a * self.sigma) * A_test / c_T
            
            results['P_values'].append(P_test)
            results['estimates'].append(lambda_test)
            results['J_values'].append(J_test)
            
            if verbose:
                print(f"  P = {P_test}: λ = {lambda_test:.6f}, J(T) = {J_test:.6e}")
        
        # Compute relative changes
        estimates = results['estimates']
        for i in range(1, len(estimates)):
            if abs(estimates[i-1]) > 1e-10:
                rel_change = abs(estimates[i] - estimates[i-1]) / abs(estimates[i-1])
            else:
                rel_change = abs(estimates[i] - estimates[i-1])
            results['relative_changes'].append(rel_change)
        
        # Check stability
        results['stable'] = all(rc < tolerance for rc in results['relative_changes'])
        
        # Recommend P if unstable
        if not results['stable'] and len(results['relative_changes']) > 0:
            # Find first stable pair
            for i, rc in enumerate(results['relative_changes']):
                if rc < tolerance:
                    results['recommended_P'] = results['P_values'][i+1]
                    break
            else:
                # No stable pair found, recommend largest P
                results['recommended_P'] = results['P_values'][-1]
        else:
            results['recommended_P'] = self.P
        
        if verbose:
            print(f"\nStability: {'✓' if results['stable'] else '✗'}")
            if not results['stable']:
                print(f"Recommended P: {results['recommended_P']}")
        
        return results
