import pytest
import numpy as np
from delta_detector.core import (
    DeltaDetector, 
    lambda_p_power,
    primes_up_to,
    get_primes_cached,
    W_weight,
    c_at_resonance
)
from delta_detector.data import sample_lambda_ST
import math


def test_lambda_p_power():
    """Test Chebyshev recursion"""
    # For λ(p) = 2, we have λ(p²) = 2² - 2 = 2
    assert lambda_p_power(2.0, 0) == 1.0
    assert lambda_p_power(2.0, 1) == 2.0
    assert lambda_p_power(2.0, 2) == 2.0
    assert lambda_p_power(2.0, 3) == 2.0
    
    # For λ(p) = 0, we have λ(p²) = 0 - 1 = -1
    assert lambda_p_power(0.0, 0) == 1.0
    assert lambda_p_power(0.0, 1) == 0.0
    assert lambda_p_power(0.0, 2) == -1.0
    assert lambda_p_power(0.0, 3) == 0.0


def test_primes_up_to():
    """Test prime generation"""
    primes = primes_up_to(30)
    expected = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert primes == expected


def test_weight_function():
    """Test weight function properties"""
    # At resonance x = a*log(p0)
    a_logp0 = 1.0 * math.log(29)
    T = 100.0
    
    # Weight at resonance should be maximum
    w_resonance = W_weight(a_logp0, a_logp0, T)
    w_off = W_weight(a_logp0 + 1, a_logp0, T)
    
    assert w_resonance > w_off
    assert 0.5 <= w_resonance <= 1.0
    assert 0 <= w_off <= 1.0
    
    # c(T) calculation
    c_T = c_at_resonance(a_logp0, T)
    assert abs(c_T - w_resonance) < 1e-10


def test_sato_tate_distribution():
    """Test Sato-Tate sampling"""
    primes = primes_up_to(1000)
    lambda_values = sample_lambda_ST(primes, seed=42)
    
    # Check bounds
    for p, lam in lambda_values.items():
        assert -2 <= lam <= 2
    
    # Check distribution (roughly)
    values = list(lambda_values.values())
    assert -0.2 < np.mean(values) < 0.2  # Should be near 0
    assert 0.8 < np.std(values) < 1.2   # Should be near 1


def test_detector_basic():
    """Test basic detector functionality"""
    detector = DeltaDetector(p0=29, sigma=1.3)
    
    # Generate test lambda values
    primes, _ = get_primes_cached(10000)
    lambda_values = sample_lambda_ST(primes, seed=42)
    detector.set_lambda_values(lambda_values)
    
    # Compute J
    T = 100.0
    J = detector.compute_J_weighted(T)
    
    # Should get a finite value
    assert np.isfinite(J)
    
    # Test estimate
    result = detector.estimate_lambda(T=T)
    
    assert result.p0 == 29
    assert result.a == 1
    assert result.sigma == 1.3
    assert result.T == T
    assert result.sign in ['+', '-', '0']
    assert 0 <= result.confidence <= 1
    assert np.isfinite(result.lambda_estimate)


def test_stability_analysis():
    """Test stability analysis"""
    detector = DeltaDetector(p0=29, sigma=1.3)
    
    # Set simple lambda values
    primes, _ = get_primes_cached(10000)
    lambda_values = {p: 1.0 for p in primes}  # All positive
    detector.set_lambda_values(lambda_values)
    
    # Analyze stability
    stability = detector.analyze_stability(T_range=[50, 100, 200])
    
    assert 'mean' in stability
    assert 'std' in stability
    assert 'cv' in stability
    assert 'sign_stable' in stability
    assert stability['sign_stable'] == True  # All should be positive


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
