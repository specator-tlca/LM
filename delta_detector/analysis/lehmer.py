import numpy as np
from typing import Dict, List, Optional
import csv
from tqdm import tqdm

from ..core import DeltaDetector, primes_up_to
from ..data import sample_lambda_ST


class LehmerTester:
    """
    Test Lehmer's hypothesis: τ(p) ≠ 0 for all primes p
    
    If λ(p) = 0, then λ(p²) = -1 (from Chebyshev recursion)
    So we look for primes where λ(p) ≈ 0 and verify λ(p²) ≈ -1
    """
    
    def __init__(
        self,
        threshold: float = 1e-6,
        sigma: float = 1.3,
        T_mul: float = 20.0,
        P: int = 1_000_000
    ):
        """
        Args:
            threshold: How close to 0 counts as "zero"
            sigma: Real part for detector
            T_mul: Window size multiplier
            P: Prime cutoff
        """
        self.threshold = threshold
        self.sigma = sigma
        self.T_mul = T_mul
        self.P = P
        self.lambda_p: Dict[int, float] = {}
        
    def set_lambda_values(self, lambda_p: Dict[int, float]):
        """Set the lambda values to use"""
        self.lambda_p = lambda_p
    
    def test_prime(self, p: int) -> Dict:
        """
        Test a single prime for Lehmer property
        
        Returns:
            Dictionary with test results
        """
        # Test λ(p) with a=1
        detector1 = DeltaDetector(
            p0=p,
            a=1,
            sigma=self.sigma,
            T_mul=self.T_mul,
            P=self.P
        )
        
        # Set lambda values
        if self.lambda_p:
            detector1.set_lambda_values(self.lambda_p)
        else:
            primes, _ = detector1.get_primes_cached(self.P)
            lambda_vals = sample_lambda_ST(primes)
            detector1.set_lambda_values(lambda_vals)
        
        # Estimate λ(p)
        result1 = detector1.estimate_lambda()
        
        # Check if potentially zero
        is_candidate = abs(result1.lambda_estimate) < self.threshold
        
        result = {
            'prime': p,
            'lambda_p': result1.lambda_estimate,
            'sign_p': result1.sign,
            'is_candidate': is_candidate,
            'lambda_p2': None,
            'sign_p2': None,
            'verified': False
        }
        
        # If candidate, verify with λ(p²)
        if is_candidate:
            detector2 = DeltaDetector(
                p0=p,
                a=2,  # Check p²
                sigma=self.sigma,
                T_mul=self.T_mul,
                P=self.P
            )
            detector2.set_lambda_values(detector1.lambda_p)  # Use same values
            
            result2 = detector2.estimate_lambda()
            
            result['lambda_p2'] = result2.lambda_estimate
            result['sign_p2'] = result2.sign
            # Verify: should be close to -1
            result['verified'] = abs(result2.lambda_estimate + 1) < 0.1
        
        return result
    
    def find_candidates(
        self,
        limit: int,
        show_progress: bool = True
    ) -> List[Dict]:
        """
        Search for Lehmer candidates up to limit
        
        Args:
            limit: Check primes up to this value
            show_progress: Show progress bar
            
        Returns:
            List of candidate primes
        """
        primes = primes_up_to(limit)
        candidates = []
        
        if show_progress:
            pbar = tqdm(primes, desc="Testing Lehmer hypothesis")
        else:
            pbar = primes
        
        for p in pbar:
            result = self.test_prime(p)
            
            if result['is_candidate']:
                candidates.append(result)
                if show_progress:
                    tqdm.write(f"Candidate found: p={p}, λ(p)≈{result['lambda_p']:.2e}")
        
        return candidates
    
    def save_candidates(self, filepath: str, candidates: List[Dict]):
        """Save candidates to CSV file"""
        if not candidates:
            return
            
        with open(filepath, 'w', newline='') as f:
            fieldnames = ['prime', 'lambda_p', 'sign_p', 'lambda_p2', 'sign_p2', 'verified']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for c in candidates:
                row = {k: c[k] for k in fieldnames}
                writer.writerow(row)
    
    def analyze_candidates(self, candidates: List[Dict]) -> Dict:
        """Analyze found candidates"""
        if not candidates:
            return {
                'count': 0,
                'verified': 0,
                'false_positives': 0
            }
        
        verified = sum(1 for c in candidates if c['verified'])
        
        # Check distribution of λ(p²) values
        lambda_p2_values = [c['lambda_p2'] for c in candidates if c['lambda_p2'] is not None]
        
        stats = {
            'count': len(candidates),
            'verified': verified,
            'false_positives': len(candidates) - verified,
            'smallest_lambda_p': min(abs(c['lambda_p']) for c in candidates),
            'largest_prime': max(c['prime'] for c in candidates)
        }
        
        if lambda_p2_values:
            stats['mean_lambda_p2'] = np.mean(lambda_p2_values)
            stats['std_lambda_p2'] = np.std(lambda_p2_values)
        
        return stats


def run_lehmer_test(
    limit: int = 1_000_000,
    threshold: float = 1e-6,
    output_file: Optional[str] = None,
    lambda_csv: Optional[str] = None
) -> Dict:
    """
    Convenience function to run complete Lehmer test
    
    Returns:
        Dictionary with candidates and statistics
    """
    tester = LehmerTester(threshold=threshold)
    
    # Load lambda values if provided
    if lambda_csv:
        from ..data import load_lambda_from_csv
        lambda_p = load_lambda_from_csv(lambda_csv)
        tester.set_lambda_values(lambda_p)
    
    print(f"Testing Lehmer hypothesis up to {limit}")
    print(f"Zero threshold: {threshold}")
    
    # Find candidates
    candidates = tester.find_candidates(limit)
    
    # Analyze
    stats = tester.analyze_candidates(candidates)
    
    print(f"\nResults:")
    print(f"Candidates found: {stats['count']}")
    if stats['count'] > 0:
        print(f"Verified (λ(p²)≈-1): {stats['verified']}")
        print(f"False positives: {stats['false_positives']}")
        print(f"Smallest |λ(p)|: {stats['smallest_lambda_p']:.2e}")
    
    # Save if requested
    if output_file and candidates:
        tester.save_candidates(output_file, candidates)
        print(f"\nCandidates saved to {output_file}")
    
    # Conclusion
    if stats['count'] == 0:
        print("\nNo candidates found - supports Lehmer's hypothesis!")
    else:
        print("\nWARNING: Found potential counterexamples - needs investigation")
    
    return {
        'candidates': candidates,
        'statistics': stats
    }
