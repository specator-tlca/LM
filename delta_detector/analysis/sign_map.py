import numpy as np
from typing import Dict, List, Optional
from concurrent.futures import ProcessPoolExecutor
import csv
from tqdm import tqdm

from ..core import DeltaDetector, primes_up_to, default_T
from ..data import sample_lambda_ST


class SignMapper:
    """Map signs of λ(p) for ranges of primes"""
    
    def __init__(
        self,
        sigma: float = 1.3,
        T_mul: float = 20.0,
        P: int = 1_000_000
    ):
        self.sigma = sigma
        self.T_mul = T_mul
        self.P = P
        self.lambda_p: Dict[int, float] = {}
        
    def set_lambda_values(self, lambda_p: Dict[int, float]):
        """Set the lambda values to use"""
        self.lambda_p = lambda_p
        
    def detect_sign(self, p: int) -> Dict:
        """Detect sign of λ(p) for a single prime"""
        detector = DeltaDetector(
            p0=p,
            sigma=self.sigma,
            T_mul=self.T_mul,
            P=self.P
        )
        
        # Use stored lambda values or generate
        if self.lambda_p:
            detector.set_lambda_values(self.lambda_p)
        else:
            primes, _ = detector.get_primes_cached(self.P)
            lambda_vals = sample_lambda_ST(primes)
            detector.set_lambda_values(lambda_vals)
        
        result = detector.estimate_lambda()
        
        return {
            'prime': p,
            'lambda': result.lambda_estimate,
            'sign': result.sign,
            'confidence': result.confidence,
            'J_T': result.J_T,
            'A_T': result.A_T
        }
    
    def compute_range(
        self,
        start: int,
        end: int,
        parallel: bool = False,
        show_progress: bool = True
    ) -> List[Dict]:
        """
        Compute signs for all primes in range [start, end]
        
        Args:
            start: Starting value (inclusive)
            end: Ending value (inclusive)
            parallel: Use multiprocessing
            show_progress: Show progress bar
            
        Returns:
            List of results
        """
        # Get primes in range
        all_primes = primes_up_to(end)
        primes_in_range = [p for p in all_primes if start <= p <= end]
        
        results = []
        
        if parallel:
            # Parallel processing
            with ProcessPoolExecutor() as executor:
                if show_progress:
                    futures = list(tqdm(
                        executor.map(self.detect_sign, primes_in_range),
                        total=len(primes_in_range),
                        desc="Computing signs"
                    ))
                else:
                    futures = list(executor.map(self.detect_sign, primes_in_range))
                results = futures
        else:
            # Sequential processing
            if show_progress:
                pbar = tqdm(primes_in_range, desc="Computing signs")
            else:
                pbar = primes_in_range
                
            for p in pbar:
                results.append(self.detect_sign(p))
        
        return results
    
    def save_results(self, filepath: str, results: Optional[List[Dict]] = None):
        """Save results to CSV file"""
        if results is None:
            results = self.results
            
        with open(filepath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'prime', 'lambda', 'sign', 'confidence', 'J_T', 'A_T'
            ])
            writer.writeheader()
            writer.writerows(results)
    
    def analyze_distribution(self, results: List[Dict]) -> Dict:
        """Analyze sign distribution"""
        signs = [r['sign'] for r in results]
        
        total = len(signs)
        positive = signs.count('+')
        negative = signs.count('-')
        zero = signs.count('0')
        
        # Compute runs (sign changes)
        runs = 0
        for i in range(1, len(signs)):
            if signs[i] != signs[i-1] and signs[i] != '0' and signs[i-1] != '0':
                runs += 1
        
        return {
            'total': total,
            'positive': positive,
            'negative': negative, 
            'zero': zero,
            'positive_ratio': positive / total if total > 0 else 0,
            'negative_ratio': negative / total if total > 0 else 0,
            'runs': runs,
            'average_run_length': total / (runs + 1) if runs > 0 else total
        }
    
    def plot_distribution(self, results: List[Dict], save_to: Optional[str] = None):
        """Plot sign distribution"""
        import matplotlib.pyplot as plt
        
        primes = [r['prime'] for r in results]
        signs = [r['sign'] for r in results]
        
        # Convert signs to numbers
        sign_values = []
        for s in signs:
            if s == '+':
                sign_values.append(1)
            elif s == '-':
                sign_values.append(-1)
            else:
                sign_values.append(0)
        
        plt.figure(figsize=(12, 6))
        
        # Top plot: sign values
        plt.subplot(211)
        plt.scatter(primes, sign_values, s=2, alpha=0.5)
        plt.ylabel('Sign of λ(p)')
        plt.yticks([-1, 0, 1], ['-', '0', '+'])
        plt.grid(True, alpha=0.3)
        plt.title('Sign Distribution of λ(p)')
        
        # Bottom plot: cumulative bias
        cumsum = np.cumsum(sign_values)
        bias = cumsum / np.arange(1, len(sign_values) + 1)
        
        plt.subplot(212)
        plt.plot(primes, bias)
        plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        plt.xlabel('Prime p')
        plt.ylabel('Cumulative bias')
        plt.grid(True, alpha=0.3)
        plt.ylim(-0.1, 0.1)
        
        plt.tight_layout()
        
        if save_to:
            plt.savefig(save_to, dpi=150)
        else:
            plt.show()


def create_sign_map(
    limit: int = 10000,
    output_file: str = 'sign_map.csv',
    plot_file: Optional[str] = None,
    lambda_csv: Optional[str] = None
) -> Dict:
    """
    Convenience function to create a complete sign map
    
    Returns:
        Dictionary with results and statistics
    """
    mapper = SignMapper()
    
    # Load lambda values if provided
    if lambda_csv:
        from ..data import load_lambda_from_csv
        lambda_p = load_lambda_from_csv(lambda_csv)
        mapper.set_lambda_values(lambda_p)
    
    # Compute signs
    print(f"Computing signs for primes up to {limit}...")
    results = mapper.compute_range(2, limit, parallel=True)
    
    # Save results
    mapper.save_results(output_file, results)
    print(f"Results saved to {output_file}")
    
    # Analyze
    stats = mapper.analyze_distribution(results)
    print(f"\nStatistics:")
    print(f"Total primes: {stats['total']}")
    print(f"Positive: {stats['positive']} ({100*stats['positive_ratio']:.1f}%)")
    print(f"Negative: {stats['negative']} ({100*stats['negative_ratio']:.1f}%)")
    print(f"Sign changes: {stats['runs']}")
    
    # Plot if requested
    if plot_file:
        mapper.plot_distribution(results, plot_file)
        print(f"Plot saved to {plot_file}")
    
    return {
        'results': results,
        'statistics': stats
    }
