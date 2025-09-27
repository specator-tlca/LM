import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, List
from ..core import default_T


def plot_convergence(detector, T_range: Optional[List[float]] = None, save_to: Optional[str] = None):
    """Plot convergence of J(T) and A(T) as functions of T"""
    
    if T_range is None:
        T_base = default_T(detector.p0, detector.a, detector.sigma, detector.T_mul)
        T_range = np.linspace(T_base * 0.1, T_base * 4, 50)
    
    J_values = []
    A_values = []
    
    for T in T_range:
        J = detector.compute_J_weighted(T)
        J_values.append(J)
        A_values.append(J / T)
    
    plt.figure(figsize=(12, 5))
    
    # Plot J(T)
    plt.subplot(121)
    plt.plot(T_range, J_values, 'b-', linewidth=2)
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    plt.xlabel('T')
    plt.ylabel('J(T)')
    plt.title(f'J(T) for p={detector.p0}')
    plt.grid(True, alpha=0.3)
    
    # Plot A(T) = J(T)/T
    plt.subplot(122)
    plt.plot(T_range, A_values, 'r-', linewidth=2)
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    plt.xlabel('T')
    plt.ylabel('A(T) = J(T)/T')
    plt.title(f'A(T) for p={detector.p0}')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_to:
        plt.savefig(save_to, dpi=150)
    else:
        plt.show()


def plot_sign_map(signs_data: List[dict], save_to: Optional[str] = None):
    """Plot sign map from sign detection results"""
    
    primes = [d['prime'] for d in signs_data]
    signs = [d['sign'] for d in signs_data]
    confidences = [d['confidence'] for d in signs_data]
    
    # Convert signs to colors
    colors = []
    for s in signs:
        if s == '+':
            colors.append('blue')
        elif s == '-':
            colors.append('red')
        else:
            colors.append('gray')
    
    plt.figure(figsize=(12, 6))
    
    # Scatter plot with confidence as size
    sizes = [100 * c for c in confidences]
    plt.scatter(primes, [1]*len(primes), c=colors, s=sizes, alpha=0.6)
    
    plt.xlabel('Prime p')
    plt.title('Sign map of λ(p)')
    plt.yticks([])
    plt.grid(True, axis='x', alpha=0.3)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', label='λ(p) > 0'),
        Patch(facecolor='red', label='λ(p) < 0'),
        Patch(facecolor='gray', label='λ(p) ≈ 0')
    ]
    plt.legend(handles=legend_elements)
    
    if save_to:
        plt.savefig(save_to, dpi=150)
    else:
        plt.show()


def plot_stability_analysis(detector, save_to: Optional[str] = None):
    """Plot stability analysis results"""
    
    results = detector.analyze_stability()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # T variation
    ax1.plot(results['T_values'], results['estimates'], 'o-', markersize=8)
    ax1.axhline(y=results['mean'], color='r', linestyle='--', 
                label=f'Mean: {results["mean"]:.6f}')
    ax1.fill_between(results['T_values'], 
                     results['mean'] - results['std'],
                     results['mean'] + results['std'],
                     alpha=0.3, color='r')
    ax1.set_xlabel('T')
    ax1.set_ylabel('λ estimate')
    ax1.set_title(f'Stability across T for p={detector.p0}')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # P variation
    P_factors = list(results['P_stability'].keys())
    P_estimates = list(results['P_stability'].values())
    
    ax2.bar(P_factors, P_estimates)
    ax2.axhline(y=results['mean'], color='r', linestyle='--')
    ax2.set_xlabel('P multiplier')
    ax2.set_ylabel('λ estimate')
    ax2.set_title('Stability across P')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    if save_to:
        plt.savefig(save_to, dpi=150)
    else:
        plt.show()
