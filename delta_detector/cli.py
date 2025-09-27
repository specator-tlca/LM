import click
import numpy as np
from typing import Optional
import sys

from .core.detector import DeltaDetector
from .core.params import Defaults
from .core.primes import get_primes_cached
from .data.loader import load_lambda_from_csv
from .data.sato_tate import sample_lambda_ST
from .analysis.sign_map import SignMapper
from .analysis.lehmer import LehmerTester
from .visualize.plots import plot_convergence


@click.group()
@click.version_option(version='0.1.0')
def main():
    """Delta Detector - Frequency-based Hecke eigenvalue detection"""
    pass


@main.command()
@click.argument('prime', type=int)
@click.option('--sigma', default=1.3, help='Real part of s')
@click.option('--a', default=1, help='Power in λ(p^a)')
@click.option('--T-mul', default=20.0, help='Window size multiplier')
@click.option('--P', default=1000000, help='Prime cutoff')
@click.option('--csv', type=click.Path(exists=True), help='CSV file with λ(p) values')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
def prime(prime, sigma, a, t_mul, p, csv, verbose):
    """Detect eigenvalue λ(p^a) for a single prime"""
    click.echo(f"Detecting λ({prime}^{a})...")
    
    # Create detector
    detector = DeltaDetector(
        p0=prime,
        a=a,
        sigma=sigma,
        T_mul=t_mul,
        P=p
    )
    
    # Load lambda values
    if csv:
        lambda_p = load_lambda_from_csv(csv)
    else:
        click.echo("No CSV provided, using Sato-Tate simulation")
        primes, _ = get_primes_cached(p)
        lambda_p = sample_lambda_ST(primes)
    
    detector.set_lambda_values(lambda_p)
    
    # Compute estimate
    result = detector.estimate_lambda()
    
    # Display results
    click.echo(f"\nResults for p={prime}, a={a}:")
    click.echo(f"λ({prime}^{a}) ≈ {result.lambda_estimate:.12f}")
    click.echo(f"Sign: {result.sign}")
    click.echo(f"Confidence: {result.confidence:.2%}")
    
    if verbose:
        click.echo(f"\nDiagnostics:")
        click.echo(f"T = {result.T:.2f}")
        click.echo(f"J(T) = {result.J_T:.12e}")
        click.echo(f"A(T) = {result.A_T:.12e}")
        if result.J_2T is not None:
            click.echo(f"J(2T) = {result.J_2T:.12e}")
            click.echo(f"A(2T) = {result.A_2T:.12e}")
        click.echo(f"c(T) = {result.c_T:.6f}")
        click.echo(f"σ = {result.sigma}, P = {result.P}")


@main.command()
@click.option('--start', default=1000, help='Start of prime range')
@click.option('--end', default=10000, help='End of prime range')
@click.option('--output', '-o', type=click.Path(), help='Output CSV file')
@click.option('--csv', type=click.Path(exists=True), help='CSV file with λ(p) values')
@click.option('--parallel', '-p', is_flag=True, help='Use parallel processing')
@click.option('--sigma', default=1.3, help='Real part of s')
def scan(start, end, output, csv, parallel, sigma):
    """Scan a range of primes and detect signs"""
    click.echo(f"Scanning primes from {start} to {end}...")
    
    mapper = SignMapper(sigma=sigma)
    
    # Load lambda values if provided
    if csv:
        lambda_p = load_lambda_from_csv(csv)
        mapper.set_lambda_values(lambda_p)
    
    # Compute signs
    results = mapper.compute_range(start, end, parallel=parallel)
    
    # Display summary
    total = len(results)
    positive = sum(1 for r in results if r['sign'] == '+')
    negative = sum(1 for r in results if r['sign'] == '-')
    zero = sum(1 for r in results if r['sign'] == '0')
    
    click.echo(f"\nProcessed {total} primes:")
    click.echo(f"Positive: {positive} ({100*positive/total:.1f}%)")
    click.echo(f"Negative: {negative} ({100*negative/total:.1f}%)")
    click.echo(f"Zero: {zero} ({100*zero/total:.1f}%)")
    
    # Save results
    if output:
        mapper.save_results(output)
        click.echo(f"\nResults saved to {output}")


@main.command()
@click.option('--limit', default=100000, help='Check primes up to this limit')
@click.option('--csv', type=click.Path(exists=True), help='CSV file with λ(p) values')
@click.option('--threshold', default=1e-6, help='Zero threshold')
@click.option('--output', '-o', type=click.Path(), help='Output file for candidates')
def lehmer(limit, csv, threshold, output):
    """Test Lehmer's hypothesis: search for τ(p) = 0"""
    click.echo(f"Testing Lehmer hypothesis up to {limit}...")
    
    tester = LehmerTester(threshold=threshold)
    
    # Load lambda values if provided
    if csv:
        lambda_p = load_lambda_from_csv(csv)
        tester.set_lambda_values(lambda_p)
    
    # Run test
    candidates = tester.find_candidates(limit)
    
    # Display results
    if candidates:
        click.echo(f"\nFound {len(candidates)} candidates:")
        for c in candidates[:10]:  # Show first 10
            click.echo(f"p = {c['prime']}: λ(p) ≈ {c['lambda_p']:.2e}, λ(p²) ≈ {c['lambda_p2']:.6f}")
        if len(candidates) > 10:
            click.echo(f"... and {len(candidates)-10} more")
    else:
        click.echo("\nNo candidates found (good - supports Lehmer!)")
    
    # Save if requested
    if output and candidates:
        tester.save_candidates(output, candidates)
        click.echo(f"\nCandidates saved to {output}")


@main.command()
@click.argument('prime', type=int)
@click.option('--type', 'plot_type', type=click.Choice(['convergence', 'stability']), 
              default='convergence', help='Type of visualization')
@click.option('--csv', type=click.Path(exists=True), help='CSV file with λ(p) values')
@click.option('--output', '-o', type=click.Path(), help='Save plot to file')
def visualize(prime, plot_type, csv, output):
    """Visualize detector behavior"""
    
    # Create detector
    detector = DeltaDetector(p0=prime)
    
    # Load lambda values
    if csv:
        lambda_p = load_lambda_from_csv(csv)
    else:
        primes, _ = get_primes_cached(1000000)
        lambda_p = sample_lambda_ST(primes)
    
    detector.set_lambda_values(lambda_p)
    
    # Generate plot
    if plot_type == 'convergence':
        plot_convergence(detector, save_to=output)
    else:  # stability
        results = detector.analyze_stability()
        click.echo(f"Stability analysis for p={prime}:")
        click.echo(f"Mean estimate: {results['mean']:.6f}")
        click.echo(f"Std deviation: {results['std']:.6f}")
        click.echo(f"CV: {results['cv']:.2%}")
        click.echo(f"Sign stable: {results['sign_stable']}")


if __name__ == '__main__':
    main()
