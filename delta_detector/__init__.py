"""
Delta Detector - Frequency-based Hecke eigenvalue detection

A breakthrough computational tool for detecting Hecke eigenvalues λ(p^a) 
from L-functions using frequency domain analysis.
"""

__version__ = "0.1.0"
__author__ = "Your Name"

from .core.detector import DeltaDetector
from .analysis.sign_map import SignMapper
from .analysis.lehmer import LehmerTester

__all__ = ["DeltaDetector", "SignMapper", "LehmerTester"]
