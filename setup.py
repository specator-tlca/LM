from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="delta-detector",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Frequency-based detector for Hecke eigenvalues from L-functions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/delta-detector",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.10",
    install_requires=[
        "numpy>=1.21.0",
        "scipy>=1.7.0", 
        "matplotlib>=3.4.0",
        "click>=8.0.0",
    ],
    extras_require={
        "fast": ["numba>=0.54.0"],  # для ускорения
        "dev": ["pytest>=7.0.0"],    # для тестов
    },
    entry_points={
        "console_scripts": [
            "delta-detect=delta_detector.cli:main",
        ],
    },
)
