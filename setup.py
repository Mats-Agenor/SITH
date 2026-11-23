from setuptools import find_packages, setup

setup(
    name="sith",
    version="2.0.0",
    description="SITH 2.0 (S-PLUS Integrated Tool for Homogeneous Photometry)",
    packages=find_packages(),
    include_package_data=True,
    package_data={"sith": ["resources/*", "external/*"]},
    install_requires=[
        "numpy>=1.24",
        "astropy>=5.3",
        "astropy-healpix>=1.0",
        "scipy>=1.10",
        "scikit-learn>=1.2",
        "matplotlib>=3.7",
    ],
    python_requires=">=3.10",
    entry_points={"console_scripts": ["sith=sith.cli:main"]},
)
