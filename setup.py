from setuptools import find_packages, setup

setup(
    name="pyranges_plot",
    version="0.0.24",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
)
