from setuptools import find_packages, setup

setup(
    name="pyranges_plot",
    version="0.1.2",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,  # Ensure package data is included
    package_data={
        "pyranges_plot": ["data/*"],  # Specify the path to include data folder contents
    },
)
