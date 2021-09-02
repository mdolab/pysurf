from setuptools import setup
import re
import os

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("tsurf/__init__.py").read(),
)[0]

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="tsurf",
    version=__version__,
    description="TSurf provides geometric operations for on triangulated surfaces",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="intersection geometry projection",
    author="",
    author_email="",
    url="https://github.com/mdolab/pysurf",
    license="Apache License 2.0",
    packages=[
        "tsurf",
    ],
    package_data={"tsurf": ["*.so"]},
    install_requires=[
        "numpy>=1.16",
        "mpi4py>=3.0",
    ],
    extras_require={"testing": ["testflo"]},
    classifiers=["Operating System :: Linux", "Programming Language :: Python, Fortran"],
)
