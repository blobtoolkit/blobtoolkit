"""
A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

import io
from os.path import dirname
from os.path import join

from setuptools import find_namespace_packages
from setuptools import setup


def read(*names, **kwargs):
    """Read file."""
    with io.open(
        join(dirname(__file__), *names), encoding=kwargs.get("encoding", "utf8")
    ) as fh:
        return fh.read()


setup(
    name="blobtoolkit-data",  # Required
    version="3.5.5",
    description="blobtoolkit-data",  # Optional
    long_description="blobtoolkit-data",  # Optional
    long_description_content_type="text/markdown",
    url="https://github.com/blobtoolkit/blobtoolkit",  # Optional
    # This should be your name or the name of the organization which owns the
    # project.
    author="blobtoolkit",  # Optional
    # This should be a valid email address corresponding to the author listed
    # above.
    author_email="blobtoolkit@genomehubs.org",  # Optional
    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 5 - Production/Stable",
        # Indicate who your project is intended for
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        # Pick your license as you wish
        "License :: OSI Approved :: MIT License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="bioinformatics",
    package_dir={"": "src"},
    py_modules=["blobtoolkit_data"],
    packages=find_namespace_packages(
        where="src",
        exclude=[],
    ),
    python_requires=">=3.7, <=3.11",
    install_requires=[
        "docopt>=0.6.2",
    ],
    entry_points={
        "console_scripts": [
            "blobtoolkit-data = blobtoolkit-data:cli",
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/blobtoolkit/blobtoolkit/issues",
        "Source": "https://github.com/blobtoolkit/blobtoolkit",
    },
    include_package_data=True,
    zip_safe=False,
)
