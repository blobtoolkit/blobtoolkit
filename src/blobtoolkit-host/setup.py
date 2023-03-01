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
    name="blobtoolkit-host",  # Required
    version="4.1.0",
    description="blobtoolkit-host",  # Optional
    long_description="blobtoolkit-host",  # Optional
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
    # This field adds keywords for your project which will appear on the
    # project page. What does your project relate to?
    #
    # Note that this is a list of additional keywords, separated
    # by commas, to be used to assist searching for the distribution in a
    # larger catalog.
    keywords="bioinformatics",  # Optional
    # When your source code is in a subdirectory under the project root, e.g.
    # `src/`, it is necessary to specify the `package_dir` argument.
    package_dir={"": "src"},  # Optional
    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    py_modules=["blobtoolkit_host"],
    #
    # packages=find_packages(where="src"),  # Required
    packages=find_namespace_packages(
        where="src",
        exclude=[],
    ),
    # Specify which Python versions you support. In contrast to the
    # 'Programming Language' classifiers above, 'pip install' will check this
    # and refuse to install the project if the version does not match. See
    # https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
    python_requires=">=3.7, <3.12",
    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        "docopt>=0.6.2",
        "psutil==5.9.4",
    ],  # Optional
    extras_require={  # Optional
        "dev": ["pycodestyle>=2.6.0", "pydocstyle>=5.0.2", "pylint>=2.5.3"],
        "test": [
            "coverage>=5.1",
            "coveralls>=2.0.0",
            "mock>=4.0.2",
            "pytest-cov>=2.10.0",
            "pytest-isort>=5",
            "pytest-mock>=3.1.1",
            "pytest>=6.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "blobtoolkit-host = lib.host:cli",
            "blobtoolkit-view = lib.view:cli",
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/blobtoolkit/blobtoolkit/issues",
        "Source": "https://github.com/blobtoolkit/blobtoolkit",
    },
    include_package_data=True,
    zip_safe=False,
)
