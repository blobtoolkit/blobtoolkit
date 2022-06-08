"""
A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

import io

# import platform
# import re
from os import makedirs
from os import remove
from os.path import abspath
from os.path import dirname
from os.path import isfile
from os.path import join

from setuptools import find_packages
from setuptools import setup

# from setuptools.command.install import install

# system = platform.system()
# # arch, _ = platform.architecture()
# if system == "Linux":
#     blobtoolkit_api = "bin/blobtoolkit-api-linux"
#     blobtoolkit_viewer = "bin/blobtoolkit-viewer-linux"
# if system == "Windows":
#     blobtoolkit_api = "bin/blobtoolkit-api-win.exe"
#     blobtoolkit_viewer = "bin/blobtoolkit-viewer-win.exe"
# if system == "Darwin":
#     blobtoolkit_api = "bin/blobtoolkit-api-macos"
#     blobtoolkit_viewer = "bin/blobtoolkit-viewer-macos"

# from https://stackoverflow.com/questions/45150304/how-to-force-a-python-wheel-to-be-platform-specific-when-building-it # noqa
# try:
#     from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

#     class bdist_wheel(_bdist_wheel):
#         def finalize_options(self):
#             _bdist_wheel.finalize_options(self)
#             # Mark us as not a pure python package (we have platform specific rust code)
#             self.root_is_pure = False

#         def get_tag(self):
#             # this set's us up to build generic wheels.
#             # note: we're only doing this for windows right now (causes packaging issues
#             # with osx)
#             tag = _bdist_wheel.get_tag(self)
#             print(tag)
#             repl = "macosx_10_6_intel.macosx_10_9_intel.macosx_10_9_x86_64"
#             if tag[2] == "macosx_10_9_x86_64":
#                 tag = (tag[0], tag[1], repl)
#                 return tag

#             if not system == "Windows":
#                 return _bdist_wheel.get_tag(self)

#             python, abi, plat = _bdist_wheel.get_tag(self)
#             python, abi = "py3", "none"
#             return python, abi, plat


# except ImportError:
#     print("Error happens")
#     bdist_wheel = None


# class PostInstallCommand(install):
#     """Post-installation for installation mode."""

#     def run(self):
#         source_dir = dirname(abspath(__file__))
#         makedirs(self.install_scripts, exist_ok=True)

#         self.copy_file(source_dir, blobtoolkit_api, "blobtoolkit-api")
#         self.copy_file(source_dir, blobtoolkit_viewer, "blobtoolkit-viewer")

#     def copy_file(self, source_dir, file_name, executable_name):
#         source = join(source_dir, file_name)
#         target = join(self.install_scripts, executable_name)
#         if isfile(target):
#             remove(target)

#         self.move_file(source, target)


def read(*names, **kwargs):
    """Read file."""
    with io.open(
        join(dirname(__file__), *names), encoding=kwargs.get("encoding", "utf8")
    ) as fh:
        return fh.read()


setup(
    name="blobtoolkit",  # Required
    version="3.2.1",
    description="blobtoolkit",  # Optional
    long_description="blobtoolkit",  # Optional
    long_description_content_type="text/markdown",
    # long_description="%s\n%s"
    # % (
    #     re.compile("^.. start-badges.*^.. end-badges", re.M | re.S).sub(
    #         "", read("README.rst")
    #     ),
    #     re.sub(":[a-z]+:`~?(.*?)`", r"``\1``", read("CHANGELOG.rst")),
    # ),
    # long_description_content_type="text/x-rst",  # Optional (see note above)
    url="https://github.com/blobtoolkit/blobtools2",  # Optional
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
    #   py_modules=["my_module"],
    #
    packages=find_packages(where="src"),  # Required
    # Specify which Python versions you support. In contrast to the
    # 'Programming Language' classifiers above, 'pip install' will check this
    # and refuse to install the project if the version does not match. See
    # https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
    python_requires=">=3.6, <4",
    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        "chromedriver-binary-auto==0.1.1",
        "docopt>=0.6.2",
        "fastjsonschema==2.15.3",
        "geckodriver-autoinstaller==0.1.0",
        "psutil==5.9.0",
        "pysam==0.18.0",
        "pyvirtualdisplay==2.2",
        "pyyaml",
        "selenium==4.1.0",
        "tolkein>=0.4.0",
        "tqdm==4.62.3",
        "ujson>=3.0.0",
        # Additional pipeline dependencies
        "defusedxml==0.7.1",
        "requests==2.27.1",
        "snakemake==6.15.3",
    ],  # Optional
    # List additional groups of dependencies here (e.g. development
    # dependencies). Users will be able to install these using the "extras"
    # syntax, for example:
    #
    #   $ pip install sampleproject[dev]
    #
    # Similar to `install_requires` above, these must be valid existing
    # projects.
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
        "console_scripts": ["blobtools = blobtools:cli", "btk = btk:cli",],
        "blobtools.subcmd": [
            "add = blobtools.lib.add:cli",
            "create = blobtools.lib.add:cli",
            "filter = blobtools.lib.filter:cli",
            "host = blobtools.lib.host:cli",
            "remove = blobtools.lib.remove:cli",
            "replace = blobtools.lib.add:cli",
            "validate = blobtools.lib.validate:cli",
            "view = blobtools.lib.view:cli",
        ],
        "btk.subcmd": [
            "pipeline = btk_pipeline:main",
            # "blobtools = blobtools:cli",
        ],
        "pipeline.subcmd": [
            "data = btk_pipeline.lib.data:main",
            "run = btk_pipeline.lib.run:main",
            "add-summary-to-metadata = btk_pipeline.lib.add_summary_to_metadata:main",
            "chunk-fasta = btk_pipeline.lib.chunk_fasta:main",
            "count-busco-genes = btk_pipeline.lib.count_busco_genes:main",
            "extract-busco-genes = btk_pipeline.lib.extract_busco_genes:main",
            "generate-config = btk_pipeline.lib.generate_config:main",
            "generate-static-images = btk_pipeline.lib.generate_static_images:main",
            "resume-pipeline = btk_pipeline.lib.resume_pipeline:main",
            "transfer-completed = btk_pipeline.lib.transfer_completed:main",
            "unchunk-blast = btk_pipeline.lib.unchunk_blast:main",
            "window-stats = btk_pipeline.lib.window_stats:main",
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/blobtoolkit/blobtools2/issues",
        "Source": "https://github.com/blobtoolkit/blobtools2",
    },
    include_package_data=True,
    zip_safe=False,
)
