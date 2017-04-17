# setup.py: based off setup.py for toil-vg, modified to install this pipeline
# instead.
import sys
import os

# Get the local version.py and not any other version module
execfile(os.path.join(os.path.dirname(os.path.realpath(__file__)), "version.py"))

from setuptools import find_packages, setup
from setuptools.command.test import test as TestCommand

kwargs = dict(
    name='hgvm-builder',
    version=version,
    description="Human Genome Variation Map construction kit",
    author='Adam Novak',
    author_email='anovak@soe.ucsc.edu',
    url="https://github.com/BD2KGenomics/hgvm-builder",
    install_requires=[package + ver for package, ver in required_versions.iteritems()],
    dependency_links = dependency_links,
    tests_require=['pytest==2.8.3'],
    package_dir={'': 'src'},
    packages=find_packages('src'),
    entry_points={
        'console_scripts': ['build-hgvm = hgvmbuilder.build:entrypoint']})

class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        # Sanitize command line arguments to avoid confusing Toil code
        # attempting to parse them
        sys.argv[1:] = []
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

kwargs['cmdclass'] = {'test': PyTest}

setup(**kwargs)

# Wen we run setup, tell the user they need a good Toil with cloud support
print("""
Thank you for installing the hgvm-builder pipeline!

If you want to run this Toil-based pipeline on a cluster in a cloud, please
install Toil with the appropriate extras. For example, To install AWS/EC2
support for example, run

pip install toil[aws,mesos]{}

on every EC2 instance. For Microsoft Azure, deploy your cluster using the Toil
template at

https://github.com/BD2KGenomics/toil/tree/master/contrib/azure

For more information, please refer to Toil's documentation at

http://toil.readthedocs.io/en/latest/installation.html

To start building HGVMs, run

build-hgvm --help 2>&1 | less
""".format(required_versions['toil']))

