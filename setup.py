import os
from setuptools import setup, Command, find_packages

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "uniquesam",
    version = "0.1.0",
    author = "dlmeduLi",
    author_email = "dlmeduLi@163.com",
    description = ("Analyse sam file and keep the unique aligment record."),
    license = "BSD",
    keywords = "bioinformatic samfile unique",
    url = "https://github.com/dlmeduLi/uniqe-sam",
    packages=find_packages(),
    install_requires=['libsam'],
    long_description=read('README.rst'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
    ],
    entry_points = {
        'console_scripts': ['unique-sam=unique_sam.unique_sam:main'],
    },
    # ... Other setup options
    cmdclass={
        'clean': CleanCommand,
    }
)