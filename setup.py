
from setuptools import setup, find_packages
import bioframes


setup(
    name = bioframes.__projectname__,
    version = bioframes.__release__,
    packages = find_packages(),
    author = bioframes.__authors__,
    author_email = bioframes.__authoremails__,
    description = bioframes.__description__,
    license = "GPLv2",
    keywords = bioframes.__keywords__,
    install_requires = [
       'pandas',
       'numpy',
       'schema',
       'pyparsing',
       'biopython',
       'pyvcf',
       'toolz',
       'sh',
       'mock'
    ],
)

