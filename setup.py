"""
this setup will check for dependencies and install pyaxon on your computer
"""
from setuptools import setup, find_packages

setup(
    name='pyaxon',
    version='1.1.1',
    url='https://github.com/mooniean/axonalGrowth.git',
    author='mooniean and phydev',
    author_email='-',
    description='Multi-phase-field model for axonal growth and mRNA transport',
    license='GNU GPLv3',
    platform='Python 3.7',
    packages=find_packages(),
    install_requires=['numpy >= 1.14.3',
                      'scipy >= 1.3.0'],
