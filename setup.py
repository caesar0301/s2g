import os
import sys
from setuptools import setup

__version__ = '0.2.5.1'

setup(
    name = "s2g",
    version = __version__,
    url = 'https://github.com/caesar0301/python-s2g',
    author = 'Xiaming Chen',
    author_email = 'chenxm35@gmail.com',
    description = 'Shapefile to graph/network converter in Python',
    long_description='''''',
    license = "MIT",
    packages = ['s2g'],
    keywords = ['shapefile', 'graph', 'network', 'GIS'],
    install_requires=[
        'fiona',
        'shapely',
        'networkx',
        'numpy',
        'progressbar2'
      ],
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: Freely Distributable',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2.6',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.2',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Topic :: Software Development :: Libraries :: Python Modules',
   ],
)
