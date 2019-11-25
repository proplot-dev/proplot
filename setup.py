from setuptools import setup
from os.path import exists

with open('requirements.txt') as f:
    install_req = [req.strip() for req in f.read().split('\n')]
install_req = [req for req in install_req if req and req[0] != '#']

classifiers = [
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3 :: Only',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    ]

if exists('README.rst'): # when does this not exist?
    with open('README.rst') as f:
        long_description = f.read()
else:
    long_description = ''

if exists('LICENSE.txt'):
    with open('LICENSE.txt') as f:
        license_text = f.read()
else:
    license_text = ''

setup(
    url = 'https://lukelbd.github.io/proplot',
    name = 'proplot',
    version = '1.0',
    author = 'Luke Davis',
    author_email = 'lukelbd@gmail.com',
    maintainer = 'Luke Davis',
    maintainer_email = 'lukelbd@gmail.com',
    python_requires = '>=3.6.0',
    project_urls={
        'Bug Tracker': 'https://github.com/lukelbd/proplot/issues',
        'Documentation':  'https://lukelbd.github.io/proplot',
        'Source Code': 'https://github.com/lukelbd/proplot'
        },
    packages = ['proplot'],
    package_data = {'': ['cmaps/*', 'fonts/*', 'colors/*', '.proplotrc']},
    include_package_data = True, # use MANIFEST.in
    classifiers = classifiers,
    install_requires = install_req,
    license = license_text,
    description = 'A comprehensive wrapper for making beautiful, publication-quality graphics.',
    long_description = long_description,
    )
