from setuptools import setup
# For including non-python data, see:
# https://stackoverflow.com/a/1857436/4970632
# NOTE: To rename the repo (which you do a lot), use this command:
# find . \( -name '*.ipynb' -o -name '*.py' -o -name '*.md' -o -name '*.txt' -o -name '.vimsession' \) -exec gsed -i 's/panplot/proplot/g' {} +
setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name         = 'proplot',
    url          = 'https://github.com/lukelbd/proplot',
    author       = 'Luke Davis',
    author_email = 'lukelbd@gmail.com',
    # Package stuff
    # Also include package data
    packages     = ['proplot'], # reserve name
    package_data = {'': ['cmaps/*', 'fonts/*', 'colors/*', '.proplotrc']},
    # Needed for dependencies
    install_requires = ['numpy', 'matplotlib'],
    # *Strongly* suggested for sharing
    version = '1.0',
    # The license can be anything you like
    license          = open('LICENSE.txt').read(),
    description      = 'Matplotlib wrapper for making clear, concise, publication-quality graphics quickly and easily.',
    long_description = open('README.rst').read(),
)
