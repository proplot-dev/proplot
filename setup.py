# Notes
# * To rename the repo, use this command:
#   find . \( -name '*.ipynb' -o -name '*.py' -o -name '*.md' -o -name '*.txt' -o -name '.vimsession' \) -exec gsed -i 's/panplot/proplot/g' {} +
# * install_requires spec is run on pip install, requirements.txt can be
#   run manually with pip install -r requirements.txt but is not run automatically.
# * For including non-python data, see:
#   https://stackoverflow.com/a/1857436/4970632
from setuptools import setup
setup(
    name = 'proplot',
    url = 'https://lukelbd.github.io/proplot',
    author = 'Luke Davis',
    version = '1.0',
    author_email = 'lukelbd@gmail.com',
    python_requires = '>=3.6.0',
    project_urls={
        'Bug Tracker': 'https://github.com/lukelbd/proplot/issues',
        'Documentation':  'https://lukelbd.github.io/proplot',
        'Source Code': 'https://github.com/lukelbd/proplot'
        },
    packages = ['proplot'], # reserve name
    package_data = {'': ['cmaps/*', 'fonts/*', 'colors/*', '.proplotrc']},
    install_requires = ['matplotlib>=2.2', 'numpy>=1.11', 'lxml>=4.0.0'],
    license = open('LICENSE.txt').read(),
    description = 'Matplotlib wrapper for making beautiful, publication-quality graphics.',
    long_description = open('README.rst').read(),
)
