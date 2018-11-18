from setuptools import setup
# For including non-python data, see:
# https://stackoverflow.com/a/1857436/4970632
setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name         = 'ProPlot',
    url          = 'https://github.com/lukelbd/explot',
    author       = 'Luke Davis',
    author_email = 'lukelbd@gmail.com',
    # Package stuff
    # Also include package data
    packages     = ['explot'],
    package_data = {'': ['cmaps/*', 'fonts/*', 'colors/*']},
    # Command-line scripts
    # scripts = ['scripts/explot_fonts'],
    # Needed for dependencies
    install_requires = ['numpy', 'matplotlib'],
    # *Strongly* suggested for sharing
    version = '1.0',
    # The license can be anything you like
    license          = open('LICENSE.txt').read(),
    description      = 'Matplotlib wrapper for making clear, concise, publication-quality graphics quickly and easily.',
    long_description = open('README.md').read(),
)
