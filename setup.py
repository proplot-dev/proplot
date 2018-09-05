from setuptools import setup
setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name         = 'PubPlot',
    url          = 'https://github.com/lukelbd/pubplot',
    author       = 'Luke Davis',
    author_email = 'lukelbd@gmail.com',
    # Needed to actually package something
    packages = ['pubplot'],
    # Needed for dependencies
    install_requires = ['numpy', 'matplotlib'],
    # *Strongly* suggested for sharing
    version = '0.1',
    # The license can be anything you like
    license          = open('LICENSE').read(),
    description      = 'Tricked out matplotlib wrapper for making clear, compact, publication-quality graphics quickly and easily.',
    long_description = open('README.md').read(),
)
