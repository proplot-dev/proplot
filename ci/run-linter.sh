#!/bin/bash
# Run the travis CI tests
# WARNING: Make sure to keep flags in sync with .pre-commit.config.yaml
set -e
set -eo pipefail

echo 'Code Styling with (flake8, isort)'

echo '[flake8]'
flake8 proplot docs --exclude .ipynb_checkpoints --max-line-length=88 --ignore=W503,E402,E741

echo '[isort]'
isort --recursive --check-only --line-width=88 --skip __init__.py --multi-line=3 --force-grid-wrap=0 --trailing-comma proplot

# echo '[black]'
# black --check -S proplot

# echo '[doc8]'
# doc8 ./*.rst docs/*.rst --ignore D001  # ignore line-too-long due to RST tables
