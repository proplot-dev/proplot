#!/bin/bash
set -e
set -eo pipefail

echo 'Code Styling with (flake8, isort)'

source activate proplot-dev

echo '[flake8]'
flake8 proplot docs --max-line-length=88 --ignore=W503,E402,E741

# TODO: Add this
# echo '[black]'
# black --check -S proplot

echo '[isort]'
isort --recursive --check-only --line-width=88 proplot

echo '[doc8]'
doc8 *.rst docs --max-line-length=88
