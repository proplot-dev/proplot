#!/usr/bin/env bash
#------------------------------------------------------------------------------#
# Function for converting notebooks into rst files
# This is in my .bashrc, but copied it to this repo so collaborators can use it
#------------------------------------------------------------------------------#
# Find notebooks
cwd=$(pwd)
notebooks=($cwd/quickstart.ipynb) # bash array of notebooks you want to convert
cd "${0%/*}"/../docs # go to docs of script
for notebook in "${notebooks[@]}"; do
  # Convert notebook
  base=${notebook%.ipynb}
  base=${base##*/}
  jupyter nbconvert --to=rst --output-dir=. $notebook

  # Change default nbconvert file location
  rm -r $base 2>/dev/null # remove old one
  mv ${base}_files $base
  gsed -i "s:${base}_files:${base}:g" ${base}.rst

  # Run special magical vim regexs that add back in sphinx
  # links and default non-figure cell output. Your only other option for
  # non-greedy regexes is perl, and fuck that. Also need to match naked
  # module (e.g. ``numpy``) literals, hence the search for specific packages.
  echo "Running vi search and replace."
  command vi -c \
     '%s/``\(\~\?\)\(datetime\|cycler\|mpl_toolkits\|numpy\|pandas\|climpy\|metpy\|scipy\|xarray\|matplotlib\|cartopy\|proplot\)\(.\{-}\)``/`\1\2\3`/ge | '"\
    "'%s/:ref:``\(.\{-}\)``/:ref:`\1`/ge | '"\
    "'%s/code::\s*$/code:: ipython3/ge | '"\
    "'%s/.. parsed-literal::\n\n.\+\_.\{-}\n\n//ge | wq' ${base}.rst &>/dev/null

  # Split notebook into files based on sections with magical awk command
  # For each line, awk runs bracket command if search is valid, if condition is
  # empty, or if variable evaluates to true (i.e. if non-empty and non-zero)
  cat ${base}.rst | awk '
  /^=/ {++count; file="'$base'"count".rst"; print file}
  file {print line > file}
  {line=$0}
  ' -
  rm ${base}.rst
done
