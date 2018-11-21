#!/bin/bash
#------------------------------------------------------------------------------#
# ***Notes on getting ttf files on Mac****
# /System/Library/Font *OR* /Library/Fonts
# * The .otf files work in addition to .ttf files; you can verify this by
#   looking at plot.fonts_files_os (results of findSystemFonts command) --
#   the list will include a bunch of .otf files.
# * Some system fonts are .dfont which are unreadable to matplotlib; use
#   fondu (download with Homebrew) to break down.
# * Sometimes fondu created .bdf files, not .ttf; use https://github.com/Tblue/mkttf
#   Requires FontForge and PoTrace
# * Install new fonts with: "brew cask install font-<name-of-font>" after using
#   "brew tap caskroom/fonts" to initialize; these will appear in ~/Library/Fonts.
# * For some conversion tools, run "pip install fonttools". Documentation
#   is here: https://github.com/fonttools/fonttools
#------------------------------------------------------------------------------#
# ***Notes on default files packaged in font directory.
# * Location will be something like: /lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf
# * 'STIX' fonts allow different LaTeX-like math modes e.g. blackboard bold
#   and caligraphy; see: https://matplotlib.org/gallery/text_labels_and_annotations/stix_fonts_demo.html
# * The 'cm'-prefix fonts seem to provide additional mathematical symbols
#   like integrals, and italized math-mode fonts.
# * We also have 'pdfcorefonts' in this directory, but I think since these
#   are afm matplotlib cannot use them? Don't know.
#------------------------------------------------------------------------------#
# This function sets up any .ttf fonts contained in the <fonts> directory to
# be detected by matplotlib, regardless of the OS.
# See: https://olgabotvinnik.com/blog/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/
# Best fonts for dyslexia: http://dyslexiahelp.umich.edu/sites/default/files/good_fonts_for_dyslexia_study.pdf
#------------------------------------------------------------------------------#
# First add the fonts
# For abspath function see: https://stackoverflow.com/a/21188136/4970632
function abspath() { # abspath that works on mac, Linux, or anything with bash
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}
dir=$(abspath ${0%/*}) # this script's directory
base=$dir/../panplot/fonts
mpldir=$(python -c "import matplotlib; print(matplotlib.matplotlib_fname())")
mfontdir=${mpldir%matplotlibrc}fonts/ttf
echo "Transfering .ttf and .otf files to $mfontdir."
for font in $base/*.[ot]tf; do
  if [ ! -r "$mfontdir/${font##*/}" ]; then
    cp "$font" "$mfontdir/"
    echo "Added font \"${font##*/}\"."
  fi
done

# Then delete the cache; echo each one
# For get_cachedir see: https://stackoverflow.com/a/24196416/4970632
shopt -s nullglob # keep empty if does not exist
cachedir=$(python -c "import matplotlib; print(matplotlib.get_cachedir())")
rm $cachedir/font* 2>/dev/null
caches=($cachedir/*.cache) # may be more than one
for cache in ${caches[@]}; do # might be more than one
  if [ ! -d "$cache" ]; then # may be a tex.cache folder
    rm "$cache"
    echo "Deleted font cache \"$cache\"."
  fi
done
