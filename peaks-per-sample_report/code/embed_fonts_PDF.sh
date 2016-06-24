#!/bin/bash

## USAGE: embed_fonts_PDF.sh /path/to/file.pdf
## this script will embed fonts into a PDF document
## this requires the fonts to be installed on your computer, and requires ghostscript to be installed and in your PATH

# ~~~~~~ script args processing ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 1)); then
  grep '^##' $0
  exit
fi

# ~~~~~~ script args ~~~~~~ #
# get the script arguments here
PDF_input="$1"; echo -e "PDF_input is $PDF_input"

echo -e "Embedding PDF fonts with ghost script:"
gs -dCompatibilityLevel=1.4 -dPDFSETTINGS=/screen -dCompressFonts=true -dSubsetFonts=true -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${PDF_input%%.pdf}_embedfonts.pdf -c ".setpdfwrite <</NeverEmbed [ ]>> setdistillerparams" -f "$PDF_input"
