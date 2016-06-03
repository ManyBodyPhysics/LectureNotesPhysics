#!/bin/bash
pdflatex lnp.tex
for auxfile in chap*.aux
do
    bibtex `basename $auxfile .aux`
done
pdflatex lnp.tex
pdflatex lnp.tex