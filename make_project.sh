#!/bin/bash
# Bash script for making the project
pdflatex -shell-escape project.tex
biber project
pdflatex -shell-escape project.tex
#bibtex project