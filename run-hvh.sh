#!/bin/bash
mkdir -p build/tables
ERROR=0.1
printf "%s\n" {2..5} |  xargs -P 2 -I {} het-vs-hom/hvh -e $ERROR -o build/tables/e$ERROR-d{}.tsv {}
