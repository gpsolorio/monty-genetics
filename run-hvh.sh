#!/bin/bash
mkdir -p tables
ITERATIONS=1e9
ERROR=0.1
printf "%s\n" {2..20} |  xargs -P 8 -I {} python3 hvh.py $ITERATIONS {} -o tables/e$ERROR-d{}.tsv