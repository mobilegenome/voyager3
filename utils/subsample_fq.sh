#!/bin/bash
# subsample from fastq.gz file
zcat $1 | head -10000000 | gzip
