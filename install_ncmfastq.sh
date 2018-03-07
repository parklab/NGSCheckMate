#!/bin/bash

# ngscheckmate_fastq
cd ngscheckmate_fastq-source && make && chmod +x ngscheckmate_fastq && cd ..
cp ngscheckmate_fastq-source/ngscheckmate_fastq .

# patterngenerator
cd patterngenerator && make && cd ..

