#!/bin/bash

#The following commands will create a symbolic link from the covid reference file to here
find ../d1_sars-cov-2_r94/ -type f -name "ref.fa" | xargs -i{} ln -s {} .
