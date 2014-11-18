#!/bin/bash

for file in alignments/sorted/*
do
    samtools index $file
done
