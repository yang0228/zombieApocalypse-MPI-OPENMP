#!/bin/bash
mpic++ -fopenmp stage2.cpp -o stage2 && sbatch run.sbatch