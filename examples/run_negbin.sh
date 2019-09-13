#!/bin/bash

NITS=$1
SEED=$2

if [ -z $NITS ]; then
    echo "./run_negbin.sh NumIters";
    exit;
fi

if [ -z $SEED ]; then
    SEED=1;
fi

mkdir -p nbin

./sim_negbin.R    $SEED

## simple examples
./negbin.R    $NITS 0.25 $SEED 1

## Rubin Gelman statistics for examples
./negbin.R    $NITS 0.5 $SEED 30

## perfect simulation examples
./negbin.R    $NITS 0.25 $SEED 0
