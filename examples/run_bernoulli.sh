#!/bin/bash

NITS=$1
SEED=$2

if [ -z $NITS ]; then
    echo "./run_bernoulli.sh NumIters";
    exit;
fi

if [ -z $SEED ]; then
    SEED=1;
fi

mkdir -p bern

./sim_bernoulli.R $SEED

## simple examples
./bernoulli.R $NITS 0.5 $SEED 1

## Rubin Gelman statistics for examples
./bernoulli.R $NITS 0.5 $SEED 30

## perfect simulation examples
./bernoulli.R $NITS 0.5 $SEED 0
