#!/bin/bash

NITS=$1

if [ -z $NITS ]; then
    echo "./run.sh NumIters";
    exit;
fi

bash run_bernoulli.sh $NITS $SEED
bash run_negbin.sh $NITS $SEED
