#!/bin/bash

INPUT="/home/agenor/splus/pids/PST_CONF_P01_reduced"
CONF="/home/agenor/splus/SITH/JEDIS.conf"
LOGDIR="./runlog"

mkdir -p "$LOGDIR"

for dir in "$INPUT"/SGR*; do
    field=$(basename "$dir")

    echo
    echo "=== Rodando $field ==="

    sith "$field" "$CONF" 2>&1 | tee "$LOGDIR/${field}.log"
done
