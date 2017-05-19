#!/bin/bash

## COINF binary model
cat $1 | sed "s/XYX//g" | /app/vw-8.2 -i hpv16.k18.s4000.coinf.binary.model -p /dev/stdout > coinf.preds.txt

## primary lineage model
cat $1 | sed "s/XYX//g" | /app/vw-8.2 -i hpv.k18.s4000.lineage.ect.model -p /dev/stdout > lineage.preds.txt
## primary sublineage model

cat $1 | sed "s/XYX//g" | /app/vw-8.2 -i hpv16.k18.s4000.sublineage.ect.model -p /dev/stdout > sublineage.preds.txt

