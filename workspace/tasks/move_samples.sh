#!/bin/bash
# Moves downloaded dataset samples from finished packets to the Kingston disk (the root disk is tight).
cd ; D=external_media
for d in results/BT-DAT* results/BT-DENT-DAT* results/BT-DATX* results/BT-DENT-DATX*; do
  [ -f $d/RESULTS.md ] || continue
  for s in $d/samples $d/data $d/downloads $d/raw; do
    [ -d $s ] && [ ! -L $s ] || continue
    mkdir -p $D/$(basename $d) && rsync -a $s $D/$(basename $d)/ && rm -rf $s && echo "samples moved to $D/$(basename $d)/$(basename $s)" >> $d/SAMPLES_MOVED.txt
  done
done
