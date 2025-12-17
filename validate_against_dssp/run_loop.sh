#!/bin/bash

for f in test_cif/*.cif; do
  base=$(basename "$f" .cif)
  mkdssp --output-format 'dssp' "$f" "dssp_res_cif/${base}.dssp"
  echo "$base"
  echo ----
done

