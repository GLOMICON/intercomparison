#!/bin/bash

source ${CONFIG}

ASV_TABLE_FILTERED="outputs/asv_table/"$PROJECT"_dada2_"$VERSION".filtered.table.tsv.gz"

files_to_clean="$ASV_TABLE_FILTERED"

for f in $files_to_clean; do
  cp "$f" "$f~" &&   
  gzip -cd "$f~" | sed 's/NNNNNNNNNN//g' | gzip > "$f"
  mv -- "$f" "${f%.tsv.gz}.tsv_without_N.gz"
  mv "$f~" "$f"
done
