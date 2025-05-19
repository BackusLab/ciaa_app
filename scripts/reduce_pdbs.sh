#!/usr/bin/env bash
mkdir -p data/reduced_files 

for i in data/pdb_files/*.pdb; 
do
  filename=$(basename "$i") 
  echo "$filename"
  grep -vwE "HELIX|SHEET|DT|DC|DG|DA" "$i" > "strip_$filename"
  grep -v "HETATM*" "strip_$filename" > "stripped_$filename"
  reduce -build "stripped_$filename" 2> reduce.log 1> "reduced_$filename"
  mv "reduced_$filename" data/reduced_files/
  rm "strip_$filename" "stripped_$filename" reduce.log
done