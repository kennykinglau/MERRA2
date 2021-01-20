#! /usr/bin/env bash

while read -r line; 
do 

  echo "Downloading..."
  curl --globoff -n -c ~/.urs_cookies -b ~/.urs_cookies -LJO --url $line
  

done < test.txt 

# subset_M2I1NXASM_5.12.4_20210119_234818.txt  subset_M2I3NPASM_5.12.4_20210119_235152.txt

# subset_M2I1NXASM_5.12.4_20210114_185325.txt  subset_M2I3NPASM_5.12.4_20210114_185319.txt
