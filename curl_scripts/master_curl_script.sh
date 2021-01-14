#! /usr/bin/env bash

while read -r line; 
do 

  echo "Downloading..."
  curl --globoff -n -c ~/.urs_cookies -b ~/.urs_cookies -LJO --url $line
  

done < test.txt 

# subset_M2I1NXASM_5.12.4_20210114_185325.txt  subset_M2I3NPASM_5.12.4_20210114_185319.txt

#subset_M2I1NXASM_5.12.4_20210114_185325.txt 


# Example command is 
# curl --globoff -n -c ~/.urs_cookies -b ~/.urs_cookies -LJO --url 'https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NPASM.5.12.4/2019/05/MERRA2_400.inst3_3d_asm_Np.20190501.nc4.nc4?EPV[0:7][0:41][179:181][144:146],H[0:7][0:41][179:181][144:146],O3[0:7][0:41][179:181][144:146],OMEGA[0:7][0:41][179:181][144:146],PHIS[0:7][179:181][144:146],PS[0:7][179:181][144:146],QI[0:7][0:41][179:181][144:146],QL[0:7][0:41][179:181][144:146],QV[0:7][0:41][179:181][144:146],RH[0:7][0:41][179:181][144:146],SLP[0:7][179:181][144:146],T[0:7][0:41][179:181][144:146],U[0:7][0:41][179:181][144:146],V[0:7][0:41][179:181][144:146],time,lev,lat[179:181],lon[144:146]'


