domain=Antarctica
grid_name_src=ANT-16KM
grid_name_tgt=ANT-32KM
nc_src=../ice_data/${domain}/${grid_name_src}/${grid_name_src}_REGIONS.nc 

cdo gencon,grid_${grid_name_tgt}.txt -setgrid,grid_${grid_name_src}.txt ${nc_src} scrip-con_${grid_name_src}_${grid_name_tgt}.nc


# To perform remapping using the weights file
# cdo remap,grid_${grid_name_tgt}.txt,scrip-con_${grid_name_src}_${grid_name_tgt}.nc ${infile} ${outfile}

