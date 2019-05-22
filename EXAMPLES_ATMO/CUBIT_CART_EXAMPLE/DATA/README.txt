awk ' BEGIN {OFS="\t";xs=395776.696100823581;ys=8327079.706598646939} {print $1,$2,$4-ys,$3-xs,"0.00","0.00"} ' output_list_stations.txt
