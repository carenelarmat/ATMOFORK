#!/bin/tcsh -f 
awk ' {print $1,$2,$3+8673261.33304,$4+391135.82663,$5,$6} ' STATIONS_start > toto

rm -rf STATIONS
set nl = `wc -l toto | awk ' {print $1} ' `
set il = 1
while ( $il <= $nl )
 echo $il
 set yx = `awk ' NR=='$il' {print $3,$4} ' toto `
 echo 1 >! in.in
 echo $yx[2] >> in.in
 echo $yx[1] >> in.in
 echo -2 >> in.in
 ../a.out < in.in > out
 set longlat = `awk ' /long,lat/ {print $2,$3} ' out `
 set stanet = `awk ' NR=='$il' {print $1,$2} ' STATIONS_start`
 set dbur = `awk ' NR=='$il' {print $5,$6} ' STATIONS_start`
 echo $stanet $longlat[2] $longlat[1] $dbur >> STATIONS

@ il += 1
end
