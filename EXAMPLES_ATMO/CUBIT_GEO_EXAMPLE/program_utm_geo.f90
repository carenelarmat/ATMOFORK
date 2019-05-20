program call_utm_geo
implicit none

double precision :: long,lat,x,y
integer :: UTM_PROJECTION_ZONE,iway


write(*,*) ' Enter 0 for long-lat to UTM; 1 for UTM to long-lat '
read(*,*) iway

if (iway == 0) then
 write(*,*) 'Enter long,lat'
 read(*,*) long,lat
else
 write(*,*) 'Enter x,y'
 read(*,*) x,y
endif

write(*,*) 'Enter UTM_PROJECTION_ZONE'
read(*,*) UTM_PROJECTION_ZONE

call utm_geo(long,lat,x,y,UTM_PROJECTION_ZONE,iway,.false.)

if (iway == 0) then
 write(*,*) 'UTM coordinates x,y ',x,y
else
 write(*,*) 'long,lat ',long,lat
endif

end program call_utm_geo
