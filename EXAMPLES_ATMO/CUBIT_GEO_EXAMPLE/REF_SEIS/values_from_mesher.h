
 !
 ! purely informative use
 !
 ! mesh statistics:
 ! ---------------
 !
 ! note: 
 !    the values are only approximate and differ for different processes
 !    because the CUBIT + SCOTCH mesh has
 !    a different number of mesh elements and points in each slice
 !
 ! number of processors =           32
 !
 ! number of ES nodes =    4.00000000    
 ! percentage of total 640 ES nodes =   0.625000000      %
 ! total memory available on these ES nodes (Gb) =    64.0000000    
 !
 ! min vector length =           25
 ! min critical vector length =           75
 !
 ! master process: total points per AB slice =        34761
 ! total elements per AB slice = (will be read in external file)
 ! total points per AB slice = (will be read in external file)
 !
 ! total for full mesh:
 ! -------------------
 !
 !
 ! number of time steps =        20000
 !
 ! time step =    1.0000000000000000E-002
 !
 ! attenuation uses:
 !  NSPEC_ATTENUATION =            1
 ! 
 ! anisotropy uses:
 !  NSPEC_ANISO =            1
 ! 
 ! adjoint uses:
 !  NSPEC_ADJOINT =            1
 !  NGLOB_ADJOINT =            1
 ! 
 ! approximate least memory needed by the solver:
 ! ----------------------------------------------
 !
 ! size of arrays for the largest slice =    7.3675155639648438       MB
 !                                      =    7.1948394179344177E-003  GB
 !
 !   (should be below 90% or so of the amount of memory available per processor core
 !   (if significantly more, the job will not run by lack of memory)
 !   (if significantly less, you waste a significant amount of memory)
 !
 ! check parameter to ensure the code has been compiled with the right values:
&MESHER
 ABSORB_FREE_SURFACE_VAL=T,
 /

