#!python
#!python
#!/usr/bin/env python
# February 2014 Carene & Lucie 
# first attempt mesh Samoa area local 

import cubit
import boundary_definition
import cubit2specfem3d

import os
import sys


##restart from cub file 
cubit.cmd('import mesh geometry "samoa.e"')
###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

#### Define material properties for the 3 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')

## "Attributes" homogeneous solid Earth 
cubit.cmd('block 1 name "elastic 1" ')        # material region
cubit.cmd('block 1 attribute count 6')
cubit.cmd('block 1 attribute index 1 1  ')      # volume 1
cubit.cmd('block 1 attribute index 2 7500 ')
cubit.cmd('block 1 attribute index 3 4300 ')
cubit.cmd('block 1 attribute index 4 3200 ')
cubit.cmd('block 1 attribute index 5 9000.0 ')
cubit.cmd('block 1 attribute index 6 0 ')     # anisotropy_flag
# prem with tomography model 
#cubit.cmd('block 1 name "elastic prem.xyz 1" ')        # elastic material region
#cubit.cmd('block 1 attribute count 2')
#cubit.cmd('block 1 attribute index 1 -1')      # flag for material: -1 for 1. undefined material
#cubit.cmd('block 1 attribute index 2 2')      # flag for tomographic model 

## "Attributes" ocean 
cubit.cmd('block 2 name "acoustic 1" ')        # material region
cubit.cmd('block 2 attribute count 4')
cubit.cmd('block 2 attribute index 1 2  ')      # volume 2
cubit.cmd('block 2 attribute index 2 1480 ')   # vp
cubit.cmd('block 2 attribute index 3 0 ')   # vs
cubit.cmd('block 2 attribute index 4 1028 ')   # rho

## "Attributes" atmosphere
#homogeneous
#cubit.cmd('block 3 name "acoustic 2" ')        # material region
#cubit.cmd('block 3 attribute count 4')
#cubit.cmd('block 3 attribute index 1 3  ')      # volume 3
#cubit.cmd('block 3 attribute index 2 330 ')   # vp
#cubit.cmd('block 3 attribute index 3 0 ')   # vs
#cubit.cmd('block 3 attribute index 4 1.300 ')   # rho
#profile
cubit.cmd('block 3 name "acoustic atmo_model_SAM_4km.xyz 1" ')        # acoustic material region
cubit.cmd('block 3 attribute count 2')
cubit.cmd('block 3 attribute index 1 -1')      # flag for material: -1 for 1. undefined material
cubit.cmd('block 3 attribute index 2 2')      # flag for tomographic model 


# optional saves backups
cubit.cmd('export mesh "top.e" dimension 3 overwrite')
cubit.cmd('save as "meshing.cub" overwrite')

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SESAME('MESH')

# all files needed by SCOTCH are now in directory MESH





