#!python
#!python
#!/usr/bin/env python
# February 2014 Carene & Lucie 
# first attempt mesh Samoa area local 
# Revised January 2016 as examples ATMOFORK and modifs for GEOCUBIT library of SPECFEM3D v3.0 

import cubit
import boundary_definition
import cubit2specfem3d

import os
import sys

#for utm example
# moves volume to UTM coordinates of topography surface
#cubit.cmd('volume 1 move x 339101 y  -1714072 z -54000')


########################## USER SETTINGS #######################
##Parameter models to be tuned 
size_side_model = 50000 #in m 
height_solid_earth = 15000 #in m 
height_ocean = 2000 # in m 
height_atmo = 30000 # in m 
elementsize = 2000.0 
#################################################################

#starting from scratch 
cubit.cmd('reset')

# three volumes solid earth; ocean; atmosphere in that order
# volumes are moved so (0,0,0) is in the middle (hor) and at sea level 

##solid earth 
cubit.cmd('brick x '+str(size_side_model)+' y '+str(size_side_model)+' z '+str(height_solid_earth))
cubit.cmd('volume 1 move  z -'+str(height_ocean+height_solid_earth/2))

##ocean 
cubit.cmd('brick x '+str(size_side_model)+' y '+str(size_side_model)+' z '+str(height_ocean))
cubit.cmd('volume 2 move z -'+str(height_ocean/2))

##atmosphere
cubit.cmd('brick x '+str(size_side_model)+' y '+str(size_side_model)+' z '+str(height_atmo))
cubit.cmd('volume 3 move z '+str(height_atmo/2))

cubit.cmd('imprint volume all')
cubit.cmd('merge all')
cubit.cmd('compress all')

# Meshing the volumes
##Meshing 1 simple
cubit.cmd('volume all size '+str(elementsize))
cubit.cmd('mesh volume all')

##More complex - need to have the numbering of surface right... 
#cubit.cmd('surface 13 scheme equal size '+str(elementsize)+' scheme map')
#cubit.cmd('mesh surface 13')
#cubit.cmd('surface 2 scheme equal size '+str(elementsize)+' scheme map')
#cubit.cmd('mesh surface 2')
#cubit.cmd('volume 1 2 3 scheme sweep vector 0 0 1 size '+str(elementsize))
#cubit.cmd('mesh volume all')

###deleting meshes
#cubit.cmd('delete mesh volume all propagate')

#### End of meshing
# optional saves backups
cubit.cmd('export mesh "samoa.e" dimension 3 overwrite')
cubit.cmd('save as "samoa.cub" overwrite')

##restart from cub file 
#cubit.cmd('import cubit "samoa.cub" merge_globally')
##restart from e file
#cubit.cmd('import mesh geometry "samoa.e"')
###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

#### Define material properties for the 3 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')

#Get numbering of blocks
blocks=cubit.get_block_id_list()

## "Attributes" homogeneous solid Earth 
id_block = blocks[0]
cubit.cmd('block '+str(id_block)+' name "elastic 1" ')        # material region
cubit.cmd('block '+str(id_block)+' attribute count 6')
cubit.cmd('block '+str(id_block)+' attribute index 1 1  ')      # volume 1
cubit.cmd('block '+str(id_block)+' attribute index 2 7500 ')
cubit.cmd('block '+str(id_block)+' attribute index 3 4300 ')
cubit.cmd('block '+str(id_block)+' attribute index 4 3200 ')
cubit.cmd('block '+str(id_block)+' attribute index 5 9000.0 ')
cubit.cmd('block '+str(id_block)+' attribute index 6 0 ')     # anisotropy_flag
# prem with tomography model 
#cubit.cmd('block '+str(id_block)+' name "elastic prem.xyz 1" ')        # elastic material region
#cubit.cmd('block '+str(id_block)+' attribute count 2')
#cubit.cmd('block '+str(id_block)+' attribute index 1 -1')      # flag for material: -1 for 1. undefined material
#cubit.cmd('block '+str(id_block)+' attribute index 2 2')      # flag for tomographic model 

## "Attributes" ocean 
id_block = blocks[1] 
cubit.cmd('block '+str(id_block)+' name "acoustic 1" ')        # material region
cubit.cmd('block '+str(id_block)+' attribute count 4')
cubit.cmd('block '+str(id_block)+' attribute index 1 2  ')      # volume 2
cubit.cmd('block '+str(id_block)+' attribute index 2 1480 ')   # vp
cubit.cmd('block '+str(id_block)+' attribute index 3 0 ')   # vs
cubit.cmd('block '+str(id_block)+' attribute index 4 1028 ')   # rho

## "Attributes" atmosphere
id_block = blocks[2]
#homogeneous
cubit.cmd('block '+str(id_block)+' name "acoustic 2" ')        # material region
cubit.cmd('block '+str(id_block)+' attribute count 4')
cubit.cmd('block '+str(id_block)+' attribute index 1 3  ')      # volume 3
cubit.cmd('block '+str(id_block)+' attribute index 2 330 ')   # vp
cubit.cmd('block '+str(id_block)+' attribute index 3 0 ')   # vs
cubit.cmd('block '+str(id_block)+' attribute index 4 1.300 ')   # rho
#profile
#cubit.cmd('block '+str(id_block)+' name "acoustic atmo_model_SAM_4km.xyz 1" ')        # acoustic material region
#cubit.cmd('block '+str(id_block)+' attribute count 2')
#cubit.cmd('block '+str(id_block)+' attribute index 1 -1')      # flag for material: -1 for 1. undefined material
#cubit.cmd('block '+str(id_block)+' attribute index 2 2')      # flag for tomographic model 


# optional saves backups
cubit.cmd('export mesh "mesh_final.e" dimension 3 overwrite')

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SPECFEM3D('MESH')

# all files needed by SCOTCH are now in directory MESH





