#### Proprietary to LANL collaborators ####
###Carene Larmat, LANL, January 3 2015

###Alex had installed cubit on the clusters here are the instructions to load the software with the command module

##instructions 
On the cluster, you may want to use cubit without GUI (Otherwise, it is very slow). In this way, we can generate some large meshes by using large memory on cluster.

1. In your personal directory (e.g. my is /users/zhoul), create a folder “modulefiles”
    cd
    mkdir modulefiles
    cd modulefiles

2. In �modulefiles, create a file cubit with something looks like:

#%Module1.0
set prefix      /turquoise/usr/projects/spewave/zlei/SofterPackages/cubit
set bindir      $prefix/bin;
set mandir      $prefix/bin/help;
set libdir      $prefix/bin/lib;

prepend-path   PATH              $bindir;
prepend-path   MANPATH           $mandir;
prepend-path   LD_LIBRARY_PATH   $libdir;

##############################################################

3. Add the following lines in your “.cshrc” file which locates in your personal directory.

 cd
 nedit .cshrc
 source .cshrc

##############################################################

setenv MODULEPATH /users/zhoul/modulefiles/:$MODULEPATH

module load cubit
alias cubit cubit -batch -nographics

##############################################################

4. use the cubit to generate the mesh

cubit test.jou


#### Carene 

NB1: make sure to have the line #%Module1.0 as first line of the script
NB2: in my case, I have a .tcshrc instead of a .cshrc 

INSTR2: Most of the scripts to generate meshes are in python not in jou so we need to use claro test.py instead of cubit test.jou as done by Alex; we need to create an alias for claro as well, add the following line to .tcshrc

alias claro claro -batch -nographics 

INSTR3: Some of my personal scripts can be run with python instead of with claro; python needs to know where to find the cubit library for this add the following line to the .tcshrc

setenv PYTHONPATH /turquoise/usr/projects/spewave/zlei/SofterPackages/cubit/bin/

INSTR4: Reading the README.md of GEOCUBIT 
Requirements met on the clusters: cubit installed by Alex; python and numpy provided by the default configuration; instructions for paths are given for batch, in the following I will continue with instructions for tcshr instead. Add the following to .tcshrc replacing the setenv PYTHONPATH of previous instruction (still good but it will be defined in the following) 

first let's modify modulefiles/cubit to define CUBITDIR when calling module load cubit; add the following after after the line setting prefix
setenv CUBITDIR $prefix

At that point; we want to add the following to .tcshr (make sure this is after the line "module load cubit")
setenv CUBITLIB $CUBITDIR/bin:$CUBITDIR/structure
setenv PYTHONPATH $PYTHONPATH:$CUBITDIR/bin:$CUBITDIR/structure
Notice that there is no need to define LD_LIBRARY_PATH nor PATH as this is done when calling "module load cubit"
NB1: the line setenv PYTHONPATH can lead to an error if the variable is not defined. In that case, we need to test the definition of the variable:
if (! $?PYTHONPATH ) then
    setenv PYTHONPATH ${CUBITDIR}/bin:${CUBITDIR}/structure
else 
    setenv PYTHONPATH ${PYTHONPATH}:${CUBITDIR}/bin:${CUBITDIR}/structure
endif


If you don't have a build in your directory GEOCUBIT, do the following 
python setup.py install

INSTR5: Setting up the path for GEOCUBIT; I define a module geocubit that will set up the paths as setpaths.sh
Add the following line to the ~/.tcshrc to have geocubit set up automatically
module load geocubit


