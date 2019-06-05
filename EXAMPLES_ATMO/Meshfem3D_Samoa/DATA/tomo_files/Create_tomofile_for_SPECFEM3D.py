#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 12:41:31 2019

depth (km) vp vs rho Qp Qs (units: km/s; km/s; kg/cm3)
30.00 0.330 0.000 0.0013 9000.0 9000.0
0.000 1.480 0.000 1.0280 9000.0 9000.0
-2.00 7.500 4.300 3.2000 9000.0 9000.0

@author: carene
"""
import numpy as np
#import matplotlib.pyplot as plt

Ztop = 30000.0
Zbottom = -16500.0
zspacing = -100.0

Model = [[30.000,0.330,0.000,0.0013,9000.0,9000.0],
         [0.0000,1.480,0.000,1.0280,9000.0,9000.0],
         [-2.000,7.500,4.300,3.2000,9000.0,9000.0]]

#Dep = np.append(0,np.array([layer[0] for layer in Model]))
Dep = np.array([layer[0] for layer in Model])

vp = np.array([layer[1] for layer in Model])
vs = np.array([layer[2] for layer in Model])
rho = np.array([layer[3] for layer in Model])
Qp = np.array([layer[4] for layer in Model])
Qs = np.array([layer[5] for layer in Model])
#plt.step(vp,Dep,where='post',label='vp')
#plt.step(vs,Dep,where='post',label='vs')
#plt.step(rho,Dep,where='post',label='rho')

#plt.ylim(0,2.0)
#plt.ylim(plt.ylim()[::-1])

#plt.legend()
#plt.xlabel('vel (km/s); density (g/cm3)')
#plt.ylabel('depth (km)')
#plt.savefig('toto.pdf')

fout = open('tomography_model.xyz','w')
orig_x = -100000.0
orig_y = -100000.0
orig_z = Ztop 
end_x = 100000.0
end_y = 100000.0
end_z = Zbottom
fout.write(str(orig_x)+' '+str(orig_y)+' '+str(orig_z)+' '+str(end_x)+' '+str(end_y)+' '+str(end_z)+'\n')

spacing_x = 200000.0 
spacing_y = 200000.0 
spacing_z = zspacing
fout.write(str(spacing_x)+' '+str(spacing_y)+' '+str(spacing_z)+'\n')

nx = 1
ny = 1
nz = np.int((Zbottom-Ztop)/spacing_z) 
fout.write(str(nx)+' '+str(ny)+' '+str(nz)+'\n')

vpmin = np.min(vp)*1000.0
vpmax = np.max(vp)*1000.0
vsmin = np.min(vs)*1000.0
vsmax = np.max(vs)*1000.0
rhomin = np.min(rho)*1000.0
rhomax = np.max(rho)*1000.0

fout.write(str(vpmin)+' '+str(vpmax)+' '+str(vsmin)+' '+str(vsmax)+' '+str(rhomin)+' '+str(rhomax)+'\n')


id = 0 
print(Dep)
for d in np.arange(Ztop,Zbottom,spacing_z):
    if (id < 2 and d/1000.0<=Dep[id+1] ): id+=1
    print(d,id)
    fout.write('0.0 0.0 '+str(d)+' '+str(vp[id]*1000.0)+' '+str(vs[id]*1000.0)+' '+str(rho[id]*1000.0)
    +' '+str(Qp[id])+' '+str(Qs[id])+'\n')
fout.close()
