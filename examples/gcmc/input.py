
#######################################################################################
#    [M] = a.m.u
#    [L] = Angs
#    [U] = [ML^2/T^2] = eV
#
#    therefore
#    [T] = sqrt(a.m.u/eV)*Angs
#
#
#    a.m.u = 1.66053904020e-27 Kg
#    Angs = 1.0e-10 m
#    eV = 1.60217656535e-19 Kg m^2/s^2
#    therefore
#    [T] = sqrt(1.66053904020/1.60217656535)*1.0e-14 s
#
#    kB = 8.617330350e-5 eV/K
#    h = 4.13566766225e-15 eVs
#    h = 4.13566766225 * 0.1 * sqrt(1.60217656535/1.66053904020)  sqrt(eV*a.m.u)*Angs
#######################################################################################

import time 
import os
import subprocess
proc = subprocess.Popen('rm -rf dumps/*', shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)

from math import sqrt
boltz = 8.617330350e-5
planck = 4.13566766225 * 0.1 * sqrt(1.60217656535/1.66053904020)



import mapp
from mapp import md
mapp.pause_slave_out();

################################################################

min=md.min_cg(e_tol=1.0e-8);
min.ntally=1;
min.ls=mapp.ls_brent();

################################################################

muvt=md.muvt(-2.33,300,0.1,'H',2411895,'Fe',-1.2,0.3);
muvt.nevery=10;
muvt.nattempts=10;
muvt.ntally=1;
muvt.export=md.export_cfg('dumps/dump',100)

################################################################



sim=md.atoms.import_cfg("configs/Fe.cfg")
sim.add_elem('H',1.007940)
sim.ff_eam_fs("potentials/FeH.eam.fs")
sim.hP=planck
sim.kB=boltz


min.run(sim,2)

#min.H_dof=[[True],[False,False],[False,False,True]]
#min.affine=True

#min.run(sim,500)

sim.create_temp(300.0,8569643);



sim.step=0

start = time.time()
muvt.run(sim,100);
print "time elapsed: %lf seconds" % (time.time()-start)

