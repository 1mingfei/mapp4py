
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

muvt=md.muvt(-2.8,300,0.1,'Ag',2411895,'Au',-2.8,0.5);
muvt.nevery=5;
muvt.nattempts=10;
muvt.ntally=1;
muvt.export=md.export_cfg('dumps/dump',1)

################################################################



sim=md.atoms.import_cfg("configs/CuAgAu.cfg")
sim.add_elem('Ag',107.8682)
sim.add_elem('Au',196.967)
sim.ff_eam_setfl("potentials/CuAgAu_Zhou04.eam.alloy")
sim.hP=planck
sim.kB=boltz


min.run(sim,2)

#min.H_dof=[[True],[False,False],[False,False,True]]
#min.affine=True

#min.run(sim,500)

sim.create_temp(300.0,8569643);



sim.step=0

start = time.time()
muvt.run(sim,1000);
print "time elapsed: %lf seconds" % (time.time()-start)

