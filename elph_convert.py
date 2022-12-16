import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import linecache as lc
import os

os.system("grep -A 1  '\-  1'  elph_anaddb.abo | cut -c 5- >  ph_eigs_atom1.dat")
os.system("grep -A 1  '\-  2'  elph_anaddb.abo | cut -c 5- >  ph_eigs_atom2.dat")
os.system("grep -A 6 ' Mode number    Frequency' elph_anaddb.abo | tail -6 > phonon_freqs.dat")

os.system("sed -r -i '/^\s*$/d'  ph_eigs_atom1.dat")
os.system("sed -r -i '/^\s*$/d'  ph_eigs_atom2.dat")

#file_name = input("Enter filename: ")
file_name = "yambo_elphon_gkk_bymode"

n_atm = int(lc.getline(file_name, 1).split()[-1])
n_ph = int(lc.getline(file_name, 2).split()[-1])
n_bnd = int(lc.getline(file_name, 3).split()[-1])
n_kpt = int(lc.getline(file_name, 4).split()[-1])
phonon_freq = np.multiply ( np.loadtxt('phonon_freqs.dat',usecols=1) , np.loadtxt('phonon_freqs.dat',usecols=1) )

print("Number of atoms: "+str(n_atm))
print("Number of phonon modes: "+str(n_ph))
print("Number of electronic bands: "+str(n_bnd))
print("Number of k-points: "+str(n_kpt))

pol_vec = np.zeros((3,n_atm,n_ph,2))

for i in range(n_atm):
    string_file = "ph_eigs_atom"+str(i+1)+".dat"
    vec_x = np.loadtxt(string_file,usecols=0)
    vec_y = np.loadtxt(string_file,usecols=1)
    vec_z = np.loadtxt(string_file,usecols=2)
    for j in range(n_ph): 
        pol_vec [:,i,j,0]  = np.array([vec_x[j*2],vec_y[j*2],vec_z[j*2]])
        pol_vec [:,i,j,1]  = np.array([vec_x[j*2],vec_y[j*2],vec_z[j*2]])

elph = np.zeros((n_kpt,n_bnd,n_bnd,n_ph,2))
sum_elph = np.zeros((n_kpt,n_bnd,n_bnd,n_ph))


for k in range(n_kpt):
    for b1 in range(n_bnd):
        for b2 in range(n_bnd):
            for p in range(n_ph):
                line_index = (4 + n_kpt + 3) + (k*(n_bnd*n_bnd*(n_ph+1)+1)) + (b1*(n_bnd*(n_ph+1))+1) + b2*(n_ph+1) + (p+1)   
                imag = float(lc.getline(file_name,line_index).split()[-1])*1.0
                real = float(lc.getline(file_name,line_index).split()[-2])*1.0
#                print(str(line_index)+" "+str(real)+" "+str(imag))
                elph[k,b1,b2,p,0]=real
                elph[k,b1,b2,p,1]=imag
                sum_elph[k,b1,b2,p] = real ** 2 + imag **2 


#print(np.sum(sum_elph[0,3,3,:]))
print(np.sum(sum_elph[0,1:,1:,:]))
rootgrp = nc.Dataset("ndb.elph_gkkp_expanded_fragment_1", "w", format="NETCDF4")

kpoints = rootgrp.createDimension("kpoints", n_kpt)
bands = rootgrp.createDimension("bands", n_bnd)
phonons = rootgrp.createDimension("phonons", n_ph)
cmplx = rootgrp.createDimension("cmplx", 2)
atoms = rootgrp.createDimension("atoms", n_atm)
axes = rootgrp.createDimension("axes", 3)

phonon_freqs = rootgrp.createVariable("PH_FREQS1","f8",("phonons"))
pol_vecs = rootgrp.createVariable("POLARIZATION_VECTORS","f8",("axes","atoms","phonons","cmplx"))
gkk_matrix = rootgrp.createVariable("ELPH_GKKP_Q1","f8",("kpoints","bands","bands","phonons","cmplx"))

gkk_matrix[:,:,:,:,:] = elph
pol_vecs[:,:,:,:] = pol_vec 
phonon_freqs[:] = phonon_freq
rootgrp.close()

rootgrp1 = nc.Dataset("ndb.elph_gkkp_expanded", "w", format="NETCDF4")
nq = rootgrp1.createDimension("num_q", 5)
pars = rootgrp1.createVariable("PARS","i8",("num_q"))
pars[:] = np.array([n_ph,n_kpt,n_kpt,n_bnd,1])

rootgrp1.close()

