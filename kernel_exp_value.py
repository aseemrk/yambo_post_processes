# Read BSE kernel and BSE eigenstates to calculate exp. value of the kernel
# One can calculate the exchange contribution by doing a separate calculation
# with BSE kernel set to hartree


# Warning: don't use symmetries, use parallel HDF5 io in Yambo

import numpy as np
import netCDF4 as nc

bse_kernel_file = 'ndb.BS_PAR_Q1' # Read kernel file 
bse_wfn_file = 'ndb.BS_diago_Q1'  # Read BSE wavefunction file

bse_wfn_db = nc.Dataset(bse_wfn_file,'r')
bse_kernel_db = nc.Dataset(bse_kernel_file,'r')

bse_wfn = bse_wfn_db['BS_EIGENSTATES'][...].data  # This gets the data as stored, [...] is for all
bse_kernel = bse_kernel_db['BSE_RESONANT'][...].data

bse_kernel = bse_kernel[...,0] + 1j*bse_kernel[...,1] # It is a complex number but stored in array in netcdf
bse_wfn = bse_wfn[...,0] + 1j*bse_wfn[...,1]

# Yambo does something weird:
# It stores only part of the data i.e. upper triangular matrix 
# but it stores it in the lower triangular part 

bse_kernel = bse_kernel.T # bring lower triangular to upper triangular 

lower_idx = np.tril_indices(n =bse_kernel.shape[0],m =bse_kernel.shape[1]) # This gets all the lower tri. indices
bse_kernel[lower_idx[0],lower_idx[1]] = np.conj(bse_kernel[lower_idx[1],lower_idx[0]]) # Assigning them to conj of upper tri.

exp_val = np.einsum('si,ij,sj->s',np.conj(bse_wfn),bse_kernel,bse_wfn,optimize=True) * 27.2114 # einsum to get the exp. value
# Above, s: index of exciton, i and j are "k,c,v" indices, also we multiply with Hartree to eV constant
# Typicaly, the kcv order is not the same for eigenstate and the kernel 
# but we can get away without worrying about it here since we haven't used symmetries at all
np.savetxt('exp_val.dat',np.c_[np.real(exp_val)])



