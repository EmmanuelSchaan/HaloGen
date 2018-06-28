import numpy as np
import idlsave


data = idlsave.read("./Bethermin2012model_grids.save")


# band passes for each experiment
band = data.get('band')
#array(['MIPS24', 'PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350',
#       'Planck857GHz', 'SPIRE500', 'Planck545GHz', 'SCUBA850',
#       'Planck353GHz', 'AzTEC', 'Planck217GHz', 'Planck143GHz',
#       'Planck100GHz', 'Radio'], dtype=object)

# corresponding central frequency
Nu = data.get('lambda') * 1.e-6 # wavelength, converted from microns to meters
Nu = 3.e8 / Nu # frequency in Hz


# counts per unit flux and redshift shell, for each frequency band
z = data.get('z')
Snu = data.get('snu') # in Jy
dNdSnudz = data.get('dndsnudz')  # counts of all sources (in gal/Jy/sr, Nband * Nredshift * Nflux array)

# remove nans
dNdSnudz = np.nan_to_num(dNdSnudz)

# save to separate files
np.savetxt("./converted/z.txt", z)
np.savetxt("./converted/Snu.txt", Snu)
np.savetxt("./converted/dNdSnudz_Planck100GHz.txt", dNdSnudz[14,:,:])
np.savetxt("./converted/dNdSnudz_Planck143GHz.txt", dNdSnudz[13,:,:])
np.savetxt("./converted/dNdSnudz_Planck217GHz.txt", dNdSnudz[12,:,:])
np.savetxt("./converted/dNdSnudz_Planck353GHz.txt", dNdSnudz[10,:,:])
np.savetxt("./converted/dNdSnudz_Planck545GHz.txt", dNdSnudz[8,:,:])
np.savetxt("./converted/dNdSnudz_Planck857GHz.txt", dNdSnudz[6,:,:])







