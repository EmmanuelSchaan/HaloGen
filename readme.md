# HaloGen for Line Intensity Mapping

Modular halo model code for 3d power spectra, and the corresponding projected 2d power spectra in the Limber and flat sky approximations.
Observables: matter density, galaxy lensing, CMB lensing, thermal Sunyaev-Zel'dovich (Hill & Pajer 2013), cosmic infrared background (Penin+12, 14), tracers with any dn/dz, b(z) and HOD.
Computes all auto and cross spectra. Computes halo model trispectrum in simple configurations.
Mass functions and halo biases from Press & Schechter, Sheth & Tormen and Tinker are available.

Requires pathos for multiprocessing (https://pypi.org/project/pathos/): 
```
pip install pathos
```
Just clone and run: 
```
python driver.py
```
Hope you find this code useful! Please cite https://arxiv.org/abs/1406.3330 and/or https://arxiv.org/abs/1802.05706 if you use this code in a publication. Do not hesitate to contact me with any question: eschaan@lbl.gov
