# SL3ME
A Python3 implementation of the Spectroscopic Limited Maximum Efficiency (SLME) analysis of solar absorbers.

Originally written and maintained by Logan Williams. Original code cleaned and made more efficient by Michael Waters prior to upload.



############### CITATIONS AND FURTHER REFERENCES #########################
The SLME method was created by Liping Yu and Alex Zunger. Please cite their work:
L. Yu, A. Zunger, Phys. Rev. Lett. 108, 068701 (2012).
https://doi.org/10.1103/PhysRevLett.108.068701

For more analysis of the method, see the book chapter by M. Bercx, R. Saniz, B. Partoens, D. Lamoen called "Exceeding the Shockley–Queisser Limit Within the Detailed Balance Framework":
Bercx M., Saniz R., Partoens B., Lamoen D. (2018) Exceeding the Shockley–Queisser Limit Within the Detailed Balance Framework. In: Angilella G., Amovilli C. (eds) Many-body Approaches at Different Scales. Springer, Cham
https://doi.org/10.1007/978-3-319-72374-7_15

also at arxiv:
https://arxiv.org/pdf/1705.07762.pdf

################ USAGE INSTRUCTIONS #####################################
The file SL3ME.py is a python3 library containing a single function called "calculate_SLME". It can also be called as a script to plot an example using mock CdTe absorbance data.

The calculate_SLME function expects a file in the directory called "am1.5G.dat". This is the AM1.5G solar spectrum, found at:
https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls

By accessing the AM1.5G solar spectrum file, you agree to the NREL data disclaimer:
https://www.nrel.gov/disclaimer.html

ASTM solar spectrum general info page and references:
https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html

To make the am1.5G.dat file, copy the "Wvlgth nm" and "Global tilt  W*m-2*nm-1" columns into a blank text file. You can leave the headers in place (but if removed, missing the first 2 lines misses very little of the solar spectrum).


The function expects input in the form of two lists: one a list of energy values in eV, the other a list of absorbances (in m^-1 or cm^-1 with a flag) at the energy values of the first list, and two values: the smallest direct allowed band gap of the material, and the smallest indirect band gap of the material.
See the __main__ script part of the library for an example on use when you have a data file of energies and absorbances.

It also can calculate at different temperatures and different absorber material thicknesses. If the library is called as a script, it will plot a SLME vs thickness plot as an example.
