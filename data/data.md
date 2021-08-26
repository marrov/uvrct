# Data for the CSTR UV reactor simulation

Values were obtained using the Spectral Atlas:

Keller-Rudek, H., Moortgat, G. K., Sander, R., and Sörensen, R.: The MPI-Mainz UV/VIS spectral atlas of gaseous molecules of atmospheric interest, Earth Syst. Sci. Data, 5, 365–373, (2013), DOI: [10.5194/essd-5-365-2013](https://doi.org/10.5194/essd-5-365-2013)

Specific references are mentioned when more than one option was available. The references employed were selected based on recency and/or if the publication was a review paper.

The species selected were the ones that appear in the `photolysis.inp` file based on the most relevant species of the mechanism proposed by [Meusinger (2013)](https://doi.org/10.1016/j.cej.2016.08.092)

## Absorbtion cross section data

All data are obtained for the species in gas phase at ~300 K and for UV radiation at a radiation wavelength of 254 nm.

| Species | Absorbtion cross section [cm^2] | Reference |
|---------|---------------------------------|-----------|
| O3      | 1.132935E-17                    | [Hodges 2019](https://doi.org/10.1088/1681-7575/ab0bdd) |
| O2      | 1.477865E-24                    | [Bogumil 2003](https://doi.org/10.1016/S1010-6030(03)00062-5) |
| N2      | Not reported above 150 nm       | [None](http://satellite.mpic.de/spectral_atlas/cross_sections/Nitrogen+compounds(N,H,O)/N2.spc) |
| SO2     | 1.551258E-19                    | [Bogumil 2003](https://doi.org/10.1016/S1010-6030(03)00062-5) |
| H2O     | Approaches 0 (see this [figure](http://joseba.mpch-mainz.mpg.de/spectral_atlas_data/cross_sections_plots/Hydrogen+water/H2O_186-230nm_log.jpg))  | [Ranjan 2020](http://satellite.mpic.de/spectral_atlas/cross_sections/Hydrogen+water/H2O_Ranjan(2020)_292K_192.057-230.413nm(extrapolated).txt) |
| HO2     | 2.63E-19                        | [Tyndall 2001](https://doi.org/10.1029/2000JD900746) |
| SO3     | 1.34E-20                        | [Burkholder 1997](https://doi.org/10.1029/97GL03255) |
| H2SO4   | 7.19E-22 (at 248 nm)            | [Farahani 2019](https://doi.org/10.1021/acs.jpca.9b04668) |

## Quantum yield

For the reaction: 

`O3 + hv(254nm) => O(1D) + O2`

the photolysis effect must be splitted up into two reactions, one for each product. This is done in Cantera so that teh quantum yields of the individual specied can be specified. It reads as follows:

`O3 + hv(254nm) => O(1D)`

`O3 + hv(254nm) => O2`

The quantum yield for the first reaction is: **Phi_O(1D) = 0.9** for low temperatures (between 298 K and 312 K) at a radiation wavelength of 250 nm  as reported by [Matsumi 2002](https://doi.org/10.1029/2001JD000510). Note that the value is constant in the range 220 nm to 300 nm. This is also reported in [Atkinson 2004](https://doi.org/10.5194/acp-4-1461-2004) in page 1497.

The quantum yield for the second reaction is: **Phi_O2 = 0.1**, this is reported also both by [Matsumi 2002](https://doi.org/10.1029/2001JD000510) and by [Atkinson 2004](https://doi.org/10.5194/acp-4-1461-2004) in page 1497.
