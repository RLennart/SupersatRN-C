# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 12:31:44 2023

@author: Simon Ternes

    Supersaturation Rate Numerical Calculator (material properties) - Calculates the dynamic thickness and supersaturation of perovskite solution films 
    
    Copyright (C) 2023  Simon Ternes

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
"""


#Commonly used constants
R_JpmolK = 8.3144598;
NA_1pmol = 6.02214 * 10**23;


#Commonly used constants
R_JpmolK = 8.3144598;
NA_1pmol = 6.02214 * 10**23;


solvents = dict();



solvents['DMSO'] = {'massDens_gpcm3':1.1004, #https://en.wikipedia.org/wiki/Dimethyl_sulfoxide
                    'molMass_gpmol':78.133, #http://webbook.nist.gov/cgi/cbook.cgi?ID=C67685&Mask=4&Type=ANTOINE&Plot=on
                    'diffVol':(15.9 * 2 + 2.31 * 6 + 6.11 + 22.9), #2*C + 6*H+1*O+1*S WLI, Da28
                    'antoineParams':[[[293-30,323],[5.23039,2239.161,-29.215],True],[[323.5,442],[4.49107,1807.002,-60.995],True]], #[0]:tempRange, [1]:A,B,C [2]unit is bar (mmHg) http://webbook.nist.gov/cgi/cbook.cgi?ID=C67685&Mask=4&Type=ANTOINE&Plot=on
                    'surfaceTensionParams':[43.54,-0.1145], #https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
                    'dynViscosity_mPas':2.271, #https://mdpi-res.com/d_attachment/applsci/applsci-12-00116/article_deploy/applsci-12-00116-v2.pdf?version=1640332022
                    'isAntiSolvent':False
                    }

solvents['DMF'] = {'massDens_gpcm3':0.948 , #https://en.wikipedia.org/wiki/Dimethylformamide
                    'molMass_gpmol':73.0938, #http://webbook.nist.gov/cgi/cbook.cgi?ID=C68122&Mask=4
                    'diffVol':(15.9 * 3 + 2.31 * 7 + 6.11 + 4.54 ), #3*C + 7*H+1*O+1*N WLI, Da28
                    'antoineParams':[[[303-30,363],[3.93068,1337.716,-82.648],True]], #http://webbook.nist.gov/cgi/cbook.cgi?ID=C68122&Mask=4
                    'surfaceTensionParams':[37.1,-0.14],  #https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
                    'dynViscosity_mPas':0.92,#https://en.wikipedia.org/wiki/Dimethylformamide
                    'isAntiSolvent':False
                    }

solvents['GBL'] = {'massDens_gpcm3':1.1296 , #https://en.wikipedia.org/wiki/Gamma-Butyrolactone
                    'molMass_gpmol':86.090, #https://en.wikipedia.org/wiki/Gamma-Butyrolactone
                    'diffVol':(15.9 * 4 + 2.31 * 6 + 6.11*2-18.3), #4*C + 6*H+2*O-1*Aromatic Ring, WLI, Da28
                    'antoineParams':[[[-43,465.85],[7.67415,2147.61,244.041],False]], #http://webbook.nist.gov/cgi/cbook.cgi?ID=C68122&Mask=4
                    'surfaceTensionParams':[37.1,-0.14], #https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
                    'dynViscosity_mPas':1.7, #https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/celc.201700138
                    'isAntiSolvent':False
                    }

solvents['2ME'] = {'massDens_gpcm3':0.965 , #https://en.wikipedia.org/wiki/2-Methoxyethanol
                    'molMass_gpmol':76.09, #https://en.wikipedia.org/wiki/2-Methoxyethanol
                    'diffVol':(15.9 * 3 + 2.31 * 8 + 6.11*2), #2*C + 8*H+2*O
                    'antoineParams':[[[-43,465.85],[5.06386,1853.556,-30.838],True]], #https://webbook.nist.gov/cgi/cbook.cgi?ID=C109864&Mask=4&Type=ANTOINE&Plot=on
                    'surfaceTensionParams':[30.84,-0.08], #https://materials.springer.com/lb/docs/sm_lbs_978-3-540-75508-1_53
                    'dynViscosity_mPas':1.7, #https://www.carlroth.com/medias/SDB-8892-AU-EN.pdf?context=bWFzdGVyfHNlY3VyaXR5RGF0YXNoZWV0c3wyNjg5MTB8YXBwbGljYXRpb24vcGRmfGg2OC9oMGMvOTE0ODEzMjQyNTc1OC9TREJfODg5Ml9BVV9FTi5wZGZ8NGM5OGRiMzJhYmJlNzhjNGFiYTU0ZjRkZmI3NTI3YjBjYmIzNzg0YzRmY2E2YmJlMDMyMzlhZTFkYmNlNmI3MA
                    'isAntiSolvent':False
                    }

solvents['CB'] = {'massDens_gpcm3': 1.11 , #https://en.wikipedia.org/wiki/Chlorobenzene
                    'molMass_gpmol': 	112.56, #https://en.wikipedia.org/wiki/Gamma-Butyrolactone
                    'diffVol':(15.9 * 6+2.31 * 5 + 21*1-18.3), #6*C + 5*H+1*Cl-1*Aromatic Ring, WLI, Da28
                    'antoineParams':[[[335-60,405],[4.11083,1435.675,-55.124],True]], #https://vacuu-lan.com/vapor-pressure-estimates-antoine-equation/
                    'surfaceTensionParams':[33.6,-0.1191], #https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
                    'dynViscosity_mPas':0.753, #https://en.wikipedia.org/wiki/Chlorobenzene_(data_page)
                    'isAntiSolvent':True
                    }


solvents['VARCB'] = {'massDens_gpcm3': 1.11 , #https://en.wikipedia.org/wiki/Chlorobenzene
                    'molMass_gpmol': 	112.56, #https://en.wikipedia.org/wiki/Gamma-Butyrolactone
                    'diffVol':(15.9 * 6+2.31 * 5 + 21*1-18.3), #6*C + 5*H+1*Cl-1*Aromatic Ring, WLI, Da28
                    'antoineParams':[[[335-60,405],[4.11083,1435.675,-55.124],True]], #https://vacuu-lan.com/vapor-pressure-estimates-antoine-equation/
                    'surfaceTensionParams':[33.6,-0.1191], #https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
                    'dynViscosity_mPas':0.753, #https://en.wikipedia.org/wiki/Chlorobenzene_(data_page)
                    'isAntiSolvent':True
                    }

dryingGases = dict(); #'RspecGas, viscGas, diffVolGas, mMassGas, thermalConduct, specifc heat

dryingGases['DRY_AIR'] = {'dynViscSutherParams':[120,291.15,18.27*10**(-6)], #https://de.wikipedia.org/wiki/Sutherland-Modell
                          'diffVol':19.7, # Da28
                          'molMass_gpmol':28.9647, #https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
                          'thermalConductivity': 26.2, #https://de.wikipedia.org/wiki/W%C3%A4rmeleitf%C3%A4higkeit
                          'specHeatCap': 1006 #https://www.engineeringtoolbox.com/air-properties-d_156.html
    }# g/mol

dryingGases['N2'] = {'dynViscSutherParams':[111,300.55,18.27*10**(-6)], #https://de.wikipedia.org/wiki/Sutherland-Modell
                     'diffVol':18.5, # Da28
                     'molMass_gpmol':14, #https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
                     'thermalConductivity': 26.0, #https://de.wikipedia.org/wiki/W%C3%A4rmeleitf%C3%A4higkeit
                     'specHeatCap': 1040 #https://www.nuclear-power.com/nitrogen-specific-heat-latent-heat-vaporization-fusion/
    }

crystals = dict(); # Temperatures, Cubic lattice constant in A

#https://www.nature.com/articles/srep35685/figures/6
crystals['MAPI'] = {'latticeConstant_A':6.25, #https://materials.hybrid3.duke.edu/materials/24
                    'molMass_gpmol': (67.52 + 207.2 + 126.9*3),
                    'eqConcentration_molpm3':{'DMF':1520,'DMSO':1520,'GBL':1520}, #Si section 1
                    'critConcentration_molpm3':{'DMF':6340,'DMSO':6340,'GBL':6340}, #Si section 1
                    'quenchingEff_molpm3':{'CB':{'DMF':2.88220474*1000,'DMSO':2.88220474*1000,'GBL':2.88220474*1000},'VARCB':{'DMF':2.88220474*1000,'DMSO':2.88220474*1000,'GBL':2.88220474*1000}}} #Si section

crystals['TWOCAT'] = {'latticeConstant_A':6.36,#https://pubs.acs.org/doi/pdf/10.1021/acs.jpclett.5b01432 for FAPI
                      'molMass_gpmol':(132.91*0.17 + 171.97*0.83 + 207.2 + 126.9*3*0.91 + 79.9*3*0.09),
                      'eqConcentration_molpm3':{'DMF':1520,'DMSO':1520,'GBL':1520}, #Si section 1
                      'critConcentration_molpm3':{'DMF':6016,'DMSO':6016,'GBL':6016}, #Si section 1
                      'quenchingEff_molpm3':{'CB':{'DMF':3.64457431*1000,'DMSO':3.64457431*1000,'GBL':3.64457431*1000},'VARCB':{'DMF':2.88220474*1000,'DMSO':2.88220474*1000,'GBL':2.88220474*1000}}} #Si section
