# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 11:10:45 2023

@author: Simon Ternes
"""

from supersatrnc import *;
from supersatrnc.materials import InertGas,SolventGas;
import pylab as plt;
import numpy as np;

allCompounds = [MAPI,TWOCAT,DMF,DMSO,GBL,DRY_AIR,N2,CB]
T = 293
p = 101325


materials = list();
liquids = list()
gases = list();

for i,c in enumerate(allCompounds):
    if(not c.sld == None):
       materials.append((c.sld,'Solid')) 
    if(not c.liq == None):
       materials.append((c.liq,'Liquid'))  
       liquids.append(c.liq)
    if(not c.gas == None):
       materials.append((c.gas,'Gas'))   
       gases.append(c.gas)
 
def test_material_properties():
    
    for i,ma in enumerate(materials):  
        m=ma[0];
        s=ma[1]
        print("-----"+m.ident_string+" "+s+"--------")
        
        if(s=="Gas"):
            assert m.massDens_kgpm3(T, p) > 0.1 and m.massDens_kgpm3(T, p) < 10
        else:
            assert m.massDens_kgpm3(T, p) > 100 and m.massDens_kgpm3(T, p) < 100000
        
        
        print("Mass Density ",m.massDens_kgpm3(T, p)/m._molMass_kgpmol, " kgpm3");
        
        print("Mol Density ",m.molDens_molpm3(T,p), " molpm3");
        
        assert np.abs(m.massDens_kgpm3(T, p)/m._molMass_kgpmol - m.molDens_molpm3(T,p)) < 1e-5
            
        
        print("Mol Density ",m.molDens_molpm3(T,p), " molpm3");
        
        assert np.abs(1/m.molDens_molpm3(T,p) - m.specVol_m3pmol(T,p)) < 1e-5
        
        print("Spec. Volume ",m.specVol_m3pmol(T,p), " m3pmol");
        
        
        if(s=='Liquid' or isinstance(m,InertGas)):
            
            print("kin. Visc ",m.kinViscosity_m2ps(T,p), " m2ps");
            
            assert np.abs(m.kinViscosity_m2ps(T,p) - m.dynViscosity_Pas(T,p)/m.massDens_kgpm3(T, p)) < 1e-5
            
            print("dyn. Visc ",m.dynViscosity_Pas(T,p), " Pas");
            
            if(s=='Liquid'):
                assert m.dynViscosity_Pas(T,p) < 0.1 and m.dynViscosity_Pas(T,p) > 0.0001
            else:
                assert m.dynViscosity_Pas(T,p) < 1e-3 and m.dynViscosity_Pas(T,p) > 1e-5
            
            print("Reynolds at 10 m/s and 10cm:",m.reynolds(T,p,10,0.1), "");
            
            assert m.reynolds(T,p,10,0.1) > 1e4
            
        if(isinstance(m,SolventGas)):
            print("partial Pressure:",m.vaporPressure_Pa(T), " Pa")
            
            assert m.vaporPressure_Pa(T) > 10 and m.vaporPressure_Pa(T) < 10000
    

def test_gaseous_diffusion():
    for i,g1 in enumerate(gases):
        for j,g2 in enumerate(gases):
            if(isinstance(g1,InertGas)):
                assert g1.diff_m2ps(T, p,g2) < 1e-4 and g1.diff_m2ps(T, p,g2) > 1e-7
                print("Diffusion of ",g2 ," in ",g1,":",g1.diff_m2ps(T, p,g2)," mps2")
                assert g1.schmidt(T, p,g2) > 0.1 and g1.schmidt(T, p,g2) < 3; 
                print("Schmidt :",g1.schmidt(T, p,g2))

def test_liquid_diffusion():
    for i,l1 in enumerate(liquids):
        for j,l2 in enumerate(liquids):
            assert l1.diff_m2ps(T, p,l2) < 1e-8 and l1.diff_m2ps(T, p,l2) > 1e-10
            print("Diffusion of ",l2 ," in ",l1,":",l1.diff_m2ps(T, p,l2)," mps2")
            assert l1.schmidt(T, p,l2) > 100 and l1.schmidt(T, p,l2) < 10000; 
            print("Schmidt :",l1.schmidt(T, p,l2))
        

def test_plot():
    plt.figure()
    
    
    corr = LaminarAverage(DRY_AIR.gas,32,5,32,0.75)
    sp = SimulationParameters(corr,20);
    film1 = Film(sp,[DMF,MAPI],4.7, [1], 0.8)
    dynsDynGas1 = GasQuenchDynamics(film1,55);
    dynsDynGas1.computeSolution();
    plt.plot(dynsDynGas1.times,dynsDynGas1.result.thicknesses_m())
    plt.show()
    assert dynsDynGas1.result.thicknesses_m()[0] < 100e-6 and dynsDynGas1.result.thicknesses_m()[0] > 100e-9


test_plot()