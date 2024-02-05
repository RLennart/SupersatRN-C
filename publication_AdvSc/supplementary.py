# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 10:48:26 2023

@author: Simon Ternes
"""

from supersatrnc import *;
import pylab as plt;
import style as st;
import numpy as np;
import os;


abs_path = os.path.dirname(os.path.abspath(__file__));

corr = LaminarAverage(DRY_AIR.gas,32,5,32,1)



sp = SimulationParameters(corr,22);

film1 = Film(sp,[DMSO,DMF,TWOCAT],5.0, [0.3,0.7], 0.86)

dynsVac1 = VacuumDynamics(film1,25);
dynsVac2 = VacuumDynamics(film1,25,activities={DMF:1,DMSO:0.5});

sq = SimulationQueue([["d_mu",1.2]], [dynsVac1,dynsVac2])


ac = ActivityFitter(sq,abs_path+'/../dyn_interf_measurements/Double_Cation_DMF_DMSO.npy',fixedParams=[DMF])

res = ac.startFitting();

print("Fitting Result",res)

fig,axs = plt.subplots(nrows=2,figsize=st.figsize())  

sp = SimulationParameters(corr,22);

film1 = Film(sp,[DMSO,DMF,TWOCAT],5.0, [0.3,0.7], 0.86)

dynsVac1 = VacuumDynamics(film1,25);
dynsVac2 = VacuumDynamics(film1,25,activities={DMF:1,DMSO:res[0][0]});

sq = SimulationQueue([["d_mu",1.2]], [dynsVac1,dynsVac2])
sq.startQueue();
sq.plotAllThicknessCompositionAndSupersaturation(axs,[abs_path+'/../dyn_interf_measurements/Double_Cation_DMF_DMSO.npy'],['Interf. Meas.']);

plt.savefig('suppl_Vacuum_DMF_DMSO.png')


corr = LaminarAverage(DRY_AIR.gas,32,5,32,1)
sp = SimulationParameters(corr,22);

film1 = Film(sp,[GBL,DMSO,DMF,TWOCAT],7.3, [1/3*1.3,0.2*2/3*1.3,0.8*2/3*0.7], 0.86)
dynsVac1 = VacuumDynamics(film1,25);
dynsVac2 = VacuumDynamics(film1,25,activities={DMF:1,DMSO:res[0][0],GBL:0.5}); #
sq = SimulationQueue([["d_mu",1.2]], [dynsVac1,dynsVac2])

ac = ActivityFitter(sq,abs_path+'/../dyn_interf_measurements/Double_Cation.npy',fixedParams=[DMF,DMSO])

res2 = ac.startFitting();

print("Fitting Result 2",res2)

dynsVac1 = VacuumDynamics(film1,25);
dynsVac2 = VacuumDynamics(film1,25,activities={DMF:1,DMSO:res[0][0],GBL:res2[0][0]});
sq = SimulationQueue([["d_mu",1.2]], [dynsVac1,dynsVac2])

fig,axs = plt.subplots(nrows=2,figsize=st.figsize())  


sq.startQueue();
sq.plotAllThicknessCompositionAndSupersaturation(axs,[abs_path+'/../dyn_interf_measurements/Double_Cation.npy'],['Interf. Meas.']);



plt.savefig('2Cat_Gas.png')




fig,axs = plt.subplots(nrows=2,figsize=st.figsize())  
corr = LaminarAverage(DRY_AIR.gas,32,5,32,1)
sp = SimulationParameters(corr,22);

film1 = Film(sp,[DMF,MAPI],3.435+0.9819, [1], 1.18)

dynsVac = VacuumDynamics(film1,12.5);

dynsVac.computeSolution();


tc = DryAndWetThicknessFitter(dynsVac, abs_path+'/../dyn_interf_measurements/MAPI2.npy',[0,20],[1])

res = tc.startFitting();

print(res)

film1 = Film(sp,[DMF,MAPI],(res[0][0]+res[0][1])*1000*1000, [1], res[0][1]*1000*1000)
dynsVac = VacuumDynamics(film1,12.5);
dynsVac.computeSolution();
dynsVac.plotFilmThicknessCompositionAndSupersaturation(axs,[abs_path+'/../dyn_interf_measurements/MAPI.npy'],['Interf. Meas.']);

data = np.load(abs_path+'/../dyn_interf_measurements/MAPI.npy');


cryst_thickness = dynsVac.result.thicknesses_m(offset=1.15250168e-06)[np.argmin(np.abs(dynsVac.times-7.4))]*1000*1000;

print(1.15e-6*MAPI.sld.molDens_molpm3(22+273.15, 1e5)/(cryst_thickness/1000/1000))
print(1.15e-6*TWOCAT.sld.molDens_molpm3(22+273.15, 1e5)/(cryst_thickness/1000/1000))

print(data[0][-1],np.argmin(np.abs(dynsVac.times-data[0][-1])),dynsVac.result.thicknesses_m(offset=0.735e-6)[np.argmin(np.abs(dynsVac.times-data[0][0]))]*1000*1000);

print(np.argmin(np.abs(dynsVac.times-data[0][-1])))
