# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 10:48:26 2023

@author: cassi
"""

from supersatrnc import *;
from pylab import plt;
import style as st;
import os;


abs_path = os.path.dirname(os.path.abspath(__file__));

fig,axs = plt.subplots(nrows=2,figsize=st.figsize())  
corr = LaminarAverage(DRY_AIR.gas,32,5,32,1)
sp = SimulationParameters(corr,22);

film1 = Film(sp,[DMF,MAPI],5.3, [1], 1.18)

dynsVac = VacuumDynamics(film1,12.5);

dynsVac.computeSolution();

dynsVac.plotFilmThicknessCompositionAndSupersaturation(axs,[abs_path+'/../dyn_interf_measurements/MAPI.npy'],['Interf. Meas.']);

plt.savefig('MAPI_VACUUM.png')

corr = LaminarAverage(DRY_AIR.gas,32,5,32,1)
fig,axs = plt.subplots(nrows=2,figsize=st.figsize())  

sp = SimulationParameters(corr,22);

film1 = Film(sp,[GBL,DMSO,DMF,TWOCAT],7.4, [1/3*1.3,0.2*2/3*1.3,0.8*2/3*0.7], 0.86)

dynsVac1 = VacuumDynamics(film1,25);
dynsVac2 = VacuumDynamics(film1,25,activities={DMF:1,DMSO:0.16479906,GBL:0.11083785});
#dynsVac2 = VacuumDynamics(film1,25,activities={DMF:1,DMSO:1,GBL:1},number_of_time_steps=1000);

sq = SimulationQueue([["d_mu",1.2]], [dynsVac1,dynsVac2])
sq.startQueue();
sq.plotAllThicknessCompositionAndSupersaturation(axs,[abs_path+'/../dyn_interf_measurements/Double_Cation.npy'],['Interf. Meas.']);

plt.savefig('2Cat_VACUUM.png')


fig,axs = plt.subplots(nrows=2,figsize=st.figsize()) 

corr = AverageRoundJet(DRY_AIR.gas,32,3,1,100)
sp = SimulationParameters(corr,22);
film1 = Film(sp,[DMF,MAPI],4.7, [1], 1.18)
dynsDynGas1 = GasQuenchDynamics(film1,280);
dynsDynGas1.computeSolution();
dynsDynGas1.plotFilmThicknessCompositionAndSupersaturation(axs);


plt.savefig('MAPI_Gas.png')


fig,axs = plt.subplots(nrows=2,figsize=st.figsize()) 
sp = SimulationParameters(corr,22);
film1 = Film(sp,[GBL,DMSO,DMF,TWOCAT],7.3, [1/3*1.3,0.2*2/3*1.3,0.8*2/3*0.7], 0.86)

dynsDynGas1 = GasQuenchDynamics(film1,5000);
dynsDynGas2 = GasQuenchDynamics(film1,5000,activities={DMF:1,DMSO:0.16479906,GBL:0.11083785});
#dynsDynGas2 = GasQuenchDynamics(film1,3200,activities={DMF:1,DMSO:1,GBL:1});

sq = SimulationQueue([["d_mu",1.2]], [dynsDynGas1,dynsDynGas2])
sq.startQueue();
sq.plotAllThicknessCompositionAndSupersaturation(axs,xytext=(10,10));

plt.savefig('2Cat_Gas.png')


fig,axs = plt.subplots(nrows=2,figsize=st.figsize())
corr = AngleImpingingSlotJet(DRY_AIR.gas,0,3,0.1,100);
sp = SimulationParameters(corr,22);

film1 = Film(sp,[DMF,MAPI],4.7, [1], 1.18)

dynsDynGas = DynamicGasQuenchDynamics(film1,35,v_mps=.00595,L_m=.15);
dynsDynGas.computeSolution();
dynsDynGas.plotFilmThicknessCompositionAndSupersaturation(axs);

plt.savefig('MAPI_dynamicsGas.png')


fig,axs = plt.subplots(nrows=2,figsize=st.figsize())
corr = AngleImpingingSlotJet(DRY_AIR.gas,0,3,.1,100);
sp = SimulationParameters(corr,22);
film1 = Film(sp,[GBL,DMSO,DMF,TWOCAT],7.3, [1/3*1.3,0.2*2/3*1.3,0.8*2/3*0.7], 0.86)
dynsDynGas1 = DynamicGasQuenchDynamics(film1,550,v_mps=.00035,L_m=.15);

dynsDynGas2 = DynamicGasQuenchDynamics(film1,550,v_mps=.00035,L_m=.15,activities={DMF:1,DMSO:0.16479906,GBL:0.11083785});


sq = SimulationQueue([["d_mu",1.2]], [dynsDynGas1,dynsDynGas2])
sq.startQueue();
sq.plotAllThicknessCompositionAndSupersaturation(axs);
plt.savefig('2Cat_dynamicsGas.png')

fig,axs = plt.subplots(nrows=2,figsize=st.figsize())
corr = LiquidRoundJet(CB.liq,32,3,1,0.1)

sp = SimulationParameters(corr,22);

film1 = Film(sp,[DMF,MAPI],2, [1], 1.18)

dyns = AntiSolventDynamics(film1,0.03,antiSolvent=CB);
dyns.computeSolution();
dyns.plotFilmThicknessCompositionAndSupersaturation(axs);
plt.savefig('MAPI_antiSolvent.png')


fig,axs = plt.subplots(nrows=2,figsize=st.figsize())
film1 = Film(sp,[GBL,DMSO,DMF,TWOCAT],2, [1/3*1.3,0.2*2/3*1.3,0.8*2/3*0.7], 0.86)

dyns = AntiSolventDynamics(film1,0.03,antiSolvent=CB);
dyns.computeSolution();
dyns.plotFilmThicknessCompositionAndSupersaturation(axs);
plt.savefig('2cat_antiSolvent.png')

plt.show()