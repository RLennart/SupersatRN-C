# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 10:48:26 2023

@author: cassi
"""

#from supersatrnc import DRY_AIR, DMF, DMSO, GBL, TWOCAT, MAPI
from all import AngleImpingingSlotJet, SimulationParameters, Film, DynamicGasQuenchDynamics, DRY_AIR, DMF, DMSO, GBL, TWOCAT, MAPI
from pylab import plt;
import style as st;
import os;

# fragen
# wieso macht ein so kleiner anteil dmso so einen großen unterschied in der optischen wahrnehmung der samples während/nach dem quenching
#

abs_path = os.path.dirname(os.path.abspath(__file__));

corr = AngleImpingingSlotJet(DRY_AIR.gas,0,3,0.1,100);
sp = SimulationParameters(corr,22);


fig,axs = plt.subplots(nrows=2,figsize=st.figsize(1.1))
film1 = Film(sp,[DMF,MAPI],4.7, [1], 1.18) # d_pm = d_final_mum
dynsDynGas = DynamicGasQuenchDynamics(film1,35,v_mps=.00595,L_m=.15); # vmps coating speed
dynsDynGas.computeSolution();
dynsDynGas.plotFilmThicknessCompositionAndSupersaturation(axs);


fig,axs = plt.subplots(nrows=2,figsize=st.figsize(1.1))
corr = AngleImpingingSlotJet(DRY_AIR.gas,0,3,.1,100);
sp = SimulationParameters(corr,22);
film1 = Film(sp,[GBL,DMSO,DMF,TWOCAT],7.3, [1/3*1.3,0.2*2/3*1.3,0.8*2/3*0.7], 0.86)
film2 = Film(sp,[GBL,DMSO,DMF,TWOCAT],7.3, [1/3*1.3,0.2*2/3*1.3,0.8*2/3*0.7], 0.86)
dynsDynGas1 = DynamicGasQuenchDynamics(film1,550,v_mps=.00035,L_m=.15);
dynsDynGas2 = DynamicGasQuenchDynamics(film2,550,v_mps=.00035,L_m=.15,activities={DMF:1,DMSO:0.16479906,GBL:0.11083785});

dynsDynGas1.computeSolution();
dynsDynGas1.plotFilmThicknessCompositionAndSupersaturation(axs);

dynsDynGas2.computeSolution();
dynsDynGas2.plotFilmThicknessCompositionAndSupersaturation(axs);

plt.show()

fig,axs = plt.subplots(nrows=2,figsize=st.figsize(1.1))
film1 = Film(sp,[DMF, DMSO,MAPI],4.7, [0.95,0.05], 1.18)
dynsDynGas = DynamicGasQuenchDynamics(film1,35,v_mps=.00595,L_m=.15);
dynsDynGas.computeSolution();
dynsDynGas.plotFilmThicknessCompositionAndSupersaturation(axs);


fig,axs = plt.subplots(nrows=2,figsize=st.figsize(1.1))
film1 = Film(sp,[TWOME,MAPI],4.7, [1], 1.18)
dynsDynGas = DynamicGasQuenchDynamics(film1,35,v_mps=.00595,L_m=.15,number_of_time_steps=200);
dynsDynGas.computeSolution();
dynsDynGas.plotFilmThicknessCompositionAndSupersaturation(axs);


fig,axs = plt.subplots(nrows=2,figsize=st.figsize(1.1))
film1 = Film(sp,[TWOME,DMSO,MAPI],4.7, [0.9,0.1], 1.18)
dynsDynGas = DynamicGasQuenchDynamics(film1,35,v_mps=.00595,L_m=.15,number_of_time_steps=200);
dynsDynGas.computeSolution();
dynsDynGas.plotFilmThicknessCompositionAndSupersaturation(axs);

plt.figure()
for d in [3,5,7,9,11,13,15,17,19]:
    film1 = Film(sp, [TWOME, MAPI], d, [1], 1.18)
    dynsDynGas = DynamicGasQuenchDynamics(film1, 35, v_mps=.00595, L_m=.15);
    dynsDynGas.computeSolution();

    plt.plot(dynsDynGas.times, dynsDynGas.result.superSatRate,label='d_init={}'.format(d))

plt.legend()
plt.show()

fig,axs = plt.subplots(nrows=2,figsize=st.figsize())
film1 = Film(sp,[TWOME,MAPI],18, [1], 1.18)
dynsDynGas = DynamicGasQuenchDynamics(film1,35,v_mps=.00595,L_m=.15);
dynsDynGas.computeSolution();
dynsDynGas.plotFilmThicknessCompositionAndSupersaturation(axs);


plt.show()




sq = SimulationQueue([["d_mu",1.2]], [dynsDynGas1,dynsDynGas2])
sq.startQueue();
sq.plotAllThicknessCompositionAndSupersaturation(axs);



fig,axs = plt.subplots(nrows=2,figsize=st.figsize(1.1))
corr = AngleImpingingSlotJet(DRY_AIR.gas,0,3,.1,300);
sp = SimulationParameters(corr,22);
film1 = Film(sp,[GBL,DMSO,DMF,TWOCAT],7.3, [1/3*1.3,0.2*2/3*1.3,0.8*2/3*0.7], 0.86)
dynsDynGas1 = DynamicGasQuenchDynamics(film1,550,v_mps=.00035,L_m=.15);
dynsDynGas1.computeSolution();
dynsDynGas1.plotFilmThicknessCompositionAndSupersaturation(axs);





