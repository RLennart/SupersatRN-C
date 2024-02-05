# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 14:14:27 2023

@author: Simon Ternes
"""
from supersatrnc import *;
import pylab as plt;
import style;
import numpy as np;
jet_chin = AngleImpingingSlotJetRefitted(DRY_AIR.gas,100,1,1,100,90)
jet_martin = AverageSlotJet(DRY_AIR.gas,100,1,1,100);
jet_round = AverageRoundJet(DRY_AIR.gas,100,1,1,100);
jet_laminar = LaminarAverage(DRY_AIR.gas,100,1,100,100);
jet_liquid = LiquidRoundJet(CB.liq,100,1,1,0.1)
T = 293;
p = 101325
xs = np.linspace(0.005,0.045,1000);



betas_chin = np.zeros(len(xs));
betas_martin = np.zeros(len(xs));
betas_laminar = np.zeros(len(xs));
betas_liquid = np.zeros(len(xs));
betas_round = np.zeros(len(xs))

for i,x in enumerate(xs):
    jet_chin.x_m = x;
    jet_martin.x_m = x;
    jet_laminar.x_m = x;
    jet_laminar.D_m = x;
    jet_liquid.x_m = x;
    jet_round.x_m = x;
    betas_chin[i] = jet_chin.beta(T, p, DMF.gas)
    betas_martin[i] = jet_martin.beta(T, p, DMF.gas)
    betas_round[i] = jet_round.beta(T, p, DMF.gas)
    betas_laminar[i] = jet_laminar.beta(T, p, DMF.gas)
    betas_liquid[i] = jet_liquid.beta(T, p, DMF.liq)
    

ind_plot = 3;



sh_chin = np.cumsum(xs*betas_chin/DRY_AIR.gas.diff_m2ps(T, p,DMF.gas))*(xs[1]-xs[0])/xs;
sh_martin = betas_martin*xs/DRY_AIR.gas.diff_m2ps(T, p,DMF.gas);
sh_laminar = betas_laminar*xs/DRY_AIR.gas.diff_m2ps(T, p,DMF.gas);
sh_round = betas_round*xs/DRY_AIR.gas.diff_m2ps(T, p,DMF.gas);
sh_liquid = betas_liquid*xs/CB.liq.diff_m2ps(T, p,DMF.liq);

plt.figure();

plt.plot(xs[ind_plot:-ind_plot],sh_chin[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_martin[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_round[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_laminar[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_liquid[ind_plot:-ind_plot]);
plt.xlabel("Substrate Width x[m]")
plt.ylabel("Avg. Sherwood number ")
plt.legend(["Angle Impinging Slot Jet (Gas dyn.)","Average Slot Jet (Gas stat.)","Average Round Jet (Gas stat.)","Average Laminar (Vacuum)","Average Liquid Round Jet (Antisolvent)"])
plt.yscale('log')

sh_chin_loc = np.gradient(xs*sh_chin,xs)

sh_martin_loc = np.gradient(xs*sh_martin,xs)
sh_round_loc = np.gradient(xs*sh_round,xs)
sh_laminar_loc = np.gradient(xs*sh_laminar,xs)
sh_liquid_loc = np.gradient(xs*sh_liquid,xs)



plt.figure();

plt.plot(xs[ind_plot:-ind_plot],sh_chin_loc[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_martin_loc[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_round_loc[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_laminar_loc[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_liquid_loc[ind_plot:-ind_plot]);
plt.xlabel("Position x[m]")
plt.ylabel("local Sherwood number ")
plt.legend(["Angle Impinging Slot Jet (Gas dyn.)","Average Slot Jet (Gas stat.)","Average Round Jet (Gas stat.)","Average Laminar (Vacuum)","Average Liquid Round Jet (Antisolvent)"])
plt.yscale('log')

plt.figure();

plt.plot(xs[ind_plot:-ind_plot],sh_chin_loc[ind_plot:-ind_plot]*DRY_AIR.gas.diff_m2ps(T, p,DMF.gas)/xs[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_martin_loc[ind_plot:-ind_plot]*DRY_AIR.gas.diff_m2ps(T, p,DMF.gas)/xs[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_round_loc[ind_plot:-ind_plot]*DRY_AIR.gas.diff_m2ps(T, p,DMF.gas)/xs[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_laminar_loc[ind_plot:-ind_plot]*DRY_AIR.gas.diff_m2ps(T, p,DMF.gas)/xs[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_liquid_loc[ind_plot:-ind_plot]*CB.liq.diff_m2ps(T, p,DMF.liq)/xs[ind_plot:-ind_plot]);
plt.xlabel("Position x[m]")
plt.ylabel("Local Mass Transfer coeff $\\beta$ [m/s] ")
plt.legend(["Angle Impinging Slot Jet (Gas dyn.)","Average Slot Jet (Gas stat.)","Average Round Jet (Gas stat.)","Average Laminar (Vacuum)","Average Liquid Round Jet (Antisolvent)"])
plt.yscale('log')



plt.figure();
plt.plot(xs[ind_plot:-ind_plot],sh_chin[ind_plot:-ind_plot]*DRY_AIR.gas.diff_m2ps(T, p,DMF.gas)/xs[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_martin[ind_plot:-ind_plot]*DRY_AIR.gas.diff_m2ps(T, p,DMF.gas)/xs[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_round[ind_plot:-ind_plot]*DRY_AIR.gas.diff_m2ps(T, p,DMF.gas)/xs[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_laminar[ind_plot:-ind_plot]*DRY_AIR.gas.diff_m2ps(T, p,DMF.gas)/xs[ind_plot:-ind_plot]);
plt.plot(xs[ind_plot:-ind_plot],sh_liquid[ind_plot:-ind_plot]*CB.liq.diff_m2ps(T, p,DMF.liq)/xs[ind_plot:-ind_plot]);
plt.xlabel("Substrate Width x[m]")
plt.ylabel("Avg. Mass Transfer coeff $\\beta$ [m/s] ")
plt.legend(["Angle Impinging Slot Jet (Gas dyn.)","Average Slot Jet (Gas stat.)","Average Round Jet (Gas stat.)","Average Laminar (Vacuum)","Average Liquid Round Jet (Antisolvent)"])
plt.yscale('log')





