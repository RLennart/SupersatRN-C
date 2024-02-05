from supersatrnc import *;
from pylab import plt;
import style as st;
import os;
import numpy as np;

abs_path = os.path.dirname(os.path.abspath(__file__));

plt.figure(figsize=(4*2.2,3*2.2))

ax = plt.gca();


corr = LaminarAverage(DRY_AIR.gas,32,5,32,1)
sp = SimulationParameters(corr,22);

film1 = Film(sp,[DMF,MAPI],5.3, [1], 1.18)

dynsVac = VacuumDynamics(film1,125);

pSweeper = ParameterSweeper(dynsVac,'x_m',np.linspace(0.01,1,100),'Substrate Size [m]');

pSweeper.calculateSweep();



pSweeper.plotResults(ax);



#plt.savefig(abs_path+'Vacuum_Size_Var.png');


corr = AverageRoundJet(DRY_AIR.gas,32,3,1,100)
sp = SimulationParameters(corr,22);
film1 = Film(sp,[DMF,MAPI],4.7, [1], 1.18)
dynsDynGas1 = GasQuenchDynamics(film1,1000);

pSweeper = ParameterSweeper(dynsDynGas1,'x_m',np.linspace(0.01,1,100),'Substrate Size [m]');

pSweeper.calculateSweep();


pSweeper.plotResults(ax);


corr = LiquidRoundJet(CB.liq,32,3,1,0.1)

sp = SimulationParameters(corr,22);

film1 = Film(sp,[DMF,MAPI],2, [1], 1.18)

dyns = AntiSolventDynamics(film1,0.05,antiSolvent=CB);


pSweeper = ParameterSweeper(dyns,'x_m',np.linspace(0.01,1,100),'Substrate Size [m]');

pSweeper.calculateSweep();



pSweeper.plotResults(ax);

plt.yscale('log')

plt.legend(['Vacuum','Static Gas','Antisolvent'])

plt.savefig(abs_path+'Sizes_Var.png');


class ViscositySweeper(ParameterSweeper):

    def setAttribute(self, p):
        self.simulation.antiSolvent._dynViscosity_Pas = p;
        self.simulation.film.sp.corr.pm._dynViscosity_Pas = p;


corr = LiquidRoundJet(VARCB.liq,32,3,1,0.1)

sp = SimulationParameters(corr,22);

film1 = Film(sp,[DMF,MAPI],2, [1], 1.18)

dyns = AntiSolventDynamics(film1,0.05,antiSolvent=VARCB);


pSweeper = ViscositySweeper(dyns,'',np.linspace(1,5,100)/3*CB.liq.dynViscosity_Pas(sp.T,sp.p),'Viscosity [Pas]');

pSweeper.calculateSweep();

plt.figure(figsize=(4*2.2,3*2.2))

ax = plt.gca();

pSweeper.plotResults(ax);

print(abs_path+'Viscosity_Var.png')
plt.savefig(abs_path+'Viscosity_Var.png');

plt.show();