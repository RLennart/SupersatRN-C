# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 11:37:24 2023

@author: Simon Ternes

    Supersaturation Rate Numerical Calculator (main module) - Calculates the dynamic thickness and supersaturation of perovskite solution films 
    
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
import numpy as np;
from scipy.integrate import solve_ivp;
from scipy.interpolate import interp1d;
from supersatrnc.materials import *;
from scipy.optimize import curve_fit;
from supersatrnc.correlations import Correlation;
import numpy.typing as npt;
from matplotlib.axes import Axes;
from matplotlib.collections import PolyCollection;
import copy;
import pylab as plt;



class SimulationParameters():

    def __init__(self, corr: Correlation, T_C: float, p_pa: float = 101325, t_s: float = 0) -> None:
        self.T = T_C + 273.15;
        self.p0 = p_pa;
        self.t0 = t_s;
        self.corr = corr;


class SimulationVariables():

    def __init__(self, sp):
        self.p = sp.p0;  # pressure outside
        self.t = sp.t0;  # start time simulation


class Film():
    
    def __init__(self,sp:SimulationParameters,compounds:list[Compound],d_init_mum:float,volume_ratios_solvent:list[float],d_final_mum:float,volume_ratios_solute:list[float]=None) -> None:
        
        self.sp = sp;

        'Key lists defining the order'
        self.slvs = list();
        self.slts = list();
        self.only_slvs = list();
        self.only_aslvs= list();
        self.compounds = list();
        
        if(volume_ratios_solute==None):
            volume_ratios_solute = [1/((len(compounds)-len(volume_ratios_solvent)))]*(len(compounds)-len(volume_ratios_solvent));
        
        
        self.Ni0s_mol = dict();
        self.Np0s_mol = dict();
        
        #self.__di0s_m = dict();
        #self.__dp0s_m = dict();
        
        self.solution = dict();
        
        self.dp_m = d_final_mum/1000/1000;
        self.d0_m = d_init_mum/1000/1000-self.dp_m;
        
        for i,c in enumerate(compounds):
            if(not c.isSolute()):
                r = volume_ratios_solvent.pop(0);
            else:
                r = volume_ratios_solute.pop(0);
                
            self.addCompound(c,r);

        # list of compounds (e.g. DMF and MAPI)
        self.compounds = self.slts+self.slvs;    
        
        self.Np_mol = np.sum(np.array(list(self.Np0s_mol.values())));
      
    # this function should be defined in film setup class, logic from above should be in here
    def addCompound(self,c:Compound,init_ratio:float) -> None:
        if(not c.isSolute()):           
            self.slvs.append(c);

            if(c.liq.isAntiSolvent()):  
                self.only_aslvs.append(c);
            else:
                self.only_slvs.append(c);
                
            self.Ni0s_mol[c] = self.d0_m*init_ratio*c.liq.molDens_molpm3(self.sp.T,self.sp.p0);
        
        else:
            self.slts.append(c);
            self.Np0s_mol[c] = self.dp_m*init_ratio*c.sld.molDens_molpm3(self.sp.T,self.sp.p0);

    def Ni0s(self) -> npt.NDArray[np.float64]:
        return np.array([self.Ni0s_mol[slv] for slv in self.slvs]);
    
    def di0_m(self,slv:Solvent) -> npt.NDArray[np.float64]:
        return self.Ni0s_mol[slv]/slv.liq.molDens_molpm3(self.sp.T,self.sp.p0);
        
    def di0s(self) -> npt.NDArray[np.float64]:
        return np.array([self.di0_m(slv) for slv in self.slvs]);
    
    def Npj0s(self) -> npt.NDArray[np.float64]:
        return np.array([self.Np0s_mol[slt] for slt in self.slts]);
    
    def dpj0_m(self,slt:Crystal) -> npt.NDArray[np.float64]:
        return self.Np0s_mol[slt]/slt.sld.molDens_molpm3(self.sp.T,self.sp.p0)
    
    def dpj0s(self) -> npt.NDArray[np.float64]:
        return np.array([self.dpj0_m(slt) for slt in self.slts]);
    
    def averageMolDensP_molpm3(self) -> npt.NDArray[np.float64]:
        return np.sum(self.dp0s()/self.dp_m*[slt.sld.molDens_molpm3(self.sp.T,self.sp.p0) for slt in self.slts])
        
    def getFilmWithAddedAntiSolvent(self,antiSolvent:Solvent) -> None:

        ret = copy.deepcopy(self);
        ret.addCompound(antiSolvent,0);
        ret.compounds = ret.slts+ret.slvs;   
    
        return ret;

    def compoundNames(self) -> list[str]:
        return [k.ident_string for k in self.compounds];
    

class CalculationResult():
    
    def __init__(self,film:Film,times:npt.NDArray[np.float64],sol_data:npt.NDArray[np.float64]) -> None:
        
        self.film = film;
        self.sp = film.sp;
        self.solution = dict();
        self.times = times;

        for i,s in enumerate(sol_data):
            self.solution[film.slvs[i]] = s;

        self.findCritSuperSatAndRate();

    def thicknesses_m(self,liquids:list[Solvent]=[],offset:float=0) -> npt.NDArray[np.float64]:
        if(liquids == []):
            liquids = self.film.slvs;
            
        ths = 0;    
        for i,s in enumerate(liquids):
            ths = ths + self.solution[s]/s.liq.molDens_molpm3(self.sp.T,self.sp.p0);            
        ths = ths + offset;
        return  ths;

        
    def allThicknessesAndColors(self,max_ind:int) -> npt.NDArray[np.float64]:
        
         ths_sld = [np.array([self.film.dpj0_m(slt)*1000*1000]*(max_ind),dtype=float) for slt in self.film.slts];
         
         ths_slv = [np.array(self.solution[k][0:max_ind]/k.liq.molDens_molpm3(self.sp.T,self.sp.p0)*1000*1000,dtype=float) for k in self.film.slvs];
         
         colors = [k.color for k in self.film.compounds]

         return ths_sld+ths_slv,colors;
    
    #
    # def filmWithShiftedInitState(self,index:int=-1) -> npt.NDArray[np.float64]:
    #
    #     d_init = self.thicknesses_m()[index]*1000*1000;
    #
    #     vrs = np.array([self.solution[s][index]/s.liq.molDens_molpm3(self.sp.T,self.sp.p0) for s in self.film.slvs])/(self.thicknesses_m()[index]-self.film.dp_m);
    #
    #     vrps = np.array(self.film.dp0s_m())/self.dp_m;
    #
    #     return Film(SimulationParameters(self.sp.corr,self.sp.T,self.sp.p0,0),self.compounds, d_init, list(vrs), self.dp_m*1000*1000, list(vrps));
    #
    #
    # def totaldisTimesWeights_m(self,liquids:list[Solvent]=[],weights:npt.NDArray[np.float64]=-1) -> npt.NDArray[np.float64]:
    #
    #
    #     if(liquids==[]):
    #         liquids=self.film.slvs;
    #
    #     if(weights==-1):
    #         weights = np.ones(len(liquids));
    #
    #     res = 0;
    #
    #     for i,slv in enumerate(liquids):
    #         res = res + weights[i]*self.solution[slv]/slv.liq.molDens_molpm3(self.sp.T,self.sp.p0);
    #
    #     return res;

    
    def c_lin(self,crit:bool=False) -> npt.NDArray[np.float64]:
        
        lenSol = len(self.solution[self.film.slvs[0]]);
        
        ret = np.zeros(lenSol)
        
        for slt in self.film.slts:
            add = np.zeros(lenSol)

            
            for slv in self.film.only_slvs:
                ### Summer über richtigen INDEX!!
                ri = self.solution[slv]/slv.liq.molDens_molpm3(self.sp.T,self.sp.p0)/self.thicknesses_m(self.film.only_slvs);          
                
                ri[np.where(np.isnan(ri))] = 1;
                
                if(not crit):
                    add = add + slt.sld.ceq[slv]*ri;
                else:
                    add = add + slt.sld.ccrit[slv]*ri;
                
                for aslv in self.film.only_aslvs:
                    
                    if(aslv in slt.sld.anti_qs):
                        
                        add = add - ri*slt.sld.anti_qs[aslv][slv]*self.solution[aslv]/aslv.liq.molDens_molpm3(self.sp.T,self.sp.p0)/self.thicknesses_m(self.film.only_slvs)
                        
            ret = ret + add*self.film.dpj0_m(slt)/self.film.dp_m; 
        
        #Concentrations are always greater zero - todo what?? why is this code required at all?
        for i,r in enumerate(ret):
            if(r<0):
                ret[i] = 0;
                
        return ret;


    def findCritSuperSatAndRate(self) -> npt.NDArray[np.float64]:
         
         self.crit_index = 0;

        # todo why errstate?
         with np.errstate(divide='ignore', invalid='ignore'):
             # print(self.film.Np_mol,self.thicknesses_m(offset=self.film.dp_m),self.c_lin())
             # return self.c_lin();
             ret = np.log(self.film.Np_mol / self.thicknesses_m(offset=self.film.dp_m) / self.c_lin());

         self.supersats = ret;
         ceqs = self.c_lin(False);
         
         crits = self.c_lin(True);

         threshold = np.log(crits[0]/ceqs[0])
         
         #print(self.crit_index<(len(supersats)-1),np.isnan(supersats[self.crit_index]),supersats[self.crit_index] < threshold,supersats[self.crit_index],threshold)
         
         # todo this is lol
         while (self.crit_index<(len(self.supersats)-1) and (np.isnan(self.supersats[self.crit_index]) or self.supersats[self.crit_index] < threshold)):
             self.crit_index = self.crit_index+1;
         
         self.superSatRate = np.gradient(self.supersats,(self.times[1]-self.times[0]));
         self.crit_superSat = self.supersats[self.crit_index];
         self.crit_superSatRate = self.superSatRate[self.crit_index];    

    
    def getFilmWithShiftedInitState(self,index:int=-1) -> Film:
        
        d_init = self.thicknesses_m(offset=self.film.dp_m)[index]*1000*1000;
        
        vrs = [self.thicknesses_m([slv])[index]/self.thicknesses_m()[index] for slv in self.film.slvs];
        
        #print("-------- VRS",vrs)
        
        vrps = self.film.dpj0s()/self.film.dp_m;

        return Film(self.film.sp,self.film.compounds, d_init, list(vrs), self.film.dp_m*1000*1000, list(vrps));
    

    # def djs(self) -> npt.NDArray[np.float64]:
    #     return self.film.d0_m-self.totaldisTimesWeights_m();
    #
    # def Vratios(self) -> npt.NDArray[np.float64]:
    #     return self.djs()/self.totaldisTimesWeights_mol();


class DynamicsSimulation:
    
    def __init__(self,film:Film,end_time_s:float,number_of_time_steps:int=10000) -> None:
        self.film = film;
        self.end_time_s = end_time_s;
        self.number_of_time_steps = number_of_time_steps;

        if(self.film.sp.t0 > self.end_time_s):
            raise AttributeError("t0 cannot be more than maximum simulation time");
        self.times = np.linspace(self.film.sp.t0,self.end_time_s,self.number_of_time_steps);

        self.result = None;
        self.sv = SimulationVariables(self.film.sp);


    def computeSolution(self) -> None:
        
        self.sv = SimulationVariables(self.film.sp);

        # this is the expensive line
        sol = solve_ivp(self.Niprimes_molps, [self.times[0], self.times[-1]],self.film.Ni0s(),dense_output=True,max_step=self.times[1]-self.times[0]);

        # this is weird. here the y-solution is extracted and then fed into Calculation Result
        soly = sol.sol(self.times); # this is likely more-dimensional if multiple solvent components inside
        
        self.result = CalculationResult(self.film, self.times, soly);
        

    def Niprimes_molps(self,time:float,N:list[float]) -> list[float]:
        ret = list();

        for i,s in enumerate(self.film.slvs):
            
            #If Ni is sufficintly small, the drying rate will be zero
            if(N[i]>1e-10*self.film.Np_mol):
                ret.append(-self.prop(s,i,N)*N[i]);
            else:
                ret.append(0);
        
        self.timeStep(time); # set correct time for time dependent properties;

        return ret;

    def timeStep(self,time:float) -> None:
            
        self.sv.t = time;
    
    
    def incrementIndex(self,max_ind:int) -> int:
        if(max_ind == -1):
            return len(self.times);
    
        return max_ind+1;

    def plotTotalFilmThickness(self,ax:Axes,interferometricDataFiles:list[str]=[],time_div:float=0,max_ind:int=-1) -> None:
        
        max_ind = self.incrementIndex(max_ind);
        
        for f in interferometricDataFiles:
            data = np.load(f)
            
            #append start and end points
            xdata = np.insert(np.append(data[0],self.times[0:max_ind][-1]),0,self.times[0:max_ind][0])
            ydata = np.insert(np.append(data[1],self.result.film.dp_m*1000*1000),0,(self.result.film.d0_m+self.result.film.dp_m)*1000*1000)
            
            ax.plot(xdata+time_div,ydata,'o',color='black')
                
        p2, = ax.plot(self.times[0:max_ind]+time_div,self.result.thicknesses_m()*1000*1000);
        ax.set_ylabel("Thickness [$\mu$m]")
        ax.set_xlabel("Time [s]")
    
    def plotThicknessComposition(self,ax:Axes,max_ind:int,time_div:float) -> list[PolyCollection]:
        
        ths, colors = self.result.allThicknessesAndColors(max_ind);
        
        #print("times",len(self.times[0:max_ind]+time_div));
        #print("ths",[len(ths[i]) for i in range(len(ths))]);
        
        polys = ax.stackplot(np.array(self.times[0:max_ind]+time_div,dtype=float),ths,colors=colors);

        #print(ths,colors)
        
        return polys;

    
    def plotFilmThicknessCompositionAndSupersaturation(self,axes:tuple[Axes],interferometricDataFiles:list[str]=[],dataLabels:list[str]=[],xytext:tuple[float]=(10,-10),max_ind:int=-1,time_div:float=0) -> None:
        
        max_ind = self.incrementIndex(max_ind);
        
        
        firstPlot = True;
        if(len(axes)==2):
            top,bottom = axes;
        else:
            if(len(axes)==3):
                firstPlot = False;
                top,bottom,par2 = axes;
            else:
                raise Exception("Pass two Axes subplot or three Axes")
            
        if(firstPlot):
            top.set_ylabel("Supersat. Rate [s$^{-1}$]")
            bottom.set_xlabel("Time t [s]")
            bottom.set_ylabel("Composition of film th. [$\mu$m]")

        p1, =  top.plot(self.times[2:max_ind]+time_div,self.result.superSatRate[2:max_ind],color='firebrick'); #first two points of drivative might be flawed.


        if(self.result.crit_index < max_ind-1):
            top.axvline(x=self.times[self.result.crit_index]+time_div,linestyle='--')
            top.plot(self.times[self.result.crit_index]+time_div,self.result.crit_superSatRate,'o',color='firebrick', markersize =15);
        
            top.annotate("$\\dot{\\sigma}_{crit.}$", (self.times[self.result.crit_index]+time_div, self.result.crit_superSatRate),textcoords="offset points",
                         xytext=xytext,color="firebrick")
            bottom.axvline(x=self.times[self.result.crit_index]+time_div,linestyle='--')
            
            print("Determined critical supersat. rate at ",self.result.crit_superSatRate," per second");

        top.set_xlim([0,self.times[-1]])
        bottom.set_xlim([0,self.times[-1]])
        
        
        top.yaxis.label.set_color(p1.get_color())
        top.tick_params(axis='y', colors=p1.get_color())

        polys = self.plotThicknessComposition(bottom,max_ind,time_div)

    
        if(firstPlot):
            par2 = top.twinx()
    
        #print(self.result.supersaturation()[0:max_ind])
    
        p2, = par2.plot(self.times[0:max_ind]+time_div,self.result.supersats[0:max_ind],color='cornflowerblue');
    
        if(firstPlot):
            par2.yaxis.label.set_color(p2.get_color())
            par2.tick_params(axis='y', colors=p2.get_color())
            par2.set_ylabel("Supersaturation")
    
    
        add_lines = list();
    
        for f in interferometricDataFiles:
            data = np.load(f)
            l, = bottom.plot(data[0],data[1],'o',color='black')
            add_lines.append(l);
    
        bottom.legend(polys+add_lines,self.result.film.compoundNames()+dataLabels,loc=1)   
    
        ylabels =  bottom.get_yticklabels()
        ylabels[0].set_alpha(0)
        top.set_ylim([-0.05*self.maximumSupersaturation(),self.maximumSupersaturation()*1.2])

        return top,bottom,par2;

    
    def maximumSupersaturation(self) -> float:
        return np.max(self.result.superSatRate)
    
    
    
class GasQuenchDynamics(DynamicsSimulation):  

    def __init__(self,film:Film,end_time_s:float,number_of_time_steps:int=10000,activities:dict[Compound,float]=dict()) -> None:    
        super().__init__(film,end_time_s,number_of_time_steps)
        self.activities = activities;
    
    def prop(self,m:Compound,i:int,Nis:list[float]) -> float:
        ret = (self.film.sp.corr.beta(self.film.sp.T,self.sv.p,m.gas)
               *m.gas.vaporPressure_Pa(self.film.sp.T)/R_JpmolK/self.film.sp.T/(np.sum(Nis)+self.film.Np_mol));

        #print(self.sv.t,self.film.sp.corr.beta(self.film.sp.T,self.sv.p,m.gas),m.gas.vaporPressure_Pa(self.film.sp.T),self.denom(i, Nis))

        # todo rather set activities to standard 1 if not given
        if(not m in self.activities.keys()):
            return ret
    
        return ret*self.activities[m];
    
    def denom(self,i:int,Nis:list[float]) -> float:
        return (np.sum(Nis)+self.film.Np_mol);


class DynamicGasQuenchDynamics(GasQuenchDynamics):
    
    def __init__(self,film:Film,end_time_s:float,number_of_time_steps:int=10000,activities:dict[Compound,float]=dict(),v_mps:float=1,L_m:float=None) -> None:
        super().__init__(film,end_time_s,number_of_time_steps,activities)
        
        self.L_m=L_m;
        if(L_m==None):
            # it does not jump in here but it would cause errors since not defined
            self.L_m = self.film.corr.D_m*10;
        
        self.v_mps = v_mps;
    
    def timeStep(self,time:float) -> None:
        super().timeStep(time);
        #self.film.sp.corr.x_m = self.v_mps*(time+self.film.sp.t0)-self.L_m;
        self.film.sp.corr.x_m = self.v_mps*(time)-self.L_m;
        
        
class VacuumDynamics(GasQuenchDynamics):
    
    def __init__(self,film:Film,end_time_s:float,number_of_time_steps:int=10000,activities:dict[Compound,float]=dict(),p_params:list[float]=[0.55143471, 7.74808886 , 0.02585276],u0_params:list[float]=[0.09453549, 5.2 ]):
        super().__init__(film,end_time_s,number_of_time_steps,activities)
        self.pp = p_params;
        self.u0p = u0_params;
        self.setPressureAndVel();

    def timeStep(self,time:float) -> None:
        super().timeStep(time);
        
        self.setPressureAndVel();

    def setPressureAndVel(self) -> None:
        self.sv.p = self.film.sp.p0*np.exp(-self.pp[0]*(self.sv.t))+self.pp[1]-self.pp[2]*(self.sv.t);
        self.film.sp.corr.u0_mps = self.u0p[1]*np.exp(-self.u0p[0]*(self.sv.t));

    
class SimulationQueue:
    def __init__(self,onsets:list[list[str,float]],simulations:list[DynamicsSimulation]) -> None:
        self.onsets = onsets;
        self.onsets.append(["t_s",simulations[-1].times[-1]])
        self.simulations = simulations;
        
        
    def startQueue(self):
        
        self.bnd_inds= list();
        self.div_times = list();
        
        bnd_time = 0;
        
        for i,s in enumerate(self.simulations):
            self.div_times.append(bnd_time);            

            o = self.onsets[i]

            s.computeSolution();
            
            if(o[0] == 't_s'):
                bnd_ind = np.argmin(np.abs(s.times-o[1]));
            else: 
                if(o[0] == 'd_mu'):
                    bnd_ind = np.argmin(np.abs(s.result.thicknesses_m(offset=s.film.dp_m)*1000*1000-o[1]));
                else:
                    raise Exception("Boundary ill specified")
            
            bnd_time = s.times[bnd_ind];
            
            next_film = s.result.getFilmWithShiftedInitState(bnd_ind);
            
            self.bnd_inds.append(bnd_ind);
    
            if(i<len(self.simulations)-1):
                #last_film.sp = self.simulations[i+1].film.sp;
                self.simulations[i+1].film = next_film;
                self.simulations[i+1].adjustStartTime(bnd_time);
                
            
            
    def plotAllThicknessCompositionAndSupersaturation(self,axes,interferometricDataFiles=[],dataLabels=[],xytext=(10,10)):        
        for i,s in enumerate(self.simulations):
            axes = s.plotFilmThicknessCompositionAndSupersaturation(axes,interferometricDataFiles,dataLabels,xytext=xytext,max_ind=self.bnd_inds[i],time_div=0);
        
        max_supersat = np.max([s.maximumSupersaturation() for s in self.simulations])
        axes[0].set_ylim([-0.05*max_supersat,max_supersat*1.2])
        
    def thicknesses_m(self,liquids=[], offset=0):

        thicknesses = self.simulations[0].result.thicknesses_m(liquids,offset)[:self.bnd_inds[0]];
        #print(self.bnd_inds,len(thicknesses))

        for i in range(0,len(self.simulations)-2):
            s= self.simulations[i+1];
            new_ths = s.result.thicknesses_m(liquids,offset)[:self.bnd_inds[i+1]];
            thicknesses = np.concatenate((thicknesses,new_ths));
        
        last_ths = self.simulations[-1].result.thicknesses_m(liquids,offset);
        thicknesses = np.concatenate((thicknesses,last_ths));

        return thicknesses;

    def times_s(self):
        times = self.simulations[0].times[:self.bnd_inds[0]];


        for i in range(0,len(self.simulations)-2):
            s= self.simulations[i+1];
            new_ts = s.times[:self.bnd_inds[i+1]];
            times = np.concatenate((times,new_ts));

        last_times = self.simulations[-1].times;

        times = np.concatenate((times,last_times));

        return times;

class DryAndWetThicknessFitter:
    
    def __init__(self,simulation:DynamicsSimulation,dataToFitFile:str,fitRange:list[int,int],init_volumeratios:list[float]):
        self.simulation = simulation;
        data = np.load(dataToFitFile);
        self.xdata = data[0][fitRange[0]:fitRange[1]]
        self.ydata = data[1][fitRange[0]:fitRange[1]]
        self.init_volumeratios = init_volumeratios;
        
    def startFitting(self):
        
        
        rs = curve_fit(self.fitFunction,self.xdata,self.ydata,p0=[self.simulation.film.d0_m+self.simulation.film.dp_m,self.simulation.film.dp_m],bounds=((0,0),(np.inf,np.inf)))
        return rs;
        
    def fitFunction(self,xdata,d0_m,dp_m):
         film = self.simulation.film;
         
         #print(film.slvs+film.slts,[i for i in self.init_volumeratios]);
         
         self.simulation.film = Film(self.simulation.film.sp,film.slvs+film.slts,(d0_m+dp_m)*1000*1000, [i for i in self.init_volumeratios], dp_m*1000*1000);
         self.simulation.computeSolution();
         
         inter_cub = interp1d( self.simulation.times, self.simulation.result.thicknesses_m(offset=0)*1000*1000, kind = 'cubic');
         
         return inter_cub(xdata);
    

# Generalize on Simulation-Queue and then determine + plot;
class ActivityFitter:
    
    def __init__(self,simulation:DynamicsSimulation,dataToFitFile:str,fixedParams:list[str]=[]):
        self.simulation = simulation;
        
        self.sq = None;
        if(not isinstance(self.simulation,DynamicsSimulation )): 
            self.sq = self.simulation;
            self.simulation = self.sq.simulations[-1]; #only fit the activities of last simulation
            self.sq.startQueue();
            self.times = self.sq.times_s();
            
        data = np.load(dataToFitFile);
            
        self.xdata = np.insert(data[0],0,self.times[-1])
        self.ydata = np.insert(data[1],0,self.simulation.film.dp_m*1000*1000)

        self.fixedParams = fixedParams;
    
        


    def startFitting(self):
        
        inter_cub = interp1d(self.xdata, self.ydata, kind = 'linear');
        
        fitted_activities = dict();
        
        keys = list(self.simulation.activities.keys());
        
        for i,a in enumerate(list(self.simulation.activities.values())):
            if(not keys[i] in self.fixedParams):
                fitted_activities[keys[i]] = np.abs(a);

        #print(self.times);

        #print(inter_cub(self.times));

        rs = curve_fit(self.fitFunction,self.times,inter_cub(self.times),p0=list(fitted_activities.values()),bounds=[[0]*len(fitted_activities),[1.0012]*len(fitted_activities)])

        # plt.figure()

        # plt.plot(self.xdata,inter_cub(self.xdata));
        # plt.plot(self.times,self.fitFunction(self.times,*rs[0]))

        # plt.show()

        return rs;
        
    def fitFunction(self,xdata,*activities):
        keys = list(self.simulation.activities.keys());
        
        new_acts = dict();
        
        ind_act = 0;
        
        for i,a in enumerate(list(self.simulation.activities.values())):
            if(not keys[i] in self.fixedParams):
                new_acts[keys[i]] = np.abs(activities[ind_act]); #only accept positive values
                ind_act = ind_act+1;
            else:
                new_acts[keys[i]] = self.simulation.activities[keys[i]];
            
        self.simulation.activities = new_acts;
        
        #print(new_acts.keys(),new_acts.values());

        if(self.sq == None):        
            self.simulation.computeSolution();
            return  self.simulation.result.thicknesses_m(offset=self.simulation.film.dp_m)*1000*1000;
        else: 
            self.sq.startQueue();
            return self.sq.thicknesses_m(offset=self.simulation.film.dp_m)*1000*1000;

class ParameterSweeper:

    def __init__(self,simulation:DynamicsSimulation,corr_argument:str,param_range:npt.NDArray[np.float64],xlabel:str=''):

        self.simulation = simulation;
        self.param_range = param_range;
        self.corr_arg = corr_argument;
        self.xlabel = xlabel;

    def calculateSweep(self):

        self.crit_supersats = np.zeros(np.shape(self.param_range));

        for i,p in enumerate(self.param_range):
            self.setAttribute(p)
            #oldfilm = copy.deepcopy(self.simulation.film)
            self.simulation.computeSolution();

            self.crit_supersats[i] = self.simulation.result.crit_superSatRate

            #self.simulation.film = oldfilm;

    def plotResults(self,ax):
        
        ax.plot(self.param_range,self.crit_supersats,'o')

        if(self.xlabel==''):
            ax.set_xlabel(self.corr_arg);
        else:
            ax.set_xlabel(self.xlabel);

        ax.set_ylabel("Critical supersaturation rate [s$^{-1}$]")



    def setAttribute(self,p):
        setattr(self.simulation.film.sp.corr,self.corr_arg,p);