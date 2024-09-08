# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 12:31:44 2023

@author: Simon Ternes

    Supersaturation Rate Numerical Calculator (materials module) - Calculates the dynamic thickness and supersaturation of perovskite solution films 
    
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
from supersatrnc.material_properties import R_JpmolK,NA_1pmol,solvents,dryingGases,crystals;
from overrides import override
####To test

class Material:
    
    def __init__(self,ident_string : str) -> None:
        self.ident_string = ident_string;
    
    def __str__(self) -> str:
        return self.ident_string;

    def kinViscosity_m2ps(self,T_K:float,p_Pa:float) -> float:
        return self.dynViscosity_Pas(T_K,p_Pa)/self.massDens_kgpm3(T_K,p_Pa);


    def molDens_molpm3(self,T_K:float,p_Pa:float) -> float:
        raise AttributeError("Density-Function not implemented. Must be overwritten.");

    def massDens_kgpm3(self,T_K:float,p_Pa:float) -> float:
        return self.molDens_molpm3(T_K, p_Pa)*self._molMass_kgpmol;
    
    def specVol_m3pmol(self,T_K:float,p_Pa:float) -> float:
        return 1/self.molDens_molpm3(T_K,p_Pa);  
    
    def dynViscosity_Pas(self,T_K:float,p_Pa:float) -> float:
         raise AttributeError("Viscosity-Function not implemented.");

class Fluid(Material):

    def __init__(self,ident_string : str) -> None:
        super().__init__(ident_string);       

    def schmidt(self,T_K:float,p_Pa:float,inMaterial:Material) -> float:
        return self.kinViscosity_m2ps(T_K,p_Pa)/self.diff_m2ps(T_K,p_Pa,inMaterial);
            
    def reynolds(self, T_K:float,p_Pa:float,vel_mps:float, width_m:float) -> float:
        return vel_mps * np.abs(width_m) / self.kinViscosity_m2ps(T_K,p_Pa);    
    
    
class IdealGas(Fluid): 

    def __init__(self,ident_string : str) -> None:
        super().__init__(ident_string);    

    @override
    def molDens_molpm3(self,T_K:float,p_Pa:float) -> float:
            return p_Pa/R_JpmolK/T_K;
    
    #Fuller-Gas Diffusion
    def diff_m2ps(self,T_K:float,p_Pa:float,gas:'IdealGas')-> float:
        return 0.01*0.01*0.00143*(T_K**1.75)*((self._molMass_kgpmol*1000)**(-1)+(gas._molMass_kgpmol*1000)**(-1))**0.5 / (p_Pa*1e-5 * np.sqrt(2) *(self._diffVol**(1.0/3)+gas._diffVol**(1.0/3))**2)

class InertGas(IdealGas):
    
    def __init__(self,ident_string : str) -> None:
        super().__init__(ident_string);
        gas = dryingGases[ident_string];
        self._molMass_kgpmol = gas['molMass_gpmol']/1000;
        self._diffVol = gas['diffVol'];
        self._sutherCoeffs = gas['dynViscSutherParams'];
    
    #Sutherland Equation for dynamic viscosity -> https://www.cfd-online.com/Wiki/Sutherland%27s_law
    @override
    def dynViscosity_Pas(self,T_K:float,p_Pa:float) -> float: # Sutherland
        C = self._sutherCoeffs[0];
        refT = self._sutherCoeffs[1];
        refDynVisc = self._sutherCoeffs[2];

        return refDynVisc*(refT+C)/(T_K+C)*(T_K/refT)**(1.5);

class SolventGas(IdealGas):
     
     def __init__(self,ident_string : str) -> None:
         super().__init__(ident_string);
         solv = solvents[ident_string];
         self._molMass_kgpmol = solv['molMass_gpmol']/1000;
         self._diffVol = solv['diffVol'];
         self._antoineParams = solv['antoineParams']
     
     #Antoine-Equation https://en.wikipedia.org/wiki/Antoine_equation
     def vaporPressure_Pa(self,T_K:float)->float:
         for tempSet in self._antoineParams:
             A = tempSet[1][0];
             B = tempSet[1][1];
             C = tempSet[1][2];
             if(tempSet[2]):
                 if(tempSet[0][0] < T_K and tempSet[0][1] > T_K):
                     return np.power(10,A-B/(T_K+C))*10**5; # from bar to pascal
             else: #Old Units °C and mmHg
                 if(tempSet[0][0] < T_K-273.15 and tempSet[0][1] > T_K-273.15):
                     return np.power(10,A-B/(T_K-273.15+C))*133.322; # from mmHg to pascal    
                          

class Solvent(Fluid):
    
    def __init__(self,ident_string : str) -> None:
        super().__init__(ident_string);
        solv = solvents[ident_string];
        self._molMass_kgpmol = solv['molMass_gpmol']/1000;
        self._isAntiSolvent = solv['isAntiSolvent'];
        self._molDens_kgpm3 = solv['massDens_gpcm3']*1000
        self._dynViscosity_Pas = solv['dynViscosity_mPas']/1000
        self._surfacetensionParams = solv['surfaceTensionParams']
    
    def isAntiSolvent(self)->bool:
        return self._isAntiSolvent;
    
    @override
    def molDens_molpm3(self,T_K:float,p_Pa:float) -> float:
        return self._molDens_kgpm3/self._molMass_kgpmol
    
    #Tyn-Clausius Liquid-Liquid Diffusion Wärematlas Da28
    def diff_m2ps(self,T_K:float,p_Pa:float,liq:'Solvent') -> float:
        return 0.01*0.01*8.93e-8*(self.specVol_m3pmol(T_K,p_Pa)*100*100*100)**(-1/3)*(liq.specVol_m3pmol(T_K,p_Pa)*100*100*100)**(1/6)*(self.paraChor_m3pmol(T_K,p_Pa)/liq.paraChor_m3pmol(T_K,p_Pa))**(0.6)*T_K*(self.dynViscosity_Pas(T_K,p_Pa)*1000)**(-1)
    
    @override
    def dynViscosity_Pas(self,T_K:float,p_Pa:float) -> float:
        return self._dynViscosity_Pas;
    
    #Linear function of temperature, Reference at 20°C https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
    def surfaceTension_Npm(self,T_K:float,p_Pa:float)-> float:
        params = self._surfacetensionParams;
        return (params[0]+params[1]*(T_K-293))/1000;
    
    #Wärematlas Da 28
    def paraChor_m3pmol(self,T_K:float,p_Pa:float) -> float:
        return self.specVol_m3pmol(T_K,p_Pa)*100*100*100*(self.surfaceTension_Npm(T_K,p_Pa)*1000)*(0.25)
        
class Crystal(Material):
    def __init__(self,ident_string : str) -> None:
        super().__init__(ident_string);
        cry = crystals[ident_string]
        
        self.ceq = cry['eqConcentration_molpm3']
        
        self.ccrit = cry['critConcentration_molpm3'];
        
        self.anti_qs = cry['quenchingEff_molpm3'];
        
        self._latticeConstant_m = cry['latticeConstant_A']*1e-10;
        
        self._molMass_kgpmol = cry['molMass_gpmol']/1000;
    
    @override
    def molDens_molpm3(self,T_K:float,p_Pa:float) -> float:
        return 1/NA_1pmol / (self._latticeConstant_m)**3
    
    def massDens_kgpm3(self,T_K:float,p_Pa:float) -> float:
        return self.molDens_molpm3(T_K,p_Pa)*self._molMass_kgpmol;

    
    
class Compound:   
    
    def __eq__(self, other)->bool:
        if isinstance(other, str): #enabling equality just by comparing ident_string
            return self.ident_string == other;
        if isinstance(other, self.__class__):
            return self.ident_string == other.ident_string;
        else:
            return False

    def __hash__(self)->float:    
        return hash(self.ident_string);
        
    def __ne__(self, other)->bool:
        return not self.__eq__(other);    
    
    def __str__(self) -> str:
        return "Compound: "+self.ident_string
    
    def __init__(self,ident_string:str,color:str)->None:
        self.color = color;
        self.ident_string = ident_string;
        if(ident_string in solvents.keys()):
            self.sld = None;
            self.liq = Solvent(ident_string);
            self.gas = SolventGas(ident_string);
        else:
            if(ident_string in crystals.keys()):
                self.liq = None;
                self.sld = Crystal(ident_string);
                self.gas= None;
            else:
                self.liq = None;
                self.sld = None;
                self.gas = InertGas(ident_string);
    
    def isSolute(self)->bool:
        return isinstance(self.sld,Crystal);


