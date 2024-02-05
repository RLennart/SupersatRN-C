# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 08:53:10 2023

@author: Simon Ternes

    Supersaturation Rate Numerical Calculator (correlations module) - Calculates the dynamic thickness and supersaturation of perovskite solution films 
    
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
import math;
import numpy as np;
import scipy.interpolate as interp;
import pylab as plt;
from supersatrnc.materials import Material,Compound;

class Correlation:
    
    def __init__(self,presentMat:Material,x_mm:float,h_mm:float,D_mm:float,u0_mps:float) -> None:
        self.pm = presentMat;
        self.x_m = x_mm/1000;
        self.h_m = h_mm/1000;
        self.D_m = D_mm/1000;
        self.u0_mps = u0_mps;
        
    def beta(self,T_K:float,p_Pa:float,diffMat:Material) -> float:
        #print("sh_per_length,diff",self.sherwood_per_length(T_K,p_Pa,diffMat),self.pm.diff_m2ps(T_K,p_Pa,diffMat),self.sherwood_per_length(T_K,p_Pa,diffMat)*self.pm.diff_m2ps(T_K,p_Pa,diffMat))
        return self.sherwood_per_length(T_K,p_Pa,diffMat)*self.pm.diff_m2ps(T_K,p_Pa,diffMat)

class RoundJetOnRotatingDisk(Correlation):
    def __init__(self,presentMat:Material,x_mm:float,h_mm:float,D_mm:float,u0_mps:float,omega_1ps:float) -> None:
        super().__init__(presentMat,x_mm,h_mm,D_mm,u0_mps);
        self.omega_1ps = omega_1ps;

    def sherwood_per_length(self,T_K:float,p_Pa:float,diffMat:Material) -> float:
        x = self.x_m;
        D = self.D_m;
        h = self.h_m;
        u0 = self.u0_mps;
        #Fundamentals of H&M page 414  
        #print("Re,Sc",self.pm.reynolds(T_K,p_Pa,u0,x),self.pm.schmidt(T_K,p_Pa,diffMat))
        return 0.083*(self.pm.reynolds(T_K,p_Pa,u0,D))**0.2895*(self.pm.reynolds(T_K,p_Pa,self.omega_1ps,x**2))**0.3006*self.pm.schmidt(T_K,p_Pa,diffMat)**(0.2663)/x; # x or D in corrleation    


class LaminarAverage(Correlation):
    def __init__(self,presentMat:Material,x_mm:float,h_mm:float,D_mm:float,u0_mps:float) -> None:
        if(not x_mm == D_mm):
            raise Exception("For the linear Average correlation, x and D should be equal!");
        super().__init__(presentMat,x_mm,h_mm,D_mm,u0_mps);

    
    def sherwood_per_length(self,T_K:float,p_Pa:float,diffMat:Material) -> float:
        x = self.x_m;
        D = self.D_m;
        h = self.h_m;
        u0 = self.u0_mps;
        #Fundamentals of H&M page 414  
        #print("Re,Sc",self.pm.reynolds(T_K,p_Pa,u0,x),self.pm.schmidt(T_K,p_Pa,diffMat))
        return 0.680*(self.pm.reynolds(T_K,p_Pa,u0,x))**0.5*self.pm.schmidt(T_K,p_Pa,diffMat)**(1/3)/x; # x or D in corrleation        



class LiquidRoundJet(Correlation):
    def __init__(self,presentMat:Material,x_mm:float,h_mm:float,D_mm:float,u0_mps:float) -> None:
        super().__init__(presentMat,x_mm,h_mm,D_mm,u0_mps);
    
    
    def sherwood_per_length(self,T_K:float,p_Pa:float,diffMat:Material) -> float:
        x = self.x_m;
        D = self.D_m;
        h = self.h_m;
        u0 = self.u0_mps;
        #Heat Transfer by Impingement of Circular Free-Surface Liquid Jets
        return 1.648*self.pm.reynolds(T_K,p_Pa,u0,x)**(0.5)*self.pm.schmidt(T_K,p_Pa,diffMat)**0.361/x;



class SpatiallyResolvedCorrelation(Correlation):
    
    
    def __init__(self,presentMat:Material,x_mm:float,h_mm:float,D_mm:float,u0_mps:float) -> None:
        super().__init__(presentMat,x_mm,h_mm,D_mm,u0_mps);
        self.oldParameters=(-1,-1,-1,-1,-1)
        self.spatial_range_D = 500;
        self.interpXrange = np.linspace(-self.D_m*self.spatial_range_D,self.D_m*self.spatial_range_D,100000);
    
    
    def fill_nan(self,A:np.array) -> None:
        inds = np.arange(A.shape[0])
        good = np.where(np.isfinite(A))
        f = interp.interp1d(inds[good], A[good],bounds_error=False,kind='cubic')
        B = np.where(np.isfinite(A),A,f(inds))
        return B
    
    def sherwood_per_length(self,T_K:float,p_Pa:float,diffMat:Material) -> float:
    
        if(np.abs(self.x_m)/self.D_m>self.spatial_range_D):  
            return self.sherwood_unsteady(T_K,p_Pa,diffMat)/np.abs(self.x_m);
    
        if(self.oldParameters==(T_K,p_Pa,self.D_m,self.h_m,self.u0_mps)):
            return self.stored_sherwood_per_length[np.argmin(np.abs(self.x_m-self.interpXrange))];
        
        self.stored_sherwood_per_length = np.zeros(len(self.interpXrange))
        
        for i,x in enumerate(self.interpXrange):
            self.stored_sherwood_per_length[i] = self.sherwood_unsteady(T_K,p_Pa,diffMat,x)/np.abs(x);
 
        self.stored_sherwood_per_length = self.fill_nan(self.stored_sherwood_per_length);
        
        self.oldParameters = (T_K,p_Pa,self.D_m,self.h_m,self.u0_mps);
    
        return self.sherwood_per_length(T_K, p_Pa,diffMat ); #Ending on the second condition.
        

    
class AngleImpingingSlotJet(SpatiallyResolvedCorrelation):

    def __init__(self,presentMat:Material,x_mm:float,h_mm:float,D_mm:float,u0_mps:float,angle_deg:float=90) -> None:
            super().__init__(presentMat,x_mm,h_mm,D_mm,u0_mps);
            self.theta = angle_deg/360*math.pi*2;
            self.params = [1.36,-0.32,0.10,0.44,-0.41,0.5,-0.19]
            self.params = self.params + [0.036,0.67,-0.20,0.42,0.061,-0.11,-0.28,0.69,-0.79,0.22]
            self.params = self.params + [-0.032,0.22,0.78,-0.12,-0.26,-0.076,0.22,-0.32,0.099]
            self.counter = 0;
    
    def sherwood_unsteady(self,T_K:float,p_Pa:float,diffMat:Material,x:float=None) -> float:
           D = self.D_m;
           h = self.h_m;
           u0 = self.u0_mps;
           if(x==None):
               x=self.x_m;
           a,alpha,beta,gamma = self.parametrization(x);
           return a *self.pm.reynolds(T_K,p_Pa,u0,x)**alpha*(np.abs(x)/D)**beta*(h/D)**gamma*self.pm.schmidt(T_K,p_Pa,diffMat)**(1/3);
        
        
    def parametrization(self,x:float) -> (float,float,float,float):
            
            params = self.params;
            
            #invert angle if moving away from center in opposite direction
            if(x<0 ):
                th = math.pi - self.theta;
            else:
                th = self.theta;
            
            
            if(np.abs(x)/self.D_m <.25):
                #Stagnation Point
                
                #Interpolate Stagnation point to avoid overfitting
                return np.nan,np.nan,np.nan,np.nan;    
            
                #not in use, unreachable
                return params[0]+params[1]*th+params[2]*th**2,params[3],0,params[4]+params[5]*th+params[6]*th**2;
                
            if(np.abs(x)/self.D_m>  1 and np.abs(x)/self.D_m<4):
                #Laminar Boundary Layer
                return params[7]+params[8]*th+params[9]*th**2,params[10]+params[11]*th, params[12]+params[13]*th,params[14]+params[15]*th+params[16]*th**2;
            if(np.abs(x)/self.D_m>8):
                #Turbulent Wall Jet
                return params[17]+params[18]*th,params[19]+params[20]*th,params[21]+params[22]*th,params[23]+params[24]*th+params[25]*th**2;
            else:
                #Interpolate between the regimes
                return np.nan,np.nan,np.nan,np.nan;    

class AngleImpingingSlotJetRefitted(AngleImpingingSlotJet):            

    
    def __init__(self,presentMat:Material,x_mm:float,h_mm:float,D_mm:float,u0_mps:float,angle_deg:float=90) -> None:
        super().__init__(presentMat,x_mm,h_mm,D_mm,u0_mps,angle_deg);
        resultOfFit = [2.50994990,-2.38237432,0.59908934,1.21436625,-0.25380268,-0.64188205,0.12061153,
                       -4.44513444,3.58226795,-0.71976560,0.06648731,0.03230517,0.84000932,-0.06024633,
                       -0.96799497,0.23788380,0.69343653,-0.88080381,0.25370242]
                
        for i,p in enumerate(resultOfFit):
            self.params[i+7] = p;
            

class AverageSlotJet(Correlation):
    def __init__(self,presentMat:Material,x_mm:float,h_mm:float,D_mm:float,u0_mps:float) -> None:
        super().__init__(presentMat,x_mm,h_mm,2*D_mm,u0_mps); #D is two times the Slot width
    
    def sherwood_per_length(self,T_K:float,p_Pa:float,diffMat:Material) -> float:
        x = self.x_m;
        D = self.D_m;
        h = self.h_m;
        u0 = self.u0_mps;
        m = 0.695-1/(x/D+(h/D)**(1.33)+3.06)
        return (1.53*(self.pm.reynolds(T_K,p_Pa,u0,D))**m/(x/D+h/D+1.39)*self.pm.schmidt(T_K,p_Pa,diffMat)**0.42)/x;
            
class AverageRoundJet(Correlation):
    def __init__(self,presentMat:Material,x_mm:float,h_mm:float,D_mm:float,u0_mps:float) -> None:
        super().__init__(presentMat,x_mm,h_mm,D_mm,u0_mps);
    
    def sherwood_per_length(self,T_K:float,p_Pa:float,diffMat:Material) -> float:
        x = self.x_m;
        D = self.D_m;
        h = self.h_m;
        u0 = self.u0_mps;
        Re=self.pm.reynolds(T_K,p_Pa,u0,D)
        F = 2*(Re*(1+0.005*Re**0.55))**0.5
        return (1-1.1/(x/D))/((x/D)+0.1*(h/D-6))*F*self.pm.schmidt(T_K,p_Pa,diffMat)**0.4/x;
            
        
