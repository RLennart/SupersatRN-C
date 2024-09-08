from supersatrnc.__main__ import DryAndWetThicknessFitter,ActivityFitter,SimulationParameters,Film,GasQuenchDynamics,DynamicGasQuenchDynamics,VacuumDynamics,SimulationQueue,ParameterSweeper;
from supersatrnc.materials import Compound,Material;
from supersatrnc.correlations import LaminarAverage,LiquidRoundJet,AngleImpingingSlotJet,AngleImpingingSlotJetRefitted,AverageSlotJet,AverageRoundJet;
##################### Initialize and set specific colors
        
MAPI = Compound('MAPI','sienna');
TWOCAT = Compound('TWOCAT','sienna');
DMF = Compound('DMF','lightskyblue');
DMSO = Compound('DMSO','gold');
GBL = Compound('GBL','palevioletred');
DRY_AIR = Compound('DRY_AIR','white');
N2 = Compound('N2','white');
CB = Compound('CB','palegreen');
VARCB = Compound('VARCB','palegreen');
TWOME = Compound('2ME','seagreen');