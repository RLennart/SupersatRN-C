import math;
import numpy as np;
import scipy.interpolate as interp;
import copy;
import numpy.typing as npt;
from scipy.integrate import solve_ivp;
from matplotlib.axes import Axes;
from matplotlib.collections import PolyCollection;
from overrides import override
import matplotlib.pyplot as plt


# from supersatrnc.materials import Material
# from supersatrnc.materials import Solvent
# from supersatrnc.materials import Compound
# from supersatrnc.materials import Crystal

#Commonly used constants
R_JpmolK = 8.3144598;
NA_1pmol = 6.02214 * 10**23;




solvents = dict();



solvents['DMSO'] = {'massDens_gpcm3':1.1004, #https://en.wikipedia.org/wiki/Dimethyl_sulfoxide
                    'molMass_gpmol':78.133, #http://webbook.nist.gov/cgi/cbook.cgi?ID=C67685&Mask=4&Type=ANTOINE&Plot=on
                    'diffVol':(15.9 * 2 + 2.31 * 6 + 6.11 + 22.9), #2*C + 6*H+1*O+1*S WLI, Da28
                    'antoineParams':[[[293-30,323],[5.23039,2239.161,-29.215],True],[[323.5,442],[4.49107,1807.002,-60.995],True]], #[0]:tempRange, [1]:A,B,C [2]unit is bar (mmHg) http://webbook.nist.gov/cgi/cbook.cgi?ID=C67685&Mask=4&Type=ANTOINE&Plot=on
                    'surfaceTensionParams':[43.54,-0.1145], #https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
                    'dynViscosity_mPas':2.271, #https://mdpi-res.com/d_attachment/applsci/applsci-12-00116/article_deploy/applsci-12-00116-v2.pdf?version=1640332022
                    'isAntiSolvent':False
                    }

solvents['DMF'] = {'massDens_gpcm3':0.948 , #https://en.wikipedia.org/wiki/Dimethylformamide
                    'molMass_gpmol':73.0938, #http://webbook.nist.gov/cgi/cbook.cgi?ID=C68122&Mask=4
                    'diffVol':(15.9 * 3 + 2.31 * 7 + 6.11 + 4.54 ), #3*C + 7*H+1*O+1*N WLI, Da28
                    'antoineParams':[[[303-30,363],[3.93068,1337.716,-82.648],True]], #http://webbook.nist.gov/cgi/cbook.cgi?ID=C68122&Mask=4
                    'surfaceTensionParams':[37.1,-0.14],  #https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
                    'dynViscosity_mPas':0.92,#https://en.wikipedia.org/wiki/Dimethylformamide
                    'isAntiSolvent':False
                    }

solvents['GBL'] = {'massDens_gpcm3':1.1296 , #https://en.wikipedia.org/wiki/Gamma-Butyrolactone
                    'molMass_gpmol':86.090, #https://en.wikipedia.org/wiki/Gamma-Butyrolactone
                    'diffVol':(15.9 * 4 + 2.31 * 6 + 6.11*2-18.3), #4*C + 6*H+2*O-1*Aromatic Ring, WLI, Da28
                    'antoineParams':[[[-43,465.85],[7.67415,2147.61,244.041],False]], #http://webbook.nist.gov/cgi/cbook.cgi?ID=C68122&Mask=4
                    'surfaceTensionParams':[37.1,-0.14], #https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
                    'dynViscosity_mPas':1.7, #https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/celc.201700138
                    'isAntiSolvent':False
                    }

solvents['2ME'] = {'massDens_gpcm3':0.965 , #https://en.wikipedia.org/wiki/2-Methoxyethanol
                    'molMass_gpmol':76.09, #https://en.wikipedia.org/wiki/2-Methoxyethanol
                    'diffVol':(15.9 * 3 + 2.31 * 8 + 6.11*2), #2*C + 8*H+2*O
                    'antoineParams':[[[-43,465.85],[5.06386,1853.556,-30.838],True]], #https://webbook.nist.gov/cgi/cbook.cgi?ID=C109864&Mask=4&Type=ANTOINE&Plot=on
                    'surfaceTensionParams':[30.84,-0.08], #https://materials.springer.com/lb/docs/sm_lbs_978-3-540-75508-1_53
                    'dynViscosity_mPas':1.7, #https://www.carlroth.com/medias/SDB-8892-AU-EN.pdf?context=bWFzdGVyfHNlY3VyaXR5RGF0YXNoZWV0c3wyNjg5MTB8YXBwbGljYXRpb24vcGRmfGg2OC9oMGMvOTE0ODEzMjQyNTc1OC9TREJfODg5Ml9BVV9FTi5wZGZ8NGM5OGRiMzJhYmJlNzhjNGFiYTU0ZjRkZmI3NTI3YjBjYmIzNzg0YzRmY2E2YmJlMDMyMzlhZTFkYmNlNmI3MA
                    'isAntiSolvent':False
                    }

solvents['CB'] = {'massDens_gpcm3': 1.11 , #https://en.wikipedia.org/wiki/Chlorobenzene
                    'molMass_gpmol': 	112.56, #https://en.wikipedia.org/wiki/Gamma-Butyrolactone
                    'diffVol':(15.9 * 6+2.31 * 5 + 21*1-18.3), #6*C + 5*H+1*Cl-1*Aromatic Ring, WLI, Da28
                    'antoineParams':[[[335-60,405],[4.11083,1435.675,-55.124],True]], #https://vacuu-lan.com/vapor-pressure-estimates-antoine-equation/
                    'surfaceTensionParams':[33.6,-0.1191], #https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
                    'dynViscosity_mPas':0.753, #https://en.wikipedia.org/wiki/Chlorobenzene_(data_page)
                    'isAntiSolvent':True
                    }


solvents['VARCB'] = {'massDens_gpcm3': 1.11 , #https://en.wikipedia.org/wiki/Chlorobenzene
                    'molMass_gpmol': 	112.56, #https://en.wikipedia.org/wiki/Gamma-Butyrolactone
                    'diffVol':(15.9 * 6+2.31 * 5 + 21*1-18.3), #6*C + 5*H+1*Cl-1*Aromatic Ring, WLI, Da28
                    'antoineParams':[[[335-60,405],[4.11083,1435.675,-55.124],True]], #https://vacuu-lan.com/vapor-pressure-estimates-antoine-equation/
                    'surfaceTensionParams':[33.6,-0.1191], #https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
                    'dynViscosity_mPas':0.753, #https://en.wikipedia.org/wiki/Chlorobenzene_(data_page)
                    'isAntiSolvent':True
                    }

dryingGases = dict(); #'RspecGas, viscGas, diffVolGas, mMassGas, thermalConduct, specifc heat

dryingGases['DRY_AIR'] = {'dynViscSutherParams':[120,291.15,18.27*10**(-6)], #https://de.wikipedia.org/wiki/Sutherland-Modell
                          'diffVol':19.7, # Da28
                          'molMass_gpmol':28.9647, #https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
                          'thermalConductivity': 26.2, #https://de.wikipedia.org/wiki/W%C3%A4rmeleitf%C3%A4higkeit
                          'specHeatCap': 1006 #https://www.engineeringtoolbox.com/air-properties-d_156.html
    }# g/mol

dryingGases['N2'] = {'dynViscSutherParams':[111,300.55,18.27*10**(-6)], #https://de.wikipedia.org/wiki/Sutherland-Modell
                     'diffVol':18.5, # Da28
                     'molMass_gpmol':14, #https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
                     'thermalConductivity': 26.0, #https://de.wikipedia.org/wiki/W%C3%A4rmeleitf%C3%A4higkeit
                     'specHeatCap': 1040 #https://www.nuclear-power.com/nitrogen-specific-heat-latent-heat-vaporization-fusion/
    }

crystals = dict(); # Temperatures, Cubic lattice constant in A

#https://www.nature.com/articles/srep35685/figures/6
crystals['MAPI'] = {'latticeConstant_A':6.25, #https://materials.hybrid3.duke.edu/materials/24
                    'molMass_gpmol': (67.52 + 207.2 + 126.9*3),
                    'eqConcentration_molpm3':{'DMF':1520,'DMSO':1520,'GBL':1520,'2ME':800}, #Si section 1
                    'critConcentration_molpm3':{'DMF':6340,'DMSO':6340,'GBL':6340,'2ME':3040}, #Si section 1
                    'quenchingEff_molpm3':{'CB':{'DMF':2.88220474*1000,'DMSO':2.88220474*1000,'GBL':2.88220474*1000,'2ME':2.88220474*1000},
                                           'VARCB':{'DMF':2.88220474*1000,'DMSO':2.88220474*1000,'GBL':2.88220474*1000}}} #Si section

crystals['TWOCAT'] = {'latticeConstant_A':6.36,#https://pubs.acs.org/doi/pdf/10.1021/acs.jpclett.5b01432 for FAPI
                      'molMass_gpmol':(132.91*0.17 + 171.97*0.83 + 207.2 + 126.9*3*0.91 + 79.9*3*0.09),
                      'eqConcentration_molpm3':{'DMF':1520,'DMSO':1520,'GBL':1520}, #Si section 1
                      'critConcentration_molpm3':{'DMF':6016,'DMSO':6016,'GBL':6016}, #Si section 1
                      'quenchingEff_molpm3':{'CB':{'DMF':3.64457431*1000,'DMSO':3.64457431*1000,'GBL':3.64457431*1000},'VARCB':{'DMF':2.88220474*1000,'DMSO':2.88220474*1000,'GBL':2.88220474*1000}}} #Si section





class Material:

    def __init__(self, ident_string: str) -> None:
        self.ident_string = ident_string;

    def __str__(self) -> str:
        return self.ident_string;

    def kinViscosity_m2ps(self, T_K: float, p_Pa: float) -> float:
        return self.dynViscosity_Pas(T_K, p_Pa) / self.massDens_kgpm3(T_K, p_Pa);

    def molDens_molpm3(self, T_K: float, p_Pa: float) -> float:
        raise AttributeError("Density-Function not implemented. Must be overwritten.");

    def massDens_kgpm3(self, T_K: float, p_Pa: float) -> float:
        return self.molDens_molpm3(T_K, p_Pa) * self._molMass_kgpmol;

    def specVol_m3pmol(self, T_K: float, p_Pa: float) -> float:
        return 1 / self.molDens_molpm3(T_K, p_Pa);

    def dynViscosity_Pas(self, T_K: float, p_Pa: float) -> float:
        raise AttributeError("Viscosity-Function not implemented.");


class Fluid(Material):

    def __init__(self, ident_string: str) -> None:
        super().__init__(ident_string);

    def schmidt(self, T_K: float, p_Pa: float, inMaterial: Material) -> float:
        return self.kinViscosity_m2ps(T_K, p_Pa) / self.diff_m2ps(T_K, p_Pa, inMaterial);

    def reynolds(self, T_K: float, p_Pa: float, vel_mps: float, width_m: float) -> float:
        return vel_mps * np.abs(width_m) / self.kinViscosity_m2ps(T_K, p_Pa);


class IdealGas(Fluid):

    def __init__(self, ident_string: str) -> None:
        super().__init__(ident_string);

    @override
    def molDens_molpm3(self, T_K: float, p_Pa: float) -> float:
        return p_Pa / R_JpmolK / T_K;

    # Fuller-Gas Diffusion
    def diff_m2ps(self, T_K: float, p_Pa: float, gas: 'IdealGas') -> float:
        return 0.01 * 0.01 * 0.00143 * (T_K ** 1.75) * (
                    (self._molMass_kgpmol * 1000) ** (-1) + (gas._molMass_kgpmol * 1000) ** (-1)) ** 0.5 / (
                    p_Pa * 1e-5 * np.sqrt(2) * (self._diffVol ** (1.0 / 3) + gas._diffVol ** (1.0 / 3)) ** 2)


class InertGas(IdealGas):

    def __init__(self, ident_string: str) -> None:
        super().__init__(ident_string);
        gas = dryingGases[ident_string];
        self._molMass_kgpmol = gas['molMass_gpmol'] / 1000;
        self._diffVol = gas['diffVol'];
        self._sutherCoeffs = gas['dynViscSutherParams'];

    # Sutherland Equation for dynamic viscosity -> https://www.cfd-online.com/Wiki/Sutherland%27s_law
    @override
    def dynViscosity_Pas(self, T_K: float, p_Pa: float) -> float:  # Sutherland
        C = self._sutherCoeffs[0];
        refT = self._sutherCoeffs[1];
        refDynVisc = self._sutherCoeffs[2];

        return refDynVisc * (refT + C) / (T_K + C) * (T_K / refT) ** (1.5);


class SolventGas(IdealGas):

    def __init__(self, ident_string: str) -> None:
        super().__init__(ident_string);
        solv = solvents[ident_string];
        self._molMass_kgpmol = solv['molMass_gpmol'] / 1000;
        self._diffVol = solv['diffVol'];
        self._antoineParams = solv['antoineParams']

    # Antoine-Equation https://en.wikipedia.org/wiki/Antoine_equation
    def vaporPressure_Pa(self, T_K: float) -> float:
        for tempSet in self._antoineParams:
            A = tempSet[1][0];
            B = tempSet[1][1];
            C = tempSet[1][2];
            if (tempSet[2]):
                if (tempSet[0][0] < T_K and tempSet[0][1] > T_K):
                    return np.power(10, A - B / (T_K + C)) * 10 ** 5;  # from bar to pascal
            else:  # Old Units °C and mmHg
                if (tempSet[0][0] < T_K - 273.15 and tempSet[0][1] > T_K - 273.15):
                    return np.power(10, A - B / (T_K - 273.15 + C)) * 133.322;  # from mmHg to pascal


class Solvent(Fluid):

    def __init__(self, ident_string: str) -> None:
        super().__init__(ident_string);
        solv = solvents[ident_string];
        self._molMass_kgpmol = solv['molMass_gpmol'] / 1000;
        self._isAntiSolvent = solv['isAntiSolvent'];
        self._molDens_kgpm3 = solv['massDens_gpcm3'] * 1000
        self._dynViscosity_Pas = solv['dynViscosity_mPas'] / 1000
        self._surfacetensionParams = solv['surfaceTensionParams']

    def isAntiSolvent(self) -> bool:
        return self._isAntiSolvent;

    @override
    def molDens_molpm3(self, T_K: float, p_Pa: float) -> float:
        return self._molDens_kgpm3 / self._molMass_kgpmol

    # Tyn-Clausius Liquid-Liquid Diffusion Wärematlas Da28
    def diff_m2ps(self, T_K: float, p_Pa: float, liq: 'Solvent') -> float:
        return 0.01 * 0.01 * 8.93e-8 * (self.specVol_m3pmol(T_K, p_Pa) * 100 * 100 * 100) ** (-1 / 3) * (
                    liq.specVol_m3pmol(T_K, p_Pa) * 100 * 100 * 100) ** (1 / 6) * (
                    self.paraChor_m3pmol(T_K, p_Pa) / liq.paraChor_m3pmol(T_K, p_Pa)) ** (0.6) * T_K * (
                    self.dynViscosity_Pas(T_K, p_Pa) * 1000) ** (-1)

    @override
    def dynViscosity_Pas(self, T_K: float, p_Pa: float) -> float:
        return self._dynViscosity_Pas;

    # Linear function of temperature, Reference at 20°C https://www.dataphysics-instruments.com/Downloads/Surface-Tensions-Energies.pdf
    def surfaceTension_Npm(self, T_K: float, p_Pa: float) -> float:
        params = self._surfacetensionParams;
        return (params[0] + params[1] * (T_K - 293)) / 1000;

    # Wärematlas Da 28
    def paraChor_m3pmol(self, T_K: float, p_Pa: float) -> float:
        return self.specVol_m3pmol(T_K, p_Pa) * 100 * 100 * 100 * (self.surfaceTension_Npm(T_K, p_Pa) * 1000) * (0.25)


class Crystal(Material):
    def __init__(self, ident_string: str) -> None:
        super().__init__(ident_string);
        cry = crystals[ident_string]

        self.ceq = cry['eqConcentration_molpm3']

        self.ccrit = cry['critConcentration_molpm3'];

        self.anti_qs = cry['quenchingEff_molpm3'];

        self._latticeConstant_m = cry['latticeConstant_A'] * 1e-10;

        self._molMass_kgpmol = cry['molMass_gpmol'] / 1000;

    @override
    def molDens_molpm3(self, T_K: float, p_Pa: float) -> float:
        return 1 / NA_1pmol / (self._latticeConstant_m) ** 3

    def massDens_kgpm3(self, T_K: float, p_Pa: float) -> float:
        return self.molDens_molpm3(T_K, p_Pa) * self._molMass_kgpmol;


class Compound:

    def __eq__(self, other) -> bool:
        if isinstance(other, str):  # enabling equality just by comparing ident_string
            return self.ident_string == other;
        if isinstance(other, self.__class__):
            return self.ident_string == other.ident_string;
        else:
            return False

    def __hash__(self) -> float:
        return hash(self.ident_string);

    def __ne__(self, other) -> bool:
        return not self.__eq__(other);

    def __str__(self) -> str:
        return "Compound: " + self.ident_string

    def __init__(self, ident_string: str, color: str) -> None:
        self.color = color;
        self.ident_string = ident_string;
        if (ident_string in solvents.keys()):
            self.sld = None;
            self.liq = Solvent(ident_string);
            self.gas = SolventGas(ident_string);
        else:
            if (ident_string in crystals.keys()):
                self.liq = None;
                self.sld = Crystal(ident_string);
                self.gas = None;
            else:
                self.liq = None;
                self.sld = None;
                self.gas = InertGas(ident_string);

    def isSolute(self) -> bool:
        return isinstance(self.sld, Crystal);


MAPI = Compound('MAPI', 'sienna');
TWOCAT = Compound('TWOCAT', 'sienna');
DMF = Compound('DMF', 'lightskyblue');
DMSO = Compound('DMSO', 'gold');
GBL = Compound('GBL', 'palevioletred');
DRY_AIR = Compound('DRY_AIR', 'white');
N2 = Compound('N2', 'white');
CB = Compound('CB', 'palegreen');
VARCB = Compound('VARCB', 'palegreen');
TWOME = Compound('2ME', 'seagreen');

class Correlation:

    def __init__(self, presentMat: Material, x_mm: float, h_mm: float, D_mm: float, u0_mps: float) -> None:
        self.pm = presentMat;
        self.x_m = x_mm / 1000;
        self.h_m = h_mm / 1000;
        self.D_m = D_mm / 1000;
        self.u0_mps = u0_mps;

    def beta(self, T_K: float, p_Pa: float, diffMat: Material) -> float:
        # print("sh_per_length,diff",self.sherwood_per_length(T_K,p_Pa,diffMat),self.pm.diff_m2ps(T_K,p_Pa,diffMat),self.sherwood_per_length(T_K,p_Pa,diffMat)*self.pm.diff_m2ps(T_K,p_Pa,diffMat))
        return self.sherwood_per_length(T_K, p_Pa, diffMat) * self.pm.diff_m2ps(T_K, p_Pa, diffMat)


class SpatiallyResolvedCorrelation(Correlation):

    def __init__(self, presentMat: Material, x_mm: float, h_mm: float, D_mm: float, u0_mps: float) -> None:
        super().__init__(presentMat, x_mm, h_mm, D_mm, u0_mps);
        self.oldParameters = (-1, -1, -1, -1, -1)
        self.spatial_range_D = 500;
        self.interpXrange = np.linspace(-self.D_m * self.spatial_range_D, self.D_m * self.spatial_range_D, 100000);

    def fill_nan(self, A: np.array) -> None:
        inds = np.arange(A.shape[0])
        good = np.where(np.isfinite(A))
        f = interp.interp1d(inds[good], A[good], bounds_error=False, kind='cubic')
        B = np.where(np.isfinite(A), A, f(inds))
        return B

    def sherwood_per_length(self, T_K: float, p_Pa: float, diffMat: Material) -> float:

        if (np.abs(self.x_m) / self.D_m > self.spatial_range_D):
            return self.sherwood_unsteady(T_K, p_Pa, diffMat) / np.abs(self.x_m);

        if (self.oldParameters == (T_K, p_Pa, self.D_m, self.h_m, self.u0_mps)):
            return self.stored_sherwood_per_length[np.argmin(np.abs(self.x_m - self.interpXrange))];

        self.stored_sherwood_per_length = np.zeros(len(self.interpXrange))

        for i, x in enumerate(self.interpXrange):
            self.stored_sherwood_per_length[i] = self.sherwood_unsteady(T_K, p_Pa, diffMat, x) / np.abs(x);

        self.stored_sherwood_per_length = self.fill_nan(self.stored_sherwood_per_length);

        self.oldParameters = (T_K, p_Pa, self.D_m, self.h_m, self.u0_mps);

        return self.sherwood_per_length(T_K, p_Pa, diffMat);  # Ending on the second condition.


class AngleImpingingSlotJet(SpatiallyResolvedCorrelation):

    def __init__(self, presentMat: Material, x_mm: float, h_mm: float, D_mm: float, u0_mps: float,
                 angle_deg: float = 90) -> None:
        super().__init__(presentMat, x_mm, h_mm, D_mm, u0_mps);
        self.theta = angle_deg / 360 * math.pi * 2;
        self.params = [1.36, -0.32, 0.10, 0.44, -0.41, 0.5, -0.19]
        self.params = self.params + [0.036, 0.67, -0.20, 0.42, 0.061, -0.11, -0.28, 0.69, -0.79, 0.22]
        self.params = self.params + [-0.032, 0.22, 0.78, -0.12, -0.26, -0.076, 0.22, -0.32, 0.099]
        self.counter = 0;

    def sherwood_unsteady(self, T_K: float, p_Pa: float, diffMat: Material, x: float = None) -> float:
        D = self.D_m;
        h = self.h_m;
        u0 = self.u0_mps;
        if (x == None):
            x = self.x_m;
        a, alpha, beta, gamma = self.parametrization(x);
        return a * self.pm.reynolds(T_K, p_Pa, u0, x) ** alpha * (np.abs(x) / D) ** beta * (
                    h / D) ** gamma * self.pm.schmidt(T_K, p_Pa, diffMat) ** (1 / 3);

    def parametrization(self, x: float) -> (float, float, float, float):

        params = self.params;

        # invert angle if moving away from center in opposite direction
        if (x < 0):
            th = math.pi - self.theta;
        else:
            th = self.theta;

        if (np.abs(x) / self.D_m < .25):
            # Stagnation Point

            # Interpolate Stagnation point to avoid overfitting
            return np.nan, np.nan, np.nan, np.nan;

            # not in use, unreachable
            return params[0] + params[1] * th + params[2] * th ** 2, params[3], 0, params[4] + params[5] * th + params[
                6] * th ** 2;

        if (np.abs(x) / self.D_m > 1 and np.abs(x) / self.D_m < 4):
            # Laminar Boundary Layer
            return params[7] + params[8] * th + params[9] * th ** 2, params[10] + params[11] * th, params[12] + params[
                13] * th, params[14] + params[15] * th + params[16] * th ** 2;
        if (np.abs(x) / self.D_m > 8):
            # Turbulent Wall Jet
            return params[17] + params[18] * th, params[19] + params[20] * th, params[21] + params[22] * th, params[
                23] + params[24] * th + params[25] * th ** 2;
        else:
            # Interpolate between the regimes
            return np.nan, np.nan, np.nan, np.nan;



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

    def __init__(self, sp: SimulationParameters, compounds: list[Compound], d_init_mum: float,
                 volume_ratios_solvent: list[float], d_final_mum: float,
                 volume_ratios_solute: list[float] = None) -> None:

        self.sp = sp;

        'Key lists defining the order'
        self.slvs = list();
        self.slts = list();
        self.only_slvs = list();
        self.only_aslvs = list();
        self.compounds = list();

        if (volume_ratios_solute == None):
            volume_ratios_solute = [1 / ((len(compounds) - len(volume_ratios_solvent)))] * (
                        len(compounds) - len(volume_ratios_solvent));

        self.Ni0s_mol = dict();
        self.Np0s_mol = dict();

        # self.__di0s_m = dict();
        # self.__dp0s_m = dict();

        self.solution = dict();

        self.dp_m = d_final_mum / 1000 / 1000;
        self.d0_m = d_init_mum / 1000 / 1000 - self.dp_m;

        for i, c in enumerate(compounds):
            if (not c.isSolute()):
                r = volume_ratios_solvent.pop(0);
            else:
                r = volume_ratios_solute.pop(0);

            self.addCompound(c, r);

        # list of compounds (e.g. DMF and MAPI)
        self.compounds = self.slts + self.slvs;

        self.Np_mol = np.sum(np.array(list(self.Np0s_mol.values())));

    # this function should be defined in film setup class, logic from above should be in here
    def addCompound(self, c: Compound, init_ratio: float) -> None:
        if (not c.isSolute()):
            self.slvs.append(c);

            if (c.liq.isAntiSolvent()):
                self.only_aslvs.append(c);
            else:
                self.only_slvs.append(c);

            self.Ni0s_mol[c] = self.d0_m * init_ratio * c.liq.molDens_molpm3(self.sp.T, self.sp.p0);

        else:
            self.slts.append(c);
            self.Np0s_mol[c] = self.dp_m * init_ratio * c.sld.molDens_molpm3(self.sp.T, self.sp.p0);

    def Ni0s(self) -> npt.NDArray[np.float64]:
        return np.array([self.Ni0s_mol[slv] for slv in self.slvs]);

    def di0_m(self, slv: Solvent) -> npt.NDArray[np.float64]:
        return self.Ni0s_mol[slv] / slv.liq.molDens_molpm3(self.sp.T, self.sp.p0);

    def di0s(self) -> npt.NDArray[np.float64]:
        return np.array([self.di0_m(slv) for slv in self.slvs]);

    def Npj0s(self) -> npt.NDArray[np.float64]:
        return np.array([self.Np0s_mol[slt] for slt in self.slts]);

    def dpj0_m(self, slt: Crystal) -> npt.NDArray[np.float64]:
        return self.Np0s_mol[slt] / slt.sld.molDens_molpm3(self.sp.T, self.sp.p0)

    def dpj0s(self) -> npt.NDArray[np.float64]:
        return np.array([self.dpj0_m(slt) for slt in self.slts]);

    def averageMolDensP_molpm3(self) -> npt.NDArray[np.float64]:
        return np.sum(self.dp0s() / self.dp_m * [slt.sld.molDens_molpm3(self.sp.T, self.sp.p0) for slt in self.slts])

    def getFilmWithAddedAntiSolvent(self, antiSolvent: Solvent) -> None:

        ret = copy.deepcopy(self);
        ret.addCompound(antiSolvent, 0);
        ret.compounds = ret.slts + ret.slvs;

        return ret;

    def compoundNames(self) -> list[str]:
        return [k.ident_string for k in self.compounds];


class CalculationResult():

    def __init__(self, film: Film, times: npt.NDArray[np.float64], sol_data: npt.NDArray[np.float64]) -> None:

        self.film = film;
        self.sp = film.sp;
        self.solution = dict();
        self.times = times;

        for i, s in enumerate(sol_data):
            self.solution[film.slvs[i]] = s;

        self.findCritSuperSatAndRate();

    def thicknesses_m(self, liquids: list[Solvent] = [], offset: float = 0) -> npt.NDArray[np.float64]:
        if (liquids == []):
            liquids = self.film.slvs;

        ths = 0;
        for i, s in enumerate(liquids):
            ths = ths + self.solution[s] / s.liq.molDens_molpm3(self.sp.T, self.sp.p0);
        ths = ths + offset;
        return ths;

    def allThicknessesAndColors(self, max_ind: int) -> npt.NDArray[np.float64]:

        ths_sld = [np.array([self.film.dpj0_m(slt) * 1000 * 1000] * (max_ind), dtype=float) for slt in self.film.slts];

        ths_slv = [np.array(self.solution[k][0:max_ind] / k.liq.molDens_molpm3(self.sp.T, self.sp.p0) * 1000 * 1000,
                            dtype=float) for k in self.film.slvs];

        colors = [k.color for k in self.film.compounds]

        return ths_sld + ths_slv, colors;

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

    def c_lin(self, crit: bool = False) -> npt.NDArray[np.float64]:

        lenSol = len(self.solution[self.film.slvs[0]]);

        ret = np.zeros(lenSol)

        for slt in self.film.slts:
            add = np.zeros(lenSol)

            for slv in self.film.only_slvs:
                ### Summer über richtigen INDEX!!
                ri = self.solution[slv] / slv.liq.molDens_molpm3(self.sp.T, self.sp.p0) / self.thicknesses_m(
                    self.film.only_slvs);

                ri[np.where(np.isnan(ri))] = 1;

                if (not crit):
                    add = add + slt.sld.ceq[slv] * ri;
                else:
                    add = add + slt.sld.ccrit[slv] * ri;

                for aslv in self.film.only_aslvs:

                    if (aslv in slt.sld.anti_qs):
                        add = add - ri * slt.sld.anti_qs[aslv][slv] * self.solution[aslv] / aslv.liq.molDens_molpm3(
                            self.sp.T, self.sp.p0) / self.thicknesses_m(self.film.only_slvs)

            ret = ret + add * self.film.dpj0_m(slt) / self.film.dp_m;

            # Concentrations are always greater zero - todo what?? why is this code required at all?
        for i, r in enumerate(ret):
            if (r < 0):
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

        threshold = np.log(crits[0] / ceqs[0])

        # print(self.crit_index<(len(supersats)-1),np.isnan(supersats[self.crit_index]),supersats[self.crit_index] < threshold,supersats[self.crit_index],threshold)

        # todo this is lol
        while (self.crit_index < (len(self.supersats) - 1) and (
                np.isnan(self.supersats[self.crit_index]) or self.supersats[self.crit_index] < threshold)):
            self.crit_index = self.crit_index + 1;

        self.superSatRate = np.gradient(self.supersats, (self.times[1] - self.times[0]));
        self.crit_superSat = self.supersats[self.crit_index];
        self.crit_superSatRate = self.superSatRate[self.crit_index];

    def getFilmWithShiftedInitState(self, index: int = -1) -> Film:

        d_init = self.thicknesses_m(offset=self.film.dp_m)[index] * 1000 * 1000;

        vrs = [self.thicknesses_m([slv])[index] / self.thicknesses_m()[index] for slv in self.film.slvs];

        # print("-------- VRS",vrs)

        vrps = self.film.dpj0s() / self.film.dp_m;

        return Film(self.film.sp, self.film.compounds, d_init, list(vrs), self.film.dp_m * 1000 * 1000, list(vrps));

    # def djs(self) -> npt.NDArray[np.float64]:
    #     return self.film.d0_m-self.totaldisTimesWeights_m();
    #
    # def Vratios(self) -> npt.NDArray[np.float64]:
    #     return self.djs()/self.totaldisTimesWeights_mol();


class DynamicsSimulation:

    def __init__(self, film: Film, end_time_s: float, number_of_time_steps: int = 10000) -> None:
        self.film = film;
        self.end_time_s = end_time_s;
        self.number_of_time_steps = number_of_time_steps;

        if (self.film.sp.t0 > self.end_time_s):
            raise AttributeError("t0 cannot be more than maximum simulation time");
        self.times = np.linspace(self.film.sp.t0, self.end_time_s, self.number_of_time_steps);

        self.result = None;
        self.sv = SimulationVariables(self.film.sp);

    def computeSolution(self) -> None:

        self.sv = SimulationVariables(self.film.sp);

        # this is the expensive line
        sol = solve_ivp(self.Niprimes_molps, [self.times[0], self.times[-1]], self.film.Ni0s(), dense_output=True,
                        max_step=self.times[1] - self.times[0]);

        # this is weird. here the y-solution is extracted and then fed into Calculation Result
        soly = sol.sol(self.times);  # this is likely more-dimensional if multiple solvent components inside

        self.result = CalculationResult(self.film, self.times, soly);

    def Niprimes_molps(self, time: float, N: list[float]) -> list[float]:
        ret = list();

        for i, s in enumerate(self.film.slvs):

            # If Ni is sufficintly small, the drying rate will be zero
            if (N[i] > 1e-10 * self.film.Np_mol):
                ret.append(-self.prop(s, i, N) * N[i]);
            else:
                ret.append(0);

        self.timeStep(time);  # set correct time for time dependent properties;

        return ret;

    def timeStep(self, time: float) -> None:

        self.sv.t = time;

    def incrementIndex(self, max_ind: int) -> int:
        if (max_ind == -1):
            return len(self.times);

        return max_ind + 1;

    def plotTotalFilmThickness(self, ax: Axes, interferometricDataFiles: list[str] = [], time_div: float = 0,
                               max_ind: int = -1) -> None:

        max_ind = self.incrementIndex(max_ind);

        for f in interferometricDataFiles:
            data = np.load(f)

            # append start and end points
            xdata = np.insert(np.append(data[0], self.times[0:max_ind][-1]), 0, self.times[0:max_ind][0])
            ydata = np.insert(np.append(data[1], self.result.film.dp_m * 1000 * 1000), 0,
                              (self.result.film.d0_m + self.result.film.dp_m) * 1000 * 1000)

            ax.plot(xdata + time_div, ydata, 'o', color='black')

        p2, = ax.plot(self.times[0:max_ind] + time_div, self.result.thicknesses_m() * 1000 * 1000);
        ax.set_ylabel("Thickness [$\mu$m]")
        ax.set_xlabel("Time [s]")

    def plotThicknessComposition(self, ax: Axes, max_ind: int, time_div: float) -> list[PolyCollection]:

        ths, colors = self.result.allThicknessesAndColors(max_ind);

        # print("times",len(self.times[0:max_ind]+time_div));
        # print("ths",[len(ths[i]) for i in range(len(ths))]);

        polys = ax.stackplot(np.array(self.times[0:max_ind] + time_div, dtype=float), ths, colors=colors);

        # print(ths,colors)

        return polys;

    def plotFilmThicknessCompositionAndSupersaturation(self, axes: tuple[Axes],
                                                       interferometricDataFiles: list[str] = [],
                                                       dataLabels: list[str] = [], xytext: tuple[float] = (10, -10),
                                                       max_ind: int = -1, time_div: float = 0) -> None:

        max_ind = self.incrementIndex(max_ind);

        firstPlot = True;
        if (len(axes) == 2):
            top, bottom = axes;
        else:
            if (len(axes) == 3):
                firstPlot = False;
                top, bottom, par2 = axes;
            else:
                raise Exception("Pass two Axes subplot or three Axes")

        if (firstPlot):
            top.set_ylabel("Supersat. Rate [s$^{-1}$]")
            bottom.set_xlabel("Time t [s]")
            bottom.set_ylabel("Composition of film th. [$\mu$m]")

        p1, = top.plot(self.times[2:max_ind] + time_div, self.result.superSatRate[2:max_ind],
                       color='firebrick');  # first two points of drivative might be flawed.

        if (self.result.crit_index < max_ind - 1):
            top.axvline(x=self.times[self.result.crit_index] + time_div, linestyle='--')
            top.plot(self.times[self.result.crit_index] + time_div, self.result.crit_superSatRate, 'o',
                     color='firebrick', markersize=15);

            top.annotate("$\\dot{\\sigma}_{crit.}$",
                         (self.times[self.result.crit_index] + time_div, self.result.crit_superSatRate),
                         textcoords="offset points",
                         xytext=xytext, color="firebrick")
            bottom.axvline(x=self.times[self.result.crit_index] + time_div, linestyle='--')

            print("Determined critical supersat. rate at ", self.result.crit_superSatRate, " per second");

        top.set_xlim([0, self.times[-1]])
        bottom.set_xlim([0, self.times[-1]])

        top.yaxis.label.set_color(p1.get_color())
        top.tick_params(axis='y', colors=p1.get_color())

        polys = self.plotThicknessComposition(bottom, max_ind, time_div)

        if (firstPlot):
            par2 = top.twinx()

        # print(self.result.supersaturation()[0:max_ind])

        p2, = par2.plot(self.times[0:max_ind] + time_div, self.result.supersats[0:max_ind], color='cornflowerblue');

        if (firstPlot):
            par2.yaxis.label.set_color(p2.get_color())
            par2.tick_params(axis='y', colors=p2.get_color())
            par2.set_ylabel("Supersaturation")

        add_lines = list();

        for f in interferometricDataFiles:
            data = np.load(f)
            l, = bottom.plot(data[0], data[1], 'o', color='black')
            add_lines.append(l);

        bottom.legend(polys + add_lines, self.result.film.compoundNames() + dataLabels, loc=1)

        ylabels = bottom.get_yticklabels()
        ylabels[0].set_alpha(0)
        top.set_ylim([-0.05 * self.maximumSupersaturation(), self.maximumSupersaturation() * 1.2])

        plt.tight_layout()
        return top, bottom, par2;

    def maximumSupersaturation(self) -> float:
        return np.max(self.result.superSatRate)


class GasQuenchDynamics(DynamicsSimulation):

    def __init__(self, film: Film, end_time_s: float, number_of_time_steps: int = 10000,
                 activities: dict[Compound, float] = dict()) -> None:
        super().__init__(film, end_time_s, number_of_time_steps)
        self.activities = activities;

    def prop(self, m: Compound, i: int, Nis: list[float]) -> float:
        ret = (self.film.sp.corr.beta(self.film.sp.T, self.sv.p, m.gas)
               * m.gas.vaporPressure_Pa(self.film.sp.T) / R_JpmolK / self.film.sp.T / (np.sum(Nis) + self.film.Np_mol));

        # print(self.sv.t,self.film.sp.corr.beta(self.film.sp.T,self.sv.p,m.gas),m.gas.vaporPressure_Pa(self.film.sp.T),self.denom(i, Nis))

        # todo rather set activities to standard 1 if not given
        if (not m in self.activities.keys()):
            return ret

        return ret * self.activities[m];

    def denom(self, i: int, Nis: list[float]) -> float:
        return (np.sum(Nis) + self.film.Np_mol);


class DynamicGasQuenchDynamics(GasQuenchDynamics):

    def __init__(self, film: Film, end_time_s: float, number_of_time_steps: int = 10000,
                 activities: dict[Compound, float] = dict(), v_mps: float = 1, L_m: float = None) -> None:
        super().__init__(film, end_time_s, number_of_time_steps, activities)

        self.L_m = L_m;
        if (L_m == None):
            # it does not jump in here but it would cause errors since not defined
            self.L_m = self.film.corr.D_m * 10;

        self.v_mps = v_mps;

    def timeStep(self, time: float) -> None:
        super().timeStep(time);
        # self.film.sp.corr.x_m = self.v_mps*(time+self.film.sp.t0)-self.L_m;
        self.film.sp.corr.x_m = self.v_mps * (time) - self.L_m;


class VacuumDynamics(GasQuenchDynamics):

    def __init__(self, film: Film, end_time_s: float, number_of_time_steps: int = 10000,
                 activities: dict[Compound, float] = dict(),
                 p_params: list[float] = [0.55143471, 7.74808886, 0.02585276],
                 u0_params: list[float] = [0.09453549, 5.2]):
        super().__init__(film, end_time_s, number_of_time_steps, activities)
        self.pp = p_params;
        self.u0p = u0_params;
        self.setPressureAndVel();

    def timeStep(self, time: float) -> None:
        super().timeStep(time);

        self.setPressureAndVel();

    def setPressureAndVel(self) -> None:
        self.sv.p = self.film.sp.p0 * np.exp(-self.pp[0] * (self.sv.t)) + self.pp[1] - self.pp[2] * (self.sv.t);
        self.film.sp.corr.u0_mps = self.u0p[1] * np.exp(-self.u0p[0] * (self.sv.t));

if __name__ == "__main__":

    # fragen
    # wieso macht ein so kleiner anteil dmso so einen großen unterschied in der optischen wahrnehmung der samples während/nach dem quenching
    #

    corr = AngleImpingingSlotJet(DRY_AIR.gas, 0, 3, 0.1, 100);
    sp = SimulationParameters(corr, 22);

    fig, axs = plt.subplots(nrows=2, figsize=(4,7))
    film1 = Film(sp, [DMF, MAPI], 4.7, [1], 1.18)  # d_pm = d_final_mum
    dynsDynGas = DynamicGasQuenchDynamics(film1, 35, v_mps=.00595, L_m=.15);  # vmps coating speed
    dynsDynGas.computeSolution();
    dynsDynGas.plotFilmThicknessCompositionAndSupersaturation(axs);

    fig, axs = plt.subplots(nrows=2, figsize=(4,7))
    corr = AngleImpingingSlotJet(DRY_AIR.gas, 0, 3, .1, 100);
    sp = SimulationParameters(corr, 22);
    film1 = Film(sp, [GBL, DMSO, DMF, TWOCAT], 7.3, [1 / 3 * 1.3, 0.2 * 2 / 3 * 1.3, 0.8 * 2 / 3 * 0.7], 0.86)
    film2 = Film(sp, [GBL, DMSO, DMF, TWOCAT], 7.3, [1 / 3 * 1.3, 0.2 * 2 / 3 * 1.3, 0.8 * 2 / 3 * 0.7], 0.86)
    dynsDynGas1 = DynamicGasQuenchDynamics(film1, 550, v_mps=.00035, L_m=.15);
    dynsDynGas2 = DynamicGasQuenchDynamics(film2, 550, v_mps=.00035, L_m=.15,
                                           activities={DMF: 1, DMSO: 0.16479906, GBL: 0.11083785});

    dynsDynGas1.computeSolution();
    dynsDynGas1.plotFilmThicknessCompositionAndSupersaturation(axs);

    dynsDynGas2.computeSolution();
    dynsDynGas2.plotFilmThicknessCompositionAndSupersaturation(axs);

    plt.show()
