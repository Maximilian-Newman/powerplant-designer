import cea
from CoolProp import CoolProp as cp # haha very funny
from CoolProp.Plots import PropertyPlot
import CoolProp
import matplotlib.pyplot as plt
import math
import copy
import numpy as np
import multiprocessing

CEA_TO_COOLPROP_NAMES = {
    "CH4" : "Methane",
    "O2" : "Oxygen",
    "H2" : "Hydrogen",
    "N2" : "Nitrogen",
    "CO2" : "CarbonDioxide",
    "CO" : "CarbonMonoxide",
    "H2O" : "Water",
    "H2O(L)" : "Water",
    "H2O(cr)" : "Water",
    "Ar" : "Argon",
    "Ne" : "Neon",
}

KNOWN_CP_UNSUPPORTED = ["H", "OH", "NO", "C(gr)", "O"]

# having these pairs of molecules in a mixture crashes coolprop:
COOLPROP_INCOMPATIBLE_PAIRS = [
    ["Argon", "Neon"],
    ["Methane", "NH3"],
    ["Nitrogen", "NH3"]
]

def cea_to_cp_name(name):
    if name in KNOWN_CP_UNSUPPORTED:
        return None
    try:
        return CEA_TO_COOLPROP_NAMES[name]
    except:
        try:
            M = cp.PropsSI("molar_mass", name)
            return name
        except:
            print("Couldn't put into coolprop mixture:", name)
            return None

class Mixture:
    def __init__(self, data=None):
        if data == None:
            data = []
        self.data = data

    def to_cea(self, temperature, products_from_reactants=False):
        names = []
        temps = []
        molFracs = []
        massFracs = []
        for i in range(0, len(self.data)):
            names.append(self.data[i][0])
            molFracs.append(self.data[i][1])
            if len(self.data[i]) == 2:
                temps.append(temperature)
            else:
                temps.append(self.data[i][2])
        
        mix = cea.Mixture(names, products_from_reactants=products_from_reactants)
        totMass = 0
        for i in range(0, len(names)):
            M = cp.PropsSI("molar_mass", "T", temperature, "P", 101325, self.to_coolprop())
            m = molFracs[i] * M
            massFracs.append(m)
            totMass += m
        for i in range(0, len(names)):
            massFracs[i] = massFracs[i] / totMass
        
        return [mix, massFracs, temps]

    def to_coolprop(self):
        if self.data == AIR().data:
            return "Air"
        
        filtered = Mixture()
        mix = "HEOS::"
        self.sort()
        for part in self.data:
            cpName = cea_to_cp_name(part[0])
            if cpName != None:
                conflicts = False
                for pair in COOLPROP_INCOMPATIBLE_PAIRS:
                    if cpName == pair[0] and pair[1] in mix:
                        conflicts = True
                    if cpName == pair[1] and pair[0] in mix:
                        conflicts = True
                
                if not conflicts:
                    filtered.data.append(part)
                    mix += cpName + "[" + str(part[1]) + "]" + "&"

        filtered.normalize()

        mix = "HEOS::"
        for part in filtered.data:
            cpName = cea_to_cp_name(part[0])
            mix += cpName + "[" + str(part[1]) + "]" + "&"
            
        return mix[:-1]

    def from_cea(self, ceaData):
        self.data = []
        for species, fraction in ceaData.mole_fractions.items():
            if fraction > 0.00001:
                self.data.append([species, fraction])
        self.normalize()

    def normalize(self):
        n_tot = 0
        for part in self.data:
            n_tot += part[1]
        for i in range(0, len(self.data)):
            self.data[i][1] = self.data[i][1] / n_tot
        if n_tot < 0.95:
            print("WARNING:", str(round(100 * (1-n_tot), 1)), "% of moles of mixture were filtered out")

    def sort(self):
        for i in range(0, len(self.data)):
            self.data[i] = [self.data[i][1], self.data[i]]
        self.data = sorted(self.data, reverse=True)
        for i in range(0, len(self.data)):
            self.data[i] = self.data[i][1]

    def pie_chart(self, title, temp=None, ax=None):
        title += " Mole Fractions"
        fracs = []
        labels = []
        for part in self.data:
            fracs.append(part[1])
            if len(part) > 2:
                labels.append(part[0] + "  T = " + str(round(part[2], 2)) + " K  (" + str(round(100*part[1], 3)) + "%)")
            elif temp != None:
                labels.append(part[0] + "  T = " + str(round(temp, 2)) + " K  (" + str(round(100*part[1], 3)) + "%)")
            else:
                labels.append(part[0] + "  (" + str(round(100*part[1], 3)) + "%)")

        if ax == None:
            ax = plt.figure().subplots()
        wedges = ax.pie(fracs)[0]
        ax.set_title(title)
        ax.legend(wedges, labels, loc="upper center", bbox_to_anchor=(0.5, 0))
        plt.tight_layout()
        plt.show(block=False)

def AIR():
    return Mixture([["N2", 0.78084], ["O2", 0.20946], ["Ar", 0.00934], ["CO2", 0.000412], ["Ne", 0.00001818]])
    # from https://www.engineeringtoolbox.com/air-composition-d_212.html

class SetOfStates():
    def __init__(self, numStates, fluid = "Water"):
        self.T = [None] * numStates
        self.P = [None] * numStates
        self.v = [None] * numStates
        self.h = [None] * numStates
        self.s = [None] * numStates
        self.X = [None] * numStates
        self.fluid = [fluid] * numStates
        self.consistent = False # True if "Air" h and s were converted to match other mixtures

    def make_consistent(self, newConsistency):
        if self.consistent == newConsistency:
            return

        if newConsistency == True:
            for i in range(0, len(self.T)):
                if self.fluid[i].data == AIR().data:
                    self.h[i] -= 125886
                    self.s[i] += 2978

        else:
            for i in range(0, len(self.T)):
                if self.fluid[i].data == AIR().data:
                    self.h[i] += 125886
                    self.s[i] -= 2978
        
        self.consistent = newConsistency

    def extend(self, otherState):
        selfConsistent = self.consistent
        otherConsistent = otherState.consistent
        self.make_consistent(True)
        otherState.make_consistent(True)
        self.T.extend(otherState.T)
        self.P.extend(otherState.P)
        self.v.extend(otherState.v)
        self.h.extend(otherState.h)
        self.s.extend(otherState.s)
        self.X.extend(otherState.X)
        self.fluid.extend(otherState.fluid)
        self.make_consistent(selfConsistent)
        otherState.make_consistent(otherConsistent)

    def get_point(self, i):
        point = SetOfStates(1)
        point.T[0] = self.T[i]
        point.P[0] = self.P[i]
        point.v[0] = self.v[i]
        point.h[0] = self.h[i]
        point.s[0] = self.s[i]
        point.X[0] = self.X[i]
        point.fluid[0] = copy.deepcopy(self.fluid[i])
        return point

    def solve_points(self):
        self.make_consistent(False)
        
        for i in range(0, len(self.T)):
            knowns = []
            if self.T[i] != None:
                knowns.append(["T", self.T[i]])
            if self.P[i] != None:
                knowns.append(["P", self.P[i]])
            if self.v[i] != None:
                knowns.append(["D", 1 / self.v[i]])
            if self.h[i] != None:
                knowns.append(["Hmass", self.h[i]])
            if self.s[i] != None:
                knowns.append(["Smass", self.s[i]])
            if self.X[i] != None:
                knowns.append(["Q", self.X[i]])

            if len(knowns) < 2:
                #print("Couldn't solve point", i)
                continue

            if self.T[i] == None:
                self.T[i] = cp.PropsSI("T", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i].to_coolprop())
            if self.P[i] == None:
                self.P[i] = cp.PropsSI("P", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i].to_coolprop())
            if self.v[i] == None:
                self.v[i] = 1 / cp.PropsSI("D", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i].to_coolprop())
            if self.h[i] == None:
                self.h[i] = cp.PropsSI("Hmass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i].to_coolprop())
            if self.s[i] == None:
                self.s[i] = cp.PropsSI("Smass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i].to_coolprop())
            if self.X[i] == None:
                self.X[i] = cp.PropsSI("Q", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i].to_coolprop())

    
    def point_is_solved(self, i):
        if self.T[i] == None: return False
        if self.P[i] == None: return False
        if self.h[i] == None: return False
        if self.s[i] == None: return False
        if self.X[i] == None: return False
        return True

    def set_is_fully_solved(self):
        for i in range(0, len(self.T)):
            if not self.point_is_solved(i):
                return False
        return True

    def polytropic_interpolation(self, a, b, numPoints = 3):
        # doesn't work well and currently not in use
        # leaving here in case it becomes useflul later and I return to it

        
        intermediates = SetOfStates(0)
        #if self.X[a] != self.X[b]:
        #    return intermediates
        
        k = math.log(self.P[a] / self.P[b]) / math.log(self.v[b] / self.v[a])
        c = self.P[a] + self.v[a] ** k

        print()
        print(self.P[a], self.P[b])
        print(self.v[a], self.v[b])
        print(k, c)

        for i in range(0, numPoints):
            v = 0
            P = 0

            if k < 1:
                v = self.v[a] + (self.v[b] - self.v[a]) * i / numPoints
                P = c / v**k
            else:
                P = self.P[a] + (self.P[b] - self.P[a]) * i / numPoints
                v = (P / c) ** (1/k)

            print(P, "\t", v)
            
            intermediates.T.append(None)
            intermediates.P.append(P)
            intermediates.v.append(v)
            intermediates.h.append(None)
            intermediates.s.append(None)
            intermediates.X.append(None)
            intermediates.fluid.append(self.fluid[a])
        try:
            intermediates.solve_points()
        except:
            return SetOfStates(0)
        return intermediates

    def get_property(self, i, symbol):
        if symbol == "T":
            return self.T[i]
        if symbol == "P":
            return self.P[i]
        if symbol == "v":
            return self.v[i]
        if symbol == "h":
            return self.h[i]
        if symbol == "s":
            return self.s[i]
        if symbol == "X":
            return self.X[i]





def intermediate_points_worker_pool(linearProperties, starts, ends, fluid1, fluid2, i, numPoints):
    linValues = []
    for j in range(0, len(linearProperties)):
        val = starts[j] + (ends[j] - starts[j]) * i / numPoints
        if linearProperties[j] == "D":
            linValues.append(1 / val)
        else:
            linValues.append(val)

    data = []
    data.append(cp.PropsSI("T", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid1))
    data.append(cp.PropsSI("P", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid1))
    data.append(1 / cp.PropsSI("D", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid1))
    data.append(cp.PropsSI("Hmass", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid1))
    data.append(cp.PropsSI("Smass", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid1))
    data.append(cp.PropsSI("Q", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid1))
    
    if fluid1 == "Air":
            data[3] -= 125866
            data[4] += 2978

    if fluid1 != fluid2: # interpolating between fluid properties if there is a chemical change during a process
                                   # noticable discrete jumps are visible on Ts and Pv diagrams if this is not done
        data2 = []
        data2.append(cp.PropsSI("T", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid2))
        data2.append(cp.PropsSI("P", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid2))
        data2.append(1 / cp.PropsSI("D", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid2))
        data2.append(cp.PropsSI("Hmass", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid2))
        data2.append(cp.PropsSI("Smass", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid2))
        data2.append(cp.PropsSI("Q", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fluid2))

        if fluid2 == "Air":
            data2[3] -= 125866
            data2[4] += 2978
        

        for k in range(0, len(data)):
            data[k] = data[k] + (data2[k] - data[k]) * i / numPoints
    
    return data


def get_intermediate_points(fromSet, fromI, toSet, toI, linearProperties, numPoints=30):
    # for linearProperties use [most changed, constant when ideal]

    fromConsistent = fromSet.consistent
    toConsistent = toSet.consistent
    fromSet.make_consistent(False)
    toSet.make_consistent(False)
    
    intermediates = SetOfStates(0)
    starts = [fromSet.get_property(fromI, symbol) for symbol in linearProperties]
    ends = [toSet.get_property(toI, symbol) for symbol in linearProperties]

    for i in range(0, len(linearProperties)):
        if linearProperties[i] == "h":
            linearProperties[i] = "Hmass"
        if linearProperties[i] == "s":
            linearProperties[i] = "Smass"
        if linearProperties[i] == "X":
            linearProperties[i] = "Q"
        if linearProperties[i] == "v":
            linearProperties[i] = "D"

    fluid1 = fromSet.fluid[fromI].to_coolprop()
    fluid2 = toSet.fluid[toI].to_coolprop()

    tasks = [[linearProperties, starts, ends, fluid1, fluid2, i, numPoints] for i in range(1, numPoints)]

    
    with multiprocessing.Pool(processes=3) as pool:
        fullData = pool.starmap(intermediate_points_worker_pool, tasks)

    intermediates.extend(fromSet.get_point(fromI))

    for point in fullData:
        intermediates.T.append(point[0])
        intermediates.P.append(point[1])
        intermediates.v.append(point[2])
        intermediates.h.append(point[3])
        intermediates.s.append(point[4])
        intermediates.X.append(point[5])
        intermediates.fluid.append(fromSet.fluid[fromI])
        
    intermediates.extend(toSet.get_point(toI))
    intermediates.consistent = True

    if fluid1 == "Air":
        intermediates.h[0] -= 125866
        intermediates.s[0] += 2978

    if fluid2 == "Air":
        intermediates.h[-1] -= 125866
        intermediates.s[-1] += 2978

    
    fromSet.make_consistent(fromConsistent)
    toSet.make_consistent(toConsistent)
    return intermediates




def solve_point_inefficient(toSet, i, fromSet, fromIndex, enthalpyMultiplier):
    knowns = []
    if toSet.T[i] != None:
        knowns.append(["T", toSet.T[i]])
    if toSet.P[i] != None:
        knowns.append(["P", toSet.P[i]])
    if toSet.v[i] != None:
        knowns.append(["D", 1 / self.v[i]])
    if toSet.X[i] != None:
        knowns.append(["Q", toSet.X[i]])


    if len(knowns) == 0:
        print("Couldn't solve point", i, "in solve_point_inefficient")
        return

    h_isen = cp.PropsSI("Hmass", knowns[0][0], knowns[0][1], "S", fromSet.s[fromIndex], toSet.fluid[i].to_coolprop())
    toSet.h[i] = fromSet.h[fromIndex] + enthalpyMultiplier * (h_isen - fromSet.h[fromIndex])
    
    knowns.append(["Hmass", toSet.h[i]])

    if toSet.T[i] == None:
        toSet.T[i] = cp.PropsSI("T", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i].to_coolprop())
    if toSet.P[i] == None:
        toSet.P[i] = cp.PropsSI("P", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i].to_coolprop())
    if toSet.v[i] == None:
        toSet.v[i] = 1 / cp.PropsSI("D", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i].to_coolprop())
    if toSet.s[i] == None:
        toSet.s[i] = cp.PropsSI("Smass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i].to_coolprop())
    if toSet.X[i] == None:
        toSet.X[i] = cp.PropsSI("Q", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i].to_coolprop())










class Cycle: # should be valid for both rankine and brayton
    def __init__(self, numComp, numTurb, fluid):
        self.preCompStates = SetOfStates(numComp, copy.deepcopy(fluid))
        self.postCompStates = SetOfStates(numComp, copy.deepcopy(fluid))
        self.heatingStates = SetOfStates(2, copy.deepcopy(fluid)) # heatingStates = [before fuel injected, before combustion but after fuel injected]
        self.preTurbStates = SetOfStates(numTurb, copy.deepcopy(fluid))
        self.postTurbStates = SetOfStates(numTurb, copy.deepcopy(fluid))
        self.turbineEfficiencies = [1] * numTurb
        self.compressorEfficiencies = [1] * numComp
        
        self.externalFuelCompressionWork = 0 # work to pressurize fuel before injecting it into combustion chamber
        self.addedFuelMassFrac = 0
        self.addefFuelInitialEnthalpy = 0
        self.recycledHeat = 0

        self.allStateSets = [self.preCompStates, self.postCompStates, self.preTurbStates, self.postTurbStates]
        self.intermediates = None
        

    def combined_states(self):
        self.preCompStates.make_consistent(False)
        self.postCompStates.make_consistent(False)
        self.heatingStates.make_consistent(False)
        self.preTurbStates.make_consistent(False)
        self.postTurbStates.make_consistent(False)
        states = SetOfStates(0)
        for i in range(0, len(self.preCompStates.T)):
            states.extend(self.preCompStates.get_point(i))
            states.extend(self.postCompStates.get_point(i))

        if self.heatingStates.set_is_fully_solved():
            states.extend(self.heatingStates)
            
        for i in range(0, len(self.preTurbStates.T)):
            states.extend(self.preTurbStates.get_point(i))
            states.extend(self.postTurbStates.get_point(i))
        return states # returned object is READ ONLY


    def solve_cp(self):
        for subset in self.allStateSets:
            subset.solve_points()

        for i in range(0, len(self.preCompStates.T)):
            if self.preCompStates.point_is_solved(i) and self.postCompStates.point_is_solved(i):
                h_isen = cp.PropsSI("Hmass", "P", self.postCompStates.P[i], "Smass", self.preCompStates.s[i], self.preCompStates.fluid[i].to_coolprop())
                self.compressorEfficiencies[i] = (h_isen - self.preCompStates.h[i]) / (self.postCompStates.h[i] - self.preCompStates.h[i])
                #print("Compressor efficiency overwritten to", self.compressorEfficiencies[i])
            elif self.preCompStates.point_is_solved(i):
                solve_point_inefficient(self.postCompStates, i, self.preCompStates, i, 1 / self.compressorEfficiencies[i])
            elif self.postCompStates.point_is_solved(i):
                solve_point_inefficient(self.preCompStates, i, self.postCompStates, i, self.compressorEfficiencies[i])

        for i in range(0, len(self.preTurbStates.T)):
            if self.preTurbStates.point_is_solved(i) and self.postTurbStates.point_is_solved(i):
                h_isen = cp.PropsSI("Hmass", "P", self.postTurbStates.P[i], "Smass", self.preTurbStates.s[i], self.preTurbStates.fluid[i].to_coolprop())
                self.turbineEfficiencies[i] = (self.postTurbStates.h[i] - self.preTurbStates.h[i]) / (h_isen - self.preTurbStates.h[i])
                #print("Turbine efficiency overwritten to", self.turbineEfficiencies[i])
            elif self.preTurbStates.point_is_solved(i):
                solve_point_inefficient(self.postTurbStates, i, self.preTurbStates, i, self.turbineEfficiencies[i])
            elif self.postTurbStates.point_is_solved(i):
                solve_point_inefficient(self.preTurbStates, i, self.postTurbStates, i, 1 / self.turbineEfficiencies[i])

    def calculate_pre_combustion_heat_exchange(self, q_in):
        self.solve_cp()
        self.heatingStates.P[0] = self.postCompStates.P[-1]
        self.heatingStates.h[0] = self.postCompStates.h[-1] + q_in
        self.heatingStates.fluid[0] = self.postCompStates.fluid[-1]
        self.heatingStates.P[1] = self.heatingStates.P[0]
        self.heatingStates.h[1] = self.heatingStates.h[0]
        self.heatingStates.fluid[1] = self.heatingStates.fluid[0]
        self.heatingStates.solve_points()
        

    def inject_fuel_after_compressor_no_propagation(self, fuelMixture, fuelTotalMassFraction, fuelCompressorEfficiency):
        self.addedFuelMassFrac = fuelTotalMassFraction
        if not self.heatingStates.set_is_fully_solved():
            self.calculate_pre_combustion_heat_exchange(0)
            print("Assumed no pre-combustion heat exchange when injecting fuel")

        # compressing fuel from standard atmospheric conditions to same pressure as air at exit of compressor
        # then finding work done and final fuel temperature during injection
        h0 = cp.PropsSI("Hmass", "P", 101325, "T", 298, fuelMixture.to_coolprop())
        s0 = cp.PropsSI("Smass", "P", 101325, "T", 298, fuelMixture.to_coolprop())
        h1 = cp.PropsSI("Hmass", "P", self.heatingStates.P[0], "S", s0, fuelMixture.to_coolprop())
        h1 = h0 + (h1 - h0) / fuelCompressorEfficiency
        fuelTemp = cp.PropsSI("T", "P", self.heatingStates.P[0], "Hmass", h1, fuelMixture.to_coolprop())
        self.externalFuelCompressionWork = (h1 - h0) * fuelTotalMassFraction
        self.addefFuelInitialEnthalpy = h0

        M_oxidiser = cp.PropsSI("molar_mass", "P", self.heatingStates.P[0], "T", self.heatingStates.T[0], self.heatingStates.fluid[0].to_coolprop())
        M_fuel = cp.PropsSI("molar_mass", "P", self.heatingStates.P[0], "T", fuelTemp, fuelMixture.to_coolprop())
        n_oxidiser = (1 - fuelTotalMassFraction) / M_oxidiser
        n_fuel = fuelTotalMassFraction / M_fuel
        fuelMoleFraction = n_fuel / (n_oxidiser + n_fuel)

        fluid = copy.deepcopy(self.heatingStates.fluid[0])
        for i in range(0, len(fluid.data)):
            fluid.data[i][1] = fluid.data[i][1] * (1 - fuelMoleFraction)
        for i in range(0, len(fuelMixture.data)):
            fluid.data.append([fuelMixture.data[i][0], fuelMixture.data[i][1] * fuelMoleFraction, fuelTemp])

        self.heatingStates.fluid[1] = fluid


    def propagate_gas_from_combustor(self):
        for i in range(0, len(self.preTurbStates.T)):
            self.preTurbStates.fluid[i] = copy.deepcopy(self.preTurbStates.fluid[0])
            self.postTurbStates.fluid[i] = copy.deepcopy(self.preTurbStates.fluid[0])


    def combust(self):
        reactants, massFracs, temps = self.heatingStates.fluid[1].to_cea(self.heatingStates.T[1])
        products = self.heatingStates.fluid[1].to_cea(self.heatingStates.T[1], products_from_reactants=True)[0]
        solver = cea.EqSolver(products, reactants=reactants, transport=True)
        solution = cea.EqSolution(solver)
        #print(temps)
        massFracs = np.array(massFracs)
        hc = reactants.calc_property(cea.ENTHALPY, massFracs, np.array(temps)) / cea.R
        solver.solve(solution, cea.HP, hc, self.heatingStates.P[1] / 100000, massFracs)
        self.preTurbStates.T[0] = solution.T
        products = Mixture()
        products.from_cea(solution)
        self.preTurbStates.fluid[0] = products
        self.heatingStates.s[1] = None
        self.heatingStates.h[1] = None
        self.heatingStates.v[1] = None
        self.heatingStates.X[1] = None
        self.heatingStates.solve_points()
        
        self.propagate_gas_from_combustor()
        self.solve_cp()
        
        

    def compressor_works(self):
        self.preCompStates.make_consistent(True)
        self.postCompStates.make_consistent(True)
        w = []
        for i in range(0, len(self.preCompStates.T)):
            w.append((self.postCompStates.h[i] - self.preCompStates.h[i]) * (1 - self.addedFuelMassFrac))
        return w

    def intercooler_heats(self):
        self.preCompStates.make_consistent(True)
        self.postCompStates.make_consistent(True)
        q = []
        for i in range(1, len(self.preCompStates.T)):
            q.append((self.postCompStates.h[i-1] - self.preCompStates.h[i]) * (1 - self.addedFuelMassFrac))
        return q

    def turbine_works(self):
        self.preTurbStates.make_consistent(True)
        self.postTurbStates.make_consistent(True)
        w = []
        for i in range(0, len(self.preTurbStates.T)):
            w.append(self.preTurbStates.h[i] - self.postTurbStates.h[i])
        return w

    def reheater_heats(self):
        self.preTurbStates.make_consistent(True)
        self.postTurbStates.make_consistent(True)
        q = []
        for i in range(1, len(self.preTurbStates.T)):
            q.append(self.preTurbStates.h[i] - self.postTurbStates.h[i-1])
        return q

    def return_heat(self):
        self.postTurbStates.make_consistent(True)
        self.preCompStates.make_consistent(True)
        fuel = (self.postTurbStates.h[-1] - self.addefFuelInitialEnthalpy) * self.addedFuelMassFrac
        air = (self.postTurbStates.h[-1] - self.preCompStates.h[0]) * (1 - self.addedFuelMassFrac)
        return fuel + air

    def heat_exchanger_heat(self):
        self.postCompStates.make_consistent(True)
        if self.addedFuelMassFrac == 0:
            self.preTurbStates.make_consistent(True)
            return self.preTurbStates.h[0] - self.postCompStates.h[-1]
        else:
            self.heatingStates.make_consistent(True)
            return (self.heatingStates.h[0] - self.postCompStates.h[-1]) * (1 - self.addedFuelMassFrac)

    def combustion_heat(self):
        if self.addedFuelMassFrac == 0:
            return 0
        else:
            self.heatingStates.make_consistent(True)
            self.preTurbStates.make_consistent(True)
            fuel = (self.preTurbStates.h[0] - self.addefFuelInitialEnthalpy - self.externalFuelCompressionWork) * self.addedFuelMassFrac
            air = (self.preTurbStates.h[0] - self.heatingStates.h[0]) * (1 - self.addedFuelMassFrac)
            return fuel + air

    def q_in(self): # assuming heat exchanger source is external, not recycling heat
        return self.heat_exchanger_heat() + self.combustion_heat() + sum(self.reheater_heats())

    def q_out(self): # assuming no heat is recycled using a heat exchanger
        return sum(self.intercooler_heats()) + self.return_heat()

    def w_in(self):
        return sum(self.compressor_works()) + self.externalFuelCompressionWork

    def w_out(self):
        return sum(self.turbine_works())

    def w_net(self):
        return self.w_out() - self.w_in()

    def efficiency(self):
        return self.w_net() / (self.w_in() + self.q_in() - self.recycledHeat)

    def back_work(self):
        return self.w_in() / self.w_out()

    def print_energy_calcs(self):
        print()
        print("---------------------------------------------------------------")
        print()
        print("Energy Calculations:")
        print("All units in J/kg noramalized using total mass of air+fuel")
        print()
        print("Heat exchanger heat:", self.heat_exchanger_heat())
        print("Combustion heat:", self.combustion_heat())
        reheats = self.reheater_heats()
        for i in range(0, len(reheats)):
            print("Reheater", i+1, "heat:", reheats[i])
        print("Total q_in:", self.q_in())
        print("Recycled heat:", self.recycledHeat)
        print("Total external q_in:", self.q_in() - self.recycledHeat)
        print()
        intercools = self.intercooler_heats()
        for i in range(0, len(intercools)):
            print("Intercooler", i+1, "heat:", intercools[i])
        print("Return heat:", self.return_heat())
        print("Total q_out:", self.q_out())
        print()
        comps = self.compressor_works()
        for i in range(0, len(comps)):
            print("Compressor", i+1, "( η =", round(self.compressorEfficiencies[i], 4), ") work:", comps[i])
        print("Fuel compression work:", self.externalFuelCompressionWork)
        print("Total w_in :", self.w_in())
        print()
        turbs = self.turbine_works()
        for i in range(0, len(turbs)):
            print("Turbine", i+1, "( η =", round(self.turbineEfficiencies[i], 4), ") work:", turbs[i])
        print("Total w_out:", self.w_out())
        print()
        print("Net work:", self.w_net())
        print("Efficiency:", self.efficiency())
        print("Back-work ratio:", self.back_work())
        if self.recycledHeat > self.heat_exchanger_heat():
            print()
            print("WARNING: recycled heat is greater than heat entering through heat exchanger")
        print()
        print("---------------------------------------------------------------")
        print()



    def update_intermediates(self):
        self.intermediates = SetOfStates(0)

        #print("test 1")
        
        for i in range(0, len(self.preCompStates.T)): # compressor
            self.intermediates.extend(get_intermediate_points(self.preCompStates, i, self.postCompStates, i, ["h", "s"]))

            if i+1 < len(self.preCompStates.T): # intercooling
                self.intermediates.extend(get_intermediate_points(self.postCompStates, i, self.preCompStates, i+1, ["h", "P"]))

        #print("test 2")
        # boiler / combustion chamber
        
        if self.heatingStates.h[0] != None and self.heatingStates.h[0] != self.postTurbStates.h[-1]: # if heat exchange before combustion
            # heat exchanger process:
            self.intermediates.extend(get_intermediate_points(self.postCompStates, -1, self.heatingStates, 0, ["s", "P"]))
            # injection process:
            self.intermediates.extend(get_intermediate_points(self.heatingStates, 0, self.heatingStates, 1, ["T", "P"], numPoints=2))
            # combustion process:
            self.intermediates.extend(get_intermediate_points(self.heatingStates, 1, self.preTurbStates, 0, ["T", "P"]))
        
        else: # single process combustion (no heat exchanges)
            try:
                self.intermediates.extend(get_intermediate_points(self.postCompStates, -1, self.preTurbStates, 0, ["s", "P"]))
                #print("success 1c")
            except:
                self.intermediates.extend(get_intermediate_points(self.postCompStates, -1, self.preTurbStates, 0, ["T", "P"]))
                #print("fail 1c")
        
        #print("test 3")
        
        for i in range(0, len(self.preTurbStates.T)): # turbine
            self.intermediates.extend(get_intermediate_points(self.preTurbStates, i, self.postTurbStates, i, ["h", "s"]))

            if i+1 < len(self.preTurbStates.T): # reheating
                self.intermediates.extend(get_intermediate_points(self.postTurbStates, i, self.preTurbStates, i+1, ["T", "P"]))


        #print("test 4")
        
        # condenser
        separateChemistryStep = False
        dontDrawComplexBraytonCooling = True

        
        if self.postTurbStates.fluid[-1].data != self.preCompStates.fluid[0].data and separateChemistryStep and not dontDrawComplexBraytonCooling:
            intermediateStep = SetOfStates(1)
            intermediateStep.P[0] = self.preCompStates.P[0]
            intermediateStep.T[0] = self.preCompStates.T[0]
            intermediateStep.fluid[0] = self.postTurbStates.fluid[-1]
            intermediateStep.solve_points()

            # cooling:
            self.intermediates.extend(get_intermediate_points(self.postTurbStates, -1, intermediateStep, 0, ["T", "P"]))
            # undoing chemical mixing
            self.intermediates.extend(get_intermediate_points(intermediateStep, 0, self.preCompStates, 0, ["T", "P"], numPoints=2))

        elif self.postTurbStates.fluid[-1].data == self.preCompStates.fluid[0].data or not dontDrawComplexBraytonCooling: # single step
            try:
                self.intermediates.extend(get_intermediate_points(self.postTurbStates, -1, self.preCompStates, 0, ["s", "P"]))
                #print("success 2c")
            except:
                self.intermediates.extend(get_intermediate_points(self.postTurbStates, -1, self.preCompStates, 0, ["T", "P"]))
                #print("fail 2c")

        self.intermediates.make_consistent(True)

    
    def ts_diagram(self, ax=None, label=None, isRankine=False):
        allPoints = self.combined_states()
        allPoints.make_consistent(True)
        
        if ax == None:
            ax = plt.figure().add_subplot()
            ax.set_xlabel("s - $J kg^{-1} K^{-1}$")
            ax.set_ylabel("T - $K$")

        x = []
        y = []
        for i in range(0, len(allPoints.T)):
            x.append(allPoints.s[i])
            y.append(allPoints.T[i])

        ax.scatter(x, y, label=label)

        
        
        if self.intermediates == None:
            self.update_intermediates()

        x = self.intermediates.s
        y = self.intermediates.T
        
        ax.plot(x, y)
        #ax.scatter(x, y, marker=".") # for debug to ensure no large jumps

        if isRankine:
            original_xlim = ax.get_xlim()
            original_ylim = ax.get_ylim()
            propertyPlot = PropertyPlot(self.postCompStates.fluid[0].to_coolprop(), 'TS', axis=ax, unit_system='SI')
            propertyPlot.calc_isolines(CoolProp.iQ, num=2)
            plt.close(propertyPlot.figure)
            propertyPlot.draw()
            ax.set_xlim(original_xlim)
            ax.set_ylim(original_ylim)
            ax.autoscale(True)
        
        if label != None:
            ax.legend()
        plt.show(block=False)
        return ax
    
    def pv_diagram(self, ax=None, label=None):
        allPoints = self.combined_states()
        allPoints.make_consistent(True)
        
        if ax == None:
            ax = plt.figure().add_subplot()
            ax.set_xlabel("v - $m^3 kg^{-1}$")
            ax.set_ylabel("P - $Pa$")

        x = []
        y = []
        for i in range(0, len(allPoints.T)):
            x.append(allPoints.v[i])
            y.append(allPoints.P[i])

        ax.scatter(x, y, label=label)

        
        
        if self.intermediates == None:
            self.update_intermediates()

        x = self.intermediates.v
        y = self.intermediates.P
        
        ax.plot(x, y)
        #ax.scatter(x, y, marker=".") # for debug to ensure no large jumps
        
        if label != None:
            ax.legend()
        plt.show(block=False)
        return ax


    def mixture_pie_charts(self, title=None):
        if self.addedFuelMassFrac != 0:
            fig = plt.figure()
            fig.suptitle(title)
            axs = fig.subplots(1, 3)
            self.preCompStates.fluid[0].pie_chart("Inlet", ax=axs[0])
            self.heatingStates.fluid[1].pie_chart("Pre-Combustion", temp=self.heatingStates.T[1], ax=axs[1])
            self.preTurbStates.fluid[0].pie_chart("Post-Combustion", temp=self.heatingStates.T[0], ax=axs[2])

        else:
            self.preCompStates.fluid[0].pie_chart(title, ax=axs[0])





def graph_efficiency_vs_fuel_mass_frac(parentCycle, fuel, fuelCompressorEfficiency, minFrac, maxFrac, numPoints):
    fracs = []
    efficiencies = []
    for i in range(0, numPoints):
        frac = minFrac + (maxFrac - minFrac) * i / numPoints
        cycle = copy.deepcopy(parentCycle)
        cycle.inject_fuel_after_compressor_no_propagation(Mixture([["CH4", 1]]), frac, fuelCompressorEfficiency)
        cycle.combust()
        efficiencies.append(cycle.efficiency())
        fracs.append(frac)
    
    plt.figure()
    plt.xlabel("Fuel Mass fraction")
    plt.ylabel("Cycle Efficiency")
    plt.scatter(fracs, efficiencies)
    plt.show(block=False)




if __name__ == "__main__":
    multiprocessing.set_start_method('spawn')
    

    testRankine = Cycle(1, 1, Mixture([["H2O", 1]]))
    testRankine.preCompStates.P = [8000]
    testRankine.preCompStates.X = [0]
    testRankine.postCompStates.P = [8000000]
    testRankine.preTurbStates.T = [600]
    testRankine.preTurbStates.P = [8000000]
    testRankine.postTurbStates.P = [8000]
    testRankine.turbineEfficiencies[0] = 0.85
    testRankine.compressorEfficiencies[0] = 0.8

    testRankine.solve_cp()
    rankineTS = testRankine.ts_diagram(label="unmodified", isRankine=True)
    rankineTS.set_title("Rankine Test Plot Ts")
    
    print("Unmodified Rankine")
    testRankine.print_energy_calcs()

    print("Rankine with 1 reheater")

    testRankine = Cycle(1, 2, Mixture([["H2O", 1]]))
    testRankine.preCompStates.P = [8000]
    testRankine.preCompStates.X = [0]
    testRankine.postCompStates.P = [8000000]
    testRankine.preTurbStates.T = [600, 600]
    testRankine.preTurbStates.P = [8000000, 4000000]
    testRankine.postTurbStates.P = [4000000, 8000]
    testRankine.turbineEfficiencies = [0.85, 0.85]
    testRankine.compressorEfficiencies = [0.8]

    testRankine.solve_cp()
    testRankine.print_energy_calcs()
    testRankine.ts_diagram(ax=rankineTS, label = "reheater")

    print("Brayton")

    testBrayton = Cycle(1, 1, AIR())
    testBrayton.preCompStates.P = [101325]
    testBrayton.postCompStates.P = [500000]
    testBrayton.preCompStates.T = [298]
    #testBrayton.preTurbStates.T = [600]
    testBrayton.preTurbStates.P = [500000]
    testBrayton.postTurbStates.P = [101325]

    testBrayton.calculate_pre_combustion_heat_exchange(100000)
    graph_efficiency_vs_fuel_mass_frac(testBrayton, Mixture([["CH4", 1]]), 0.9, 0.04, 0.12, 30)
    #import sys
    #sys.exit()

    testBrayton.inject_fuel_after_compressor_no_propagation(Mixture([["CH4", 1]]), 0.08, 0.9)
    testBrayton.combust()

    #testBrayton.solve_cp()
    testBrayton.print_energy_calcs()
    
    braytonTS = testBrayton.ts_diagram(label="isentropic")
    braytonPV = testBrayton.pv_diagram(label="isentropic")
    #print(testBrayton.efficiency())

    braytonTS.set_title("Brayton Test Plot Ts")
    braytonPV.set_title("Brayton Test Plot Pv")

    testBrayton.mixture_pie_charts("Test Brayton")


    """
    testBrayton = Cycle(1, 1, Mixture([["AIR", 1]]))
    testBrayton.preCompStates.P = [101325]
    testBrayton.postCompStates.P = [500000]
    testBrayton.preCompStates.T = [298]
    testBrayton.preTurbStates.T = [600]
    testBrayton.preTurbStates.P = [500000]
    testBrayton.postTurbStates.P = [101325]
    testBrayton.turbineEfficiencies = [0.7]
    testBrayton.compressorEfficiencies = [0.7]

    testBrayton.solve_cp()
    print(testBrayton.efficiency())
    testBrayton.ts_diagram(ax=braytonTS, label="70% efficient turbine&compressor")
    testBrayton.pv_diagram(ax=braytonPV, label="70% efficient turbine&compressor")


    testBrayton = Cycle(2, 2, Mixture([["AIR", 1]]))
    testBrayton.preCompStates.P = [101325, 250000]
    testBrayton.postCompStates.P = [250000, 500000]
    testBrayton.preCompStates.T = [298, 298]
    testBrayton.preTurbStates.T = [600, 600]
    testBrayton.preTurbStates.P = [500000, 250000]
    testBrayton.postTurbStates.P = [250000, 101325]
    testBrayton.turbineEfficiencies = [0.7, 0.7]
    testBrayton.compressorEfficiencies = [0.7, 0.7]

    testBrayton.solve_cp()
    print(testBrayton.efficiency())
    testBrayton.ts_diagram(ax=braytonTS, label="70% efficient turbine&compressor with intercooler and reheating")
    testBrayton.pv_diagram(ax=braytonPV, label="70% efficient turbine&compressor with intercooler and reheating")
    """




