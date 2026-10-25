import cea
from CoolProp import CoolProp as cp # haha very funny
from CoolProp.Plots import PropertyPlot
import CoolProp
import matplotlib.pyplot as plt
import math
import copy

CEA_TO_COOLPROP_NAMES = {
    "CH4" : "Methane",
    "O2" : "Oxygen",
    "H2" : "Hydrogen",
    "N2" : "Nitrogen",
    "CO2" : "CarbonDioxide",
    "CO" : "CarbonMonoxide",
    "H2O" : "Water",
    "AIR" : "Air",
    "Propane" : "Propane",
    "Propylene" : "Propylene",
    "Ethanol" : "Ethanol",
}

class Mixture:
    def __init__(self, data=None):
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
            M = cp.PropsSI("molar_mass", self.to_coolprop())
            m = molFracs[i] * M
            massFracs.append(m)
            totMass += m
        for i in range(0, len(names)):
            massFracs[i] = massFracs[i] / totMass
        
        return [mix, massFracs, temps]

    def to_coolprop(self):
        mix = "HEOS::"
        for part in self.data:
            mix += CEA_TO_COOLPROP_NAMES[part[0]] + "[" + str(part[1]) + "]" + "&"
        return mix[:-1]

    def from_cea(self, ceaData):
        self.data = []
        for species, fraction in ceaData.mole_fractions.items():
            self.data.append([species, fraction[2]])

class SetOfStates():
    def __init__(self, numStates, fluid = "Water"):
        self.T = [None] * numStates
        self.P = [None] * numStates
        self.v = [None] * numStates
        self.h = [None] * numStates
        self.s = [None] * numStates
        self.X = [None] * numStates
        self.fluid = [fluid] * numStates

    def extend(self, otherState):
        self.T.extend(otherState.T)
        self.P.extend(otherState.P)
        self.v.extend(otherState.v)
        self.h.extend(otherState.h)
        self.s.extend(otherState.s)
        self.X.extend(otherState.X)
        self.fluid.extend(otherState.fluid)

    def get_point(self, i):
        point = SetOfStates(1)
        point.T[0] = self.T[i]
        point.P[0] = self.P[i]
        point.v[0] = self.v[i]
        point.h[0] = self.h[i]
        point.s[0] = self.s[i]
        point.X[0] = self.X[i]
        point.fluid[0] = self.fluid[i]
        return point

    def solve_points(self):
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





def get_intermediate_points(fromSet, fromI, toSet, toI, linearProperties, numPoints=100):
    # for linearProperties use [most changed, constant when ideal]
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
    
    for i in range(0, numPoints):
        linValues = []
        for j in range(0, len(linearProperties)):
            val = starts[j] + (ends[j] - starts[j]) * i / numPoints
            if linearProperties[j] == "D":
                linValues.append(1 / val)
            else:
                linValues.append(val)
        
        intermediates.T.append(cp.PropsSI("T", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI].to_coolprop()))
        intermediates.P.append(cp.PropsSI("P", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI].to_coolprop()))
        intermediates.v.append(1 / cp.PropsSI("D", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI].to_coolprop()))
        intermediates.h.append(cp.PropsSI("Hmass", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI].to_coolprop()))
        intermediates.s.append(cp.PropsSI("Smass", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI].to_coolprop()))
        intermediates.X.append(cp.PropsSI("Q", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI].to_coolprop()))
        intermediates.fluid.append(fromSet.fluid[fromI])
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

        self.allStateSets = [self.preCompStates, self.postCompStates, self.preTurbStates, self.postTurbStates]
        self.intermediates = None
        

    def combined_states(self):
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
                print("Compressor efficiency overwritten to", self.compressorEfficiencies[i])
            elif self.preCompStates.point_is_solved(i):
                solve_point_inefficient(self.postCompStates, i, self.preCompStates, i, 1 / self.compressorEfficiencies[i])
            elif self.postCompStates.point_is_solved(i):
                solve_point_inefficient(self.preCompStates, i, self.postCompStates, i, self.compressorEfficiencies[i])

        for i in range(0, len(self.preTurbStates.T)):
            if self.preTurbStates.point_is_solved(i) and self.postTurbStates.point_is_solved(i):
                h_isen = cp.PropsSI("Hmass", "P", self.postTurbStates.P[i], "Smass", self.preTurbStates.s[i], self.preTurbStates.fluid[i].to_coolprop())
                self.turbineEfficiencies[i] = (self.postTurbStates.h[i] - self.preTurbStates.h[i]) / (h_isen - self.preTurbStates.h[i])
                print("Turbine efficiency overwritten to", self.turbineEfficiencies[i])
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

        print(fluid.data)
        print("fuel compression work:", self.externalFuelCompressionWork)
        print("fuel injection temperature:", fuelTemp)
        print("fuel mole fraction:", fuelMoleFraction)
        print("molar masses:", M_oxidiser, M_fuel)


    def propagate_gas_from_combustor(self):
        for i in range(0, len(self.preTurbStates.T)):
            self.preTurbStates.fluid[i] = copy.deepcopy(self.preTurbStates.fluid[0])
            self.postTurbStates.fluid[i] = copy.deepcopy(self.preTurbStates.fluid[0])


    def combust(self):
        reactants, massFracs, temps = self.heatingStates.fluid[1].to_cea(self.heatingStates.T)
        products = self.heatingStates.fluid[1].to_cea(self.heatingStates.T[1], products_from_reactants=True)[0]
        solver = cea.EqSolver(products, reactants=reactants, transport=True)
        solution = cea.EqSolution(solver)
        hc = reactants.calc_property(cea.ENTHALPY, massFracs, temps) / cea.R
        solver.solve(solution, massFracs, self.heatingStates.P[1], hc=hc, mode="HP")

        
        #solution = cea.EqSolution(solver, pressure=self.heatingStates.P[1], mode='hp')
        self.preTurbStates.T[0] = solution.temperature
        products = Mixture()
        products.from_cea(solution)
        self.preTurbStates.fluid[0] = products

        self.propagate_gas_from_combustor()
        self.solve_cp()
        
        
        

    def q_in(self):
        q = abs(self.preTurbStates.h[0] - self.postCompStates.h[-1])
        for i in range(1, len(self.preTurbStates.T)): # reheating q_in
            q += abs(self.preTurbStates.h[i] - self.postTurbStates.h[i-1])
        return q

    def q_out(self):
        q = abs(self.preCompStates.h[0] - self.postTurbStates.h[-1])
        for i in range(1, len(self.preCompStates.T)): # intercooling q_out
            q += abs(self.preCompStates.h[i] - self.postCompStates.h[i-1])
        return q

    def w_in(self):
        w = 0
        for i in range(0, len(self.preCompStates.T)):
            w += abs(self.postCompStates.h[i] - self.preCompStates.h[i]) * (1 - self.addedFuelMassFrac)
        w += self.externalFuelCompressionWork
        return w

    def w_out(self):
        w = 0
        for i in range(0, len(self.preTurbStates.T)):
            w += abs(self.postTurbStates.h[i] - self.preTurbStates.h[i])
        return w

    def efficiency(self):
        return self.w_out() / (self.w_in() + self.q_in())



    def update_intermediates(self):
        self.intermediates = SetOfStates(0)
        
        for i in range(0, len(self.preCompStates.T)): # compressor
            self.intermediates.extend(get_intermediate_points(self.preCompStates, i, self.postCompStates, i, ["h", "s"]))

            if i+1 < len(self.preCompStates.T): # intercooling
                self.intermediates.extend(get_intermediate_points(self.postCompStates, i, self.preCompStates, i+1, ["T", "P"]))

        # boiler / combustion chamber
        self.intermediates.extend(get_intermediate_points(self.postCompStates, -1, self.preTurbStates, 0, ["h", "P"]))
        
        for i in range(0, len(self.preTurbStates.T)): # turbine
            self.intermediates.extend(get_intermediate_points(self.preTurbStates, i, self.postTurbStates, i, ["h", "s"]))

            if i+1 < len(self.preTurbStates.T): # reheating
                self.intermediates.extend(get_intermediate_points(self.postTurbStates, i, self.preTurbStates, i+1, ["T", "P"]))
        
        # condenser
        self.intermediates.extend(get_intermediate_points(self.postTurbStates, -1, self.preCompStates, 0, ["s", "P"]))

    
    def ts_diagram(self, ax=None, label=None, isRankine=False):
        allPoints = self.combined_states()
        
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
    
    def pv_diagram(self, ax=None, label=None, isRankine=False):
        allPoints = self.combined_states()
        
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

        """if isRankine:
            original_xlim = ax.get_xlim()
            original_ylim = ax.get_ylim()
            propertyPlot = PropertyPlot(self.postCompStates.fluid[0].to_coolprop(), 'PV', axis=ax, unit_system='SI')
            propertyPlot.calc_isolines(CoolProp.iQ, num=2)
            plt.close(propertyPlot.figure)
            propertyPlot.draw()
            ax.set_xlim(original_xlim)
            ax.set_ylim(original_ylim)
            ax.autoscale(True)"""
        
        if label != None:
            ax.legend()
        plt.show(block=False)
        return ax




testRankine = Cycle(1, 1, Mixture([["H2O", 1]]))
testRankine.preCompStates.P = [8000]
testRankine.preCompStates.X = [0]
testRankine.postCompStates.P = [8000000]
testRankine.preTurbStates.T = [600]
testRankine.preTurbStates.P = [8000000]
testRankine.postTurbStates.P = [8000]
testRankine.turbineEfficiencies[0] = 1
testRankine.compressorEfficiencies[0] = 0.7

testRankine.solve_cp()
print(testRankine.efficiency())
rankineTS = testRankine.ts_diagram(label="unmodified", isRankine=True)
rankinePV = testRankine.pv_diagram(label="unmodified", isRankine=True)

rankineTS.set_title("Rankine Test Plot Ts")
rankinePV.set_title("Rankine Test Plot Pv")

testRankine = Cycle(1, 2, Mixture([["H2O", 1]]))
testRankine.preCompStates.P = [8000]
testRankine.preCompStates.X = [0]
testRankine.postCompStates.P = [8000000]
testRankine.preTurbStates.T = [600, 600]
testRankine.preTurbStates.P = [8000000, 4000000]
testRankine.postTurbStates.P = [4000000, 8000]
testRankine.turbineEfficiencies = [0.7, 0.7]

testRankine.solve_cp()
print(testRankine.efficiency())
testRankine.ts_diagram(ax=rankineTS, label = "reheater")
testRankine.pv_diagram(ax=rankinePV, label = "reheater")


testBrayton = Cycle(1, 1, Mixture([["O2", 1]]))
testBrayton.preCompStates.P = [101325]
testBrayton.postCompStates.P = [500000]
testBrayton.preCompStates.T = [298]
#testBrayton.preTurbStates.T = [600]
testBrayton.preTurbStates.P = [500000]
testBrayton.postTurbStates.P = [101325]
#testBrayton.turbineEfficiencies = [0.9]
#testBrayton.compressorEfficiencies = [0.9]

testBrayton.calculate_pre_combustion_heat_exchange(30000)
testBrayton.inject_fuel_after_compressor_no_propagation(Mixture([["CH4", 1]]), 0.1, 0.9)
testBrayton.combust()

testBrayton.solve_cp()
braytonTS = testBrayton.ts_diagram(label="isentropic")
braytonPV = testBrayton.pv_diagram(label="isentropic")
print(testBrayton.efficiency())

braytonTS.set_title("Brayton Test Plot Ts")
braytonPV.set_title("Brayton Test Plot Pv")


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





