import cea
from CoolProp import CoolProp as cp # haha very funny
from CoolProp.Plots import PropertyPlot
import CoolProp
import matplotlib.pyplot as plt
import math
import copy
import numpy as np
import multiprocessing

STATE_POINT_SIZE = 80

def parse(text):
    if text == "None":
        return None

    f = float(text)
    i = int(f)
    if f == i:
        return i
    return f

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
    "C2H4" : "Ethylene",
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

    def to_save_str(self):
        text = ""
        for part in self.data:
            text += part[0] + "," + str(part[1])
            if len(part) > 2:
                text += "," + str(part[2])
            text += "+"
        return text[:-1]

    def from_save_str(self, text):
        self.data = []
        text = text.split("+")
        for line in text:
            line = line.split(",")
            line[1] = float(line[1])
            if len(line) > 2:
                line[2] = float(line[2])
            self.data.append(line)

def AIR():
    return Mixture([["N2", 0.78084], ["O2", 0.20946], ["Ar", 0.00934], ["CO2", 0.000412], ["Ne", 0.00001818]])
    # from https://www.engineeringtoolbox.com/air-composition-d_212.html

class SetOfStates():
    def __init__(self, numStates, fluid = None):
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

    def pop(self, i):
        if i == None:
            return
        self.T.pop(i)
        self.P.pop(i)
        self.v.pop(i)
        self.h.pop(i)
        self.s.pop(i)
        self.X.pop(i)
        self.fluid.pop(i)

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

    def solve_points(self, usePengRobinson=False):
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

            fluid = self.fluid[i].to_coolprop()
            if usePengRobinson:
                fluid = fluid.replace("HEOS::", "PR::")

            while len(knowns) >= 2:
                try:
                #if True:
                    if self.T[i] == None:
                        self.T[i] = cp.PropsSI("T", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], fluid)
                    if self.P[i] == None:
                        self.P[i] = cp.PropsSI("P", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], fluid)
                    if self.v[i] == None:
                        self.v[i] = 1 / cp.PropsSI("D", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], fluid)
                    if self.h[i] == None:
                        self.h[i] = cp.PropsSI("Hmass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], fluid)
                    if self.s[i] == None:
                        self.s[i] = cp.PropsSI("Smass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], fluid)
                    if self.X[i] == None:
                        self.X[i] = cp.PropsSI("Q", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], fluid)
                    knowns = []
                except:
                    knowns.pop(0)
    
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

    def to_save_str(self):
        text = "T="
        for val in self.T:
            text += str(val) + ","
        text = text[:-1] + "\nP="
        for val in self.P:
            text += str(val) + ","
        text = text[:-1] + "\nv="
        for val in self.v:
            text += str(val) + ","
        text = text[:-1] + "\nh="
        for val in self.h:
            text += str(val) + ","
        text = text[:-1] + "\ns="
        for val in self.s:
            text += str(val) + ","
        text = text[:-1] + "\nX="
        for val in self.X:
            text += str(val) + ","

        text = text[:-1] + "\nfluid="
        for fluid in self.fluid:
            text += fluid.to_save_str() + ";"
        
        return text[:-1]

    def from_save_str(self, text):
        for line in text.split("\n"):
            if line == "":
                continue
            prop, data = line.split("=")
            data = data.split(",")

            if prop != "fluid":
                for i in range(0, len(data)):
                    data[i] = parse(data[i])
            
            if prop == "T":
                self.T = data
            elif prop == "P":
                self.P = data
            elif prop == "v":
                self.v = data
            elif prop == "h":
                self.h = data
            elif prop == "s":
                self.s = data
            elif prop == "X":
                self.X = data

            elif prop == "fluid":
                self.fluid = []
                data = ",".join(data)
                data = data.split(";")
                for f in data:
                    mix = Mixture()
                    mix.from_save_str(f)
                    self.fluid.append(mix)

            else:
                print("ERROR: unknown property '", prop, "' from save file")
        





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
    #if toSet.X[i] != None:
    #    knowns.append(["Q", toSet.X[i]])


    if len(knowns) == 0:
        #print("Couldn't solve point", i, "in solve_point_inefficient")
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









def load_cycle(path):
    if not path.lower().endswith(".ppdc"):
        path += ".ppdc" # .ppdc = 'PowerPlant Designer Cycle'

    file = open(path, "r")
    data = file.read().split("\n<><>\n")
    file.close()

    misc = data[0]
    states = data[1:]

    cycle = Cycle(0, 0, None)

    for line in misc.split("\n"):
        prop, val = line.split("=")

        if prop == "externalFuelCompressionWork":
            cycle.externalFuelCompressionWork = parse(val)
        elif prop == "addedFuelMassFrac":
            cycle.addedFuelMassFrac = parse(val)
        elif prop == "addedFuelInitialEnthalpy":
            cycle.addedFuelInitialEnthalpy = parse(val)
        elif prop == "recycledHeat":
            cycle.recycledHeat = parse(val)
        elif prop == "closedFeedwaterMassRatio":
            cycle.closedFeedwaterMassRatio = parse(val)
        elif prop == "turbineMassRatios":
            cycle.turbineMassRatios = []
            for m in val.split(","):
                cycle.turbineMassRatios.append(parse(m))
        elif prop == "closedFeedwaterExitEnthalpy":
            cycle.closedFeedwaterExitEnthalpy = parse(val)
        elif prop == "closedFeedwaterHeatTransfer":
            cycle.closedFeedwaterHeatTransfer = parse(val)
        elif prop == "closedFeedwaterTurbineIndex":
            cycle.closedFeedwaterTurbineIndex = parse(val)
        elif prop == "turbineEfficiencies":
            cycle.turbineEfficiencies = []
            for e in val.split(","):
                cycle.turbineEfficiencies.append(parse(e))
        elif prop == "compressorEfficiencies":
            cycle.compressorEfficiencies = []
            for e in val.split(","):
                cycle.compressorEfficiencies.append(parse(e))

    for stateList in states:
        name, data = stateList.split("<>\n")
        if name == "preCompStates":
            cycle.preCompStates.from_save_str(data)
        if name == "postCompStates":
            cycle.postCompStates.from_save_str(data)
        if name == "heatingStates":
            cycle.heatingStates.from_save_str(data)
        if name == "preTurbStates":
            cycle.preTurbStates.from_save_str(data)
        if name == "postTurbStates":
            cycle.postTurbStates.from_save_str(data)
        if name == "coolingStates":
            cycle.coolingStates.from_save_str(data)



    return cycle



class Cycle: # should be valid for both rankine and brayton
    def __init__(self, numComp, numTurb, fluid):
        self.preCompStates = SetOfStates(numComp, copy.deepcopy(fluid))
        self.postCompStates = SetOfStates(numComp, copy.deepcopy(fluid))
        self.heatingStates = SetOfStates(2, copy.deepcopy(fluid)) # heatingStates = [before fuel injected, before combustion but after fuel injected]
        self.preTurbStates = SetOfStates(numTurb, copy.deepcopy(fluid))
        self.postTurbStates = SetOfStates(numTurb, copy.deepcopy(fluid))
        self.coolingStates = SetOfStates(0)
        self.turbineEfficiencies = [1] * numTurb
        self.compressorEfficiencies = [1] * numComp
        
        self.externalFuelCompressionWork = 0 # work to pressurize fuel before injecting it into combustion chamber
        self.addedFuelMassFrac = 0
        self.addedFuelInitialEnthalpy = 0
        self.recycledHeat = 0
        self.closedFeedwaterMassRatio = 0
        self.turbineMassRatios = [1] * numTurb # fraction left at each stage after removal from feedwater heater
        self.closedFeedwaterExitEnthalpy = 0
        self.closedFeedwaterHeatTransfer = 0
        self.closedFeedwaterTurbineIndex = None

        self.allStateSets = [self.preCompStates, self.postCompStates, self.preTurbStates, self.postTurbStates]
        self.intermediates = None

    def copy(self):
        return copy.deepcopy(self)

    def save(self, path):
        if not path.lower().endswith(".ppdc"):
            path += ".ppdc" # .ppdc = 'Powerplant Designer Cycle'
        data = ""

        data += "externalFuelCompressionWork=" + str(self.externalFuelCompressionWork) + "\n"
        data += "addedFuelMassFrac=" + str(self.addedFuelMassFrac) + "\n"
        data += "addedFuelInitialEnthalpy=" + str(self.addedFuelInitialEnthalpy) + "\n"
        data += "recycledHeat=" + str(self.recycledHeat) + "\n"
        data += "closedFeedwaterMassRatio=" + str(self.closedFeedwaterMassRatio) + "\n"
        data += "turbineMassRatios="
        for r in self.turbineMassRatios:
            data += str(r) + ","
        data = data[:-1] + "\n"
        data += "closedFeedwaterExitEnthalpy=" + str(self.closedFeedwaterExitEnthalpy) + "\n"
        data += "closedFeedwaterHeatTransfer=" + str(self.closedFeedwaterHeatTransfer) + "\n"
        data += "closedFeedwaterTurbineIndex=" + str(self.closedFeedwaterTurbineIndex) + "\n"
        data += "turbineEfficiencies="
        for e in self.turbineEfficiencies:
            data += str(e) + ","
        data = data[:-1] + "\n"
        data += "compressorEfficiencies="
        for e in self.compressorEfficiencies:
            data += str(e) + ","
        data = data[:-1] + "\n"

        data += "<><>\npreCompStates<>\n"
        data += self.preCompStates.to_save_str() + "\n"
        data += "<><>\npostCompStates<>\n"
        data += self.postCompStates.to_save_str() + "\n"
        data += "<><>\nheatingStates<>\n"
        data += self.heatingStates.to_save_str() + "\n"
        data += "<><>\npreTurbStates<>\n"
        data += self.preTurbStates.to_save_str() + "\n"
        data += "<><>\npostTurbStates<>\n"
        data += self.postTurbStates.to_save_str() + "\n"

        if len(self.coolingStates.T) > 0:
            data += "<><>\ncoolingStates<>\n"
            data += self.coolingStates.to_save_str() + "\n"

        file = open(path, "w")
        file.write(data[:-1])
        file.close()
        
        

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

        if self.closedFeedwaterHeatTransfer != 0:
            states.extend(self.feedwater_state())

        if self.heatingStates.set_is_fully_solved():
            states.extend(self.heatingStates)
            
        for i in range(0, len(self.preTurbStates.T)):
            states.extend(self.preTurbStates.get_point(i))
            states.extend(self.postTurbStates.get_point(i))
            
        states.extend(self.coolingStates)
        
        return states # returned object is READ ONLY

    def feedwater_state(self):
        fwState = SetOfStates(1)
        fwState.h[0] = self.postCompStates.h[-1] + self.closedFeedwaterHeatTransfer
        fwState.P[0] = self.postCompStates.P[-1]
        fwState.fluid[0] = self.postCompStates.fluid[-1]
        fwState.solve_points()
        return fwState

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

    def add_pre_combustion_heat_exchange(self, q_in):
        self.solve_cp()
        
        if self.heatingStates.h[0] != None:
            q_in += self.heatingStates.h[0] - self.postCompStates.h[-1]
        
        self.heatingStates.P[0] = self.postCompStates.P[-1]
        self.heatingStates.h[0] = self.postCompStates.h[-1] + q_in
        self.heatingStates.fluid[0] = self.postCompStates.fluid[-1]
        self.heatingStates.P[1] = self.heatingStates.P[0]
        self.heatingStates.h[1] = self.heatingStates.h[0]
        self.heatingStates.fluid[1] = self.heatingStates.fluid[0]
        self.heatingStates.solve_points()
        

    def inject_fuel_after_compressor(self, fuelMixture, fuelTotalMassFraction, fuelCompressorEfficiency):
        self.addedFuelMassFrac = fuelTotalMassFraction
        if not self.heatingStates.set_is_fully_solved():
            self.add_pre_combustion_heat_exchange(0)
            print("Assumed no pre-combustion heat exchange when injecting fuel")

        # compressing fuel from standard atmospheric conditions to same pressure as air at exit of compressor
        # then finding work done and final fuel temperature during injection
        h0 = cp.PropsSI("Hmass", "P", 101325, "T", 298, fuelMixture.to_coolprop())
        s0 = cp.PropsSI("Smass", "P", 101325, "T", 298, fuelMixture.to_coolprop())
        h1 = cp.PropsSI("Hmass", "P", self.heatingStates.P[0], "S", s0, fuelMixture.to_coolprop())
        h1 = h0 + (h1 - h0) / fuelCompressorEfficiency
        fuelTemp = cp.PropsSI("T", "P", self.heatingStates.P[0], "Hmass", h1, fuelMixture.to_coolprop())
        self.externalFuelCompressionWork = (h1 - h0) * fuelTotalMassFraction
        self.addedFuelInitialEnthalpy = h0

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

    def set_closed_feedwater_heater_from_fraction(self, fraction, turbineIndex, removeReheater=True):
        for i in range(turbineIndex + 1, len(self.preTurbStates.T)):
            self.turbineMassRatios[i] -= fraction
            #if self.turbineMassRatios[i] <= 0:
                #print("ERROR: turbine", i, "mass ratio =", self.turbineMassRatios[i])

        self.solve_cp()

        if removeReheater:
            self.preTurbStates.T[turbineIndex+1] = self.postTurbStates.T[turbineIndex]
            self.preTurbStates.P[turbineIndex+1] = self.postTurbStates.P[turbineIndex]
            self.preTurbStates.v[turbineIndex+1] = self.postTurbStates.v[turbineIndex]
            self.preTurbStates.h[turbineIndex+1] = self.postTurbStates.h[turbineIndex]
            self.preTurbStates.s[turbineIndex+1] = self.postTurbStates.s[turbineIndex]
            self.preTurbStates.X[turbineIndex+1] = self.postTurbStates.X[turbineIndex]

        
        P = self.postTurbStates.P[turbineIndex]
        h2 = cp.PropsSI("Hmass", "P", P, "Q", 0, self.postTurbStates.fluid[turbineIndex].to_coolprop())
        heat = (self.postTurbStates.h[turbineIndex] - h2) * fraction
        #print("Closed Feedwater Heat:", heat)
        self.recycledHeat += heat
        self.closedFeedwaterHeatTransfer = heat
        self.closedFeedwaterExitEnthalpy = h2
        self.closedFeedwaterMassRatio = fraction
        self.closedFeedwaterTurbineIndex = turbineIndex

    def set_closed_feedwater_heater_from_heat(self, heat, turbineIndex, removeReheater=True):
        self.solve_cp()
        P = self.postTurbStates.P[turbineIndex]
        h2 = cp.PropsSI("Hmass", "P", P, "Q", 0, self.postTurbStates.fluid[turbineIndex].to_coolprop())
        fraction = heat / (self.postTurbStates.h[turbineIndex] - h2)
        #print("Closed Feedwater Fraction:", fraction)
        self.set_closed_feedwater_heater_from_fraction(fraction, turbineIndex, removeReheater=removeReheater)

    def set_closed_feedwater_heater_from_final_temp(self, temp, turbineIndex, removeReheater=True):
        self.solve_cp()
        h1 = self.postCompStates.h[-1]
        h2 = cp.PropsSI("Hmass", "P", self.postCompStates.P[-1], "T", temp, self.postCompStates.fluid[-1].to_coolprop())
        self.set_closed_feedwater_heater_from_heat(h2 - h1, turbineIndex, removeReheater=removeReheater)

    def set_closed_feedwater_heater_from_temp_diff(self, tempDiff, turbineIndex, removeReheater=True):
        self.solve_cp()
        #temp = self.postTurbStates.T[turbineIndex] - temp
        P = self.postTurbStates.P[turbineIndex]
        temp = cp.PropsSI("T", "P", P, "Q", 0, self.postTurbStates.fluid[turbineIndex].to_coolprop()) - tempDiff
        self.set_closed_feedwater_heater_from_final_temp(temp, turbineIndex, removeReheater=removeReheater)
        

    def propagate_assuming_no_pressure_losses(self):
        if self.preCompStates.P[0] == None:
            self.preCompStates.P[0] = self.postTurbStates.P[-1]
        if self.postTurbStates.P[-1] == None:
            self.postTurbStates.P[-1] = self.preCompStates.P[0]

        if self.preTurbStates.P[0] == None:
            self.preTurbStates.P[0] = self.postCompStates.P[-1]
        if self.postCompStates.P[-1] == None:
            self.postCompStates.P[-1] = self.preTurbStates.P[0]

        for i in range(1, len(self.preCompStates.P)):
            if self.preCompStates.P[i] == None:
                self.preCompStates.P[i] = self.postCompStates.P[i-1]
            if self.postCompStates.P[i-1] == None:
                self.postCompStates.P[i-1] = self.preCompStates.P[i]

        for i in range(1, len(self.preTurbStates.P)):
            if self.preTurbStates.P[i] == None:
                self.preTurbStates.P[i] = self.postTurbStates.P[i-1]
            if self.postTurbStates.P[i-1] == None:
                self.postTurbStates.P[i-1] = self.preTurbStates.P[i]

    def solve_rankine_reheater_pressures(self, untilTurbineIndex, targetQuality):
        if untilTurbineIndex == -1:
            untilTurbineIndex = len(self.preTurbStates.P) - 1
        
        for j in range(0, untilTurbineIndex):
            if self.postTurbStates.P[j] == None:
                reheaterPressure = optimize_turbine_exit_pressure_for_quality(self, j, targetQuality)
                self.postTurbStates.P[j] = reheaterPressure
                self.propagate_assuming_no_pressure_losses()
                self.solve_cp()

    def remove_reheat_step(self, index): # also works for feedwater exit
        self.preTurbStates.pop(index + 1)
        self.postTurbStates.pop(index)

        # merging adjecent turbine efficiencies with average, if both have a value
        if self.turbineEfficiencies[index] != None and self.turbineEfficiencies[index + 1] != None:
            self.turbineEfficiencies[index] = (self.turbineEfficiencies[index] + self.turbineEfficiencies[index + 1]) / 2
        elif self.turbineEfficiencies[index + 1] != None:
            self.turbineEfficiencies[index] = self.turbineEfficiencies[index + 1]
            
        self.turbineEfficiencies.pop(1)
        



    def rankine_optimize_fixed_compression(self, turbineExitQualities, feedwaterArgs=None):
        for i in range(0, len(self.preTurbStates.P)):
            if self.preTurbStates.T[i] != None:
                if feedwaterArgs == None:
                    self.preTurbStates.T[i] = rankine_optimize_turbine_inlet_temperature(self, i, turbineExitQualities, feedwaterArgs=feedwaterArgs)
                elif self.preTurbStates.T[i] != feedwaterArgs[1]:
                    self.preTurbStates.T[i] = rankine_optimize_turbine_inlet_temperature(self, i, turbineExitQualities, feedwaterArgs=feedwaterArgs)
        
        self.rankine_optimize_no_temp_changes(turbineExitQualities, feedwaterArgs=feedwaterArgs)

        

    def rankine_optimize_no_temp_changes(self, turbineExitQualities, feedwaterArgs=None):
        if feedwaterArgs != None:
            self.solve_rankine_reheater_pressures(feedwaterArgs[1], turbineExitQualities)
            self.propagate_assuming_no_pressure_losses()
            self.solve_cp()

            P = optimize_feedwater_pressure(self, feedwaterArgs[0], feedwaterArgs[1], removeReheater=feedwaterArgs[2])
            self.postTurbStates.P[1] = P
            self.set_closed_feedwater_heater_from_temp_diff(feedwaterArgs[0], feedwaterArgs[1], removeReheater=feedwaterArgs[2])
            self.solve_cp()
        
        self.solve_rankine_reheater_pressures(-1, turbineExitQualities)

    def rankine_optimize_fixed_target_quality(self, turbineExitQualities, feedwaterArgs=None):
        boilerPressure = rankine_optimize_boiler_pressure(self, turbineExitQualities, feedwaterArgs=feedwaterArgs)
        self.postCompStates.P[-1] = boilerPressure
        self.propagate_assuming_no_pressure_losses()
        self.rankine_optimize_fixed_compression(turbineExitQualities, feedwaterArgs=feedwaterArgs)

    def rankine_optimize(self, turbineExitQualities, feedwaterArgs=None):
        turbineExitQualities = rankine_optimize_target_quality(self, turbineExitQualities, feedwaterArgs=feedwaterArgs)
        self.rankine_optimize_fixed_target_quality(turbineExitQualities, feedwaterArgs=feedwaterArgs)




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
            w.append((self.preTurbStates.h[i] - self.postTurbStates.h[i]) * self.turbineMassRatios[i])
        return w

    def reheater_heats(self):
        self.preTurbStates.make_consistent(True)
        self.postTurbStates.make_consistent(True)
        q = []
        for i in range(1, len(self.preTurbStates.T)):
            q.append((self.preTurbStates.h[i] - self.postTurbStates.h[i-1]) * self.turbineMassRatios[i])
        return q

    def return_heat(self):
        self.postTurbStates.make_consistent(True)
        self.preCompStates.make_consistent(True)
        fuel = (self.postTurbStates.h[-1] - self.addedFuelInitialEnthalpy) * self.addedFuelMassFrac * self.turbineMassRatios[-1]
        air = (self.postTurbStates.h[-1] - self.preCompStates.h[0]) * (1 - self.addedFuelMassFrac) * self.turbineMassRatios[-1]
        feedwater = (self.closedFeedwaterExitEnthalpy - self.preCompStates.h[0]) * self.closedFeedwaterMassRatio
        
        return fuel + air + feedwater

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
            fuel = (self.preTurbStates.h[0] - self.addedFuelInitialEnthalpy - self.externalFuelCompressionWork) * self.addedFuelMassFrac
            air = (self.preTurbStates.h[0] - self.heatingStates.h[0]) * (1 - self.addedFuelMassFrac)
            return fuel + air

    def get_q_in_and_temps(self): # for combined cycle calculations
        out = [[self.preTurbStates.T[0], self.postCompStates.T[-1], self.heat_exchanger_heat() - self.recycledHeat]]
        if self.closedFeedwaterHeatTransfer != 0:
            out = [[self.preTurbStates.T[0], self.feedwater_state().T[-1], self.heat_exchanger_heat() - self.recycledHeat]]
        
        for i in range(1, len(self.preTurbStates.T)):
            q = (self.preTurbStates.h[i] - self.postTurbStates.h[i-1]) * self.turbineMassRatios[i]
            T1 = self.postTurbStates.T[i-1]
            T2 = self.preTurbStates.T[i]
            if q != 0:
                out.append([T2, T1, q])

        return sorted(out, reverse=True)

    def q_in(self): # assuming heat exchanger source is external, not recycling heat
        return self.heat_exchanger_heat() + self.combustion_heat() + sum(self.reheater_heats())

    def q_out(self): # assuming no heat is recycled using a heat exchanger
        return sum(self.intercooler_heats()) + self.return_heat()

    def w_in(self):
        return sum(self.compressor_works()) + self.externalFuelCompressionWork

    def w_out(self):
        return sum(self.turbine_works())

    def w_net(self):
        if self.closedFeedwaterMassRatio >= 1:
            return 0
        return self.w_out() - self.w_in()

    def efficiency(self):
        return self.w_net() / (self.q_in() - self.recycledHeat)

    def back_work(self):
        return self.w_in() / self.w_out()

    def efficiency_2nd_law(self):
        T_H = max(self.preTurbStates.T)
        T_C = min(self.preCompStates.T)
        carnot = 1 - T_C / T_H
        return self.efficiency() / carnot

    def fuel_efficiency(self):
        return self.w_net() / self.addedFuelMassFrac

    def print_design_info(self):
        print()
        print("---------------------------------------------------------------")
        print()
        print("Lowest pressure:", self.preCompStates.P[0], "Pa")
        print()
        for i in range(0, len(self.postCompStates.P)):
            print("Compressor", i, "pressure ratio:", self.postCompStates.P[i] / self.preCompStates.P[i])
        print()
        print("Highest pressure:", self.postCompStates.P[-1], "Pa")
        print()
        for i in range(0, len(self.postTurbStates.P)):
            print("Turbine", i, "pressure ratio:", self.preTurbStates.P[i] / self.postTurbStates.P[i])
            if self.postTurbStates.X[i] == -1:
                print("          exit quality: superheated vapour")
            else:
                print("          exit quality:", self.postTurbStates.X[i])
        if self.closedFeedwaterMassRatio > 0:
            print()
            print("Closed feedwater heater flow ratio:", self.closedFeedwaterMassRatio)
            print("Closed feedwater heater heat transfer:", self.closedFeedwaterHeatTransfer, "J/kg")
            fwState = SetOfStates(1)
            fwState.h[0] = self.postCompStates.h[-1] + self.closedFeedwaterHeatTransfer
            fwState.P[0] = self.postCompStates.P[-1]
            fwState.fluid[0] = self.postCompStates.fluid[-1]
            fwState.solve_points()
            print("Temperature after feedwater heat exchange:", fwState.T[0], "K")
        if self.addedFuelMassFrac != 0:
            print()
            print("Injected fuel mass fraction:", self.addedFuelMassFrac)
        print()
        print("---------------------------------------------------------------")
        print()
        

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
        print("2nd Law Efficiency:", self.efficiency_2nd_law())
        if self.addedFuelMassFrac != 0:
            print("Fuel Efficiency (w_net / m_fuel):", self.fuel_efficiency())
        print("Back-work ratio:", self.back_work())
        if self.recycledHeat > self.heat_exchanger_heat():
            print()
            print("WARNING: recycled heat is greater than heat entering through heat exchanger")
        if self.closedFeedwaterMassRatio >= 1:
            print()
            print("WARNING: feedwater mass fraction > 1")
        if self.efficiency_2nd_law() > 1:
            print()
            print("WARNING: somehow beat Carnot efficiency")
        print()
        print("---------------------------------------------------------------")
        print()





    def print_all_states(self):
        allPoints = self.combined_states()
        allPoints.make_consistent(True)
        if self.closedFeedwaterTurbineIndex != None:
            feedwaterPoints = self.get_closed_feedwater_states()
            feedwaterPoints.pop(0)
            allPoints.extend(feedwaterPoints)
        
        print()
        print("---------------------------------------------------------------")
        print()
        print("All cycle states:")
        print()
        for i in range(0, len(allPoints.T)):
            print()
            print("State", i)
            print("T =", allPoints.T[i], "K")
            print("P =", allPoints.P[i], "Pa")
            print("v =", allPoints.v[i], "m^3/kg")
            print("h =", allPoints.h[i], "kJ/kg")
            print("s =", allPoints.s[i], "kJ/kgK")
            if allPoints.X[i] != -1:
                print("X =", allPoints.X[i])
            else:
                print("X = superheated or subcooled")
            print("fluid =", allPoints.fluid[i].to_save_str())
        
        print()
        print("---------------------------------------------------------------")
        print()





    def update_intermediates(self):
        self.intermediates = SetOfStates(0)
        
        for i in range(0, len(self.preCompStates.T)): # compressor
            self.intermediates.extend(get_intermediate_points(self.preCompStates, i, self.postCompStates, i, ["h", "s"]))

            if i+1 < len(self.preCompStates.T): # intercooling
                self.intermediates.extend(get_intermediate_points(self.postCompStates, i, self.preCompStates, i+1, ["h", "P"]))

        # boiler / combustion chamber
        
        if self.heatingStates.h[0] != None and self.heatingStates.h[0] != self.postTurbStates.h[-1]: # if heat exchange before combustion
            # heat exchanger process:
            self.intermediates.extend(get_intermediate_points(self.postCompStates, -1, self.heatingStates, 0, ["s", "P"]))
            # injection process:
            self.intermediates.extend(get_intermediate_points(self.heatingStates, 0, self.heatingStates, 1, ["T", "P"], numPoints=2))
            # combustion process:
            self.intermediates.extend(get_intermediate_points(self.heatingStates, 1, self.preTurbStates, 0, ["T", "P"]))
        
        else: # single process combustion (no heat exchangers)
            try:
                self.intermediates.extend(get_intermediate_points(self.postCompStates, -1, self.preTurbStates, 0, ["s", "P"]))
            except:
                self.intermediates.extend(get_intermediate_points(self.postCompStates, -1, self.preTurbStates, 0, ["T", "P"]))
        
        
        for i in range(0, len(self.preTurbStates.T)): # turbine
            self.intermediates.extend(get_intermediate_points(self.preTurbStates, i, self.postTurbStates, i, ["h", "s"]))

            if i+1 < len(self.preTurbStates.T): # reheating
                if self.postTurbStates.T[i] != self.preTurbStates.T[i+1]:
                    try:
                        self.intermediates.extend(get_intermediate_points(self.postTurbStates, i, self.preTurbStates, i+1, ["s", "P"]))
                    except:
                        self.intermediates.extend(get_intermediate_points(self.postTurbStates, i, self.preTurbStates, i+1, ["T", "P"]))


        

        if self.postTurbStates.fluid[-1].data == self.preCompStates.fluid[0].data: # single step return to compressor
            try:
                self.intermediates.extend(get_intermediate_points(self.postTurbStates, -1, self.preCompStates, 0, ["s", "P"]))
            except:
                self.intermediates.extend(get_intermediate_points(self.postTurbStates, -1, self.preCompStates, 0, ["T", "P"]))
        
        elif len(self.coolingStates.T) > 0: # brayton cooling processes in heat exchangers
            try:
                self.intermediates.extend(get_intermediate_points(self.postTurbStates, -1, self.coolingStates, -1, ["s", "P"]))
            except:
                self.intermediates.extend(get_intermediate_points(self.postTurbStates, -1, self.coolingStates, -1, ["T", "P"]))

        self.intermediates.make_consistent(True)


    def get_closed_feedwater_states(self):
        i = self.closedFeedwaterTurbineIndex
        feedwaterStates = SetOfStates(3, self.postTurbStates.fluid[i])
        
        feedwaterStates.P[0] = self.postTurbStates.P[i]
        feedwaterStates.T[0] = self.postTurbStates.T[i]
        feedwaterStates.v[0] = self.postTurbStates.v[i]
        feedwaterStates.h[0] = self.postTurbStates.h[i]
        feedwaterStates.s[0] = self.postTurbStates.s[i]
        feedwaterStates.X[0] = self.postTurbStates.X[i]
        
        feedwaterStates.P[1] = self.postTurbStates.P[i]
        feedwaterStates.X[1] = 0

        feedwaterStates.P[2] = self.postTurbStates.P[-1]
        feedwaterStates.h[2] = self.closedFeedwaterExitEnthalpy
        
        feedwaterStates.solve_points()
        return feedwaterStates

    
    def closed_feedwater_ts(self, ax):
        feedwaterStates = self.get_closed_feedwater_states()
        feedwaterStates.solve_points()
        intermediates = SetOfStates(3, feedwaterStates.fluid[0])
        intermediates.extend(get_intermediate_points(feedwaterStates, 0, feedwaterStates, 1, ["P", "s"]))
        intermediates.extend(get_intermediate_points(feedwaterStates, 1, feedwaterStates, 2, ["h", "P"], numPoints=5))

        ax.scatter(feedwaterStates.s, feedwaterStates.T, marker="o", s=STATE_POINT_SIZE)
        ax.plot(intermediates.s, intermediates.T, marker=".")
        for i in range(1, len(feedwaterStates.s)):
            ax.text(feedwaterStates.s[i], feedwaterStates.T[i], str(i-1 + len(self.combined_states().T)), fontsize=20)
        

    
    def ts_diagram(self, ax=None, label=None, isRankine=False):
        allPoints = self.combined_states()
        allPoints.make_consistent(True)
        
        if ax is None:
            ax = plt.figure().add_subplot()
            ax.set_xlabel("Specific Entropy - $J kg^{-1} K^{-1}$")
            ax.set_ylabel("Temperature - $K$")

        x = []
        y = []
        stateLabels = []
        for i in range(0, len(allPoints.T)):
            x.append(allPoints.s[i])
            y.append(allPoints.T[i])
            stateLabels.append(str(i))

        ax.scatter(x, y, label=label, marker="o", s=STATE_POINT_SIZE)

        for i in range(0, len(x)):
            ax.text(x[i], y[i], stateLabels[i], fontsize=20)

        
        
        if self.intermediates == None:
            self.update_intermediates()

        x = self.intermediates.s
        y = self.intermediates.T
        
        ax.plot(x, y, marker=".")

        if isRankine:
            original_xlim = ax.get_xlim()
            original_ylim = ax.get_ylim()
            original_xlabel = ax.get_xlabel()
            original_ylabel = ax.get_ylabel()
            propertyPlot = PropertyPlot(self.postCompStates.fluid[0].to_coolprop(), 'TS', axis=ax, unit_system='SI')
            propertyPlot.calc_isolines(CoolProp.iQ, num=2)
            plt.close(propertyPlot.figure)
            propertyPlot.draw()
            ax.set_xlim(original_xlim)
            ax.set_ylim(original_ylim)
            ax.set_xlabel(original_xlabel)
            ax.set_ylabel(original_ylabel)
            ax.autoscale(True)

        if self.closedFeedwaterTurbineIndex != None:
            self.closed_feedwater_ts(ax)
        
        if label != None:
            ax.legend()
        plt.show(block=False)
        return ax
    
    def pv_diagram(self, ax=None, label=None):
        allPoints = self.combined_states()
        allPoints.make_consistent(True)
        
        if ax is None:
            ax = plt.figure().add_subplot()
            ax.set_xlabel("Specific Volume - $m^3 kg^{-1}$")
            ax.set_ylabel("Pressure - $Pa$")

        x = []
        y = []
        stateLabels = []
        for i in range(0, len(allPoints.T)):
            x.append(allPoints.v[i])
            y.append(allPoints.P[i])
            stateLabels.append(str(i))

        ax.scatter(x, y, label=label, marker="o", s=STATE_POINT_SIZE)

        for i in range(0, len(x)):
            ax.text(x[i], y[i], stateLabels[i], fontsize=20)

        
        
        if self.intermediates == None:
            self.update_intermediates()

        x = self.intermediates.v
        y = self.intermediates.P
        
        ax.plot(x, y, marker=".")
        
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


class CombinedBraytonRankineCycle:
    def __init__(self, brayton, rankine):
        self.brayton = brayton
        self.rankine = rankine
        self.braytonFlowRate = 1 # brayton flow rate / rankine flow rate
        self.transferredHeat = 0  # heat transefer from brayton to rankine

    def q_in(self): # assuming no heat recycling within each cycle
        return self.brayton.q_in() * self.braytonFlowRate + self.rankine.q_in() - self.transferredHeat

    def q_out(self): # assuming no heat recycling within each cycle
        return self.brayton.q_out() * self.braytonFlowRate + self.rankine.q_out() - self.transferredHeat

    def external_q_in(self):
        return self.q_in() - self.brayton.recycledHeat * self.braytonFlowRate - self.rankine.recycledHeat
    
    def external_q_out(self):
        return self.q_out() - self.brayton.recycledHeat * self.braytonFlowRate - self.rankine.recycledHeat

    def w_in(self):
        return self.brayton.w_in() * self.braytonFlowRate + self.rankine.w_in()

    def w_out(self):
        return self.brayton.w_out() * self.braytonFlowRate + self.rankine.w_out()

    def w_net(self):
        return self.w_out() - self.w_in()

    def back_work(self):
        return self.w_in() / self.w_out()

    def efficiency(self):
        return self.w_net() / (self.q_in() - self.brayton.recycledHeat * self.braytonFlowRate - self.rankine.recycledHeat)

    def fuel_efficiency(self):
        return self.w_net() / (self.brayton.addedFuelMassFrac * self.braytonFlowRate)

    def efficiency_2nd_law(self):
        hotTemps = self.brayton.preTurbStates.T.copy()
        hotTemps.extend(self.rankine.preTurbStates.T.copy())
        T_H = max(hotTemps)
        coldTemps = self.brayton.preCompStates.T.copy()
        coldTemps.extend(self.rankine.preCompStates.T.copy())
        T_C = min(coldTemps)
        carnot = 1 - T_C / T_H
        return self.efficiency() / carnot

    def add_brayton_cooling_states(self):
        transfers = self.rankine.get_q_in_and_temps()
        for i in range(0, len(transfers)):
            transfers[i][2] = transfers[i][2] / (self.braytonFlowRate * self.brayton.turbineMassRatios[-1])

        lastH = self.brayton.postTurbStates.h[-1]
        self.brayton.coolingStates = SetOfStates(len(transfers), self.brayton.postTurbStates.fluid[-1])
        for i in range(0, len(transfers)):
            lastH -= transfers[i][2]
            self.brayton.coolingStates.h[i] = lastH
            self.brayton.coolingStates.P[i] = self.brayton.postTurbStates.P[-1]
        
        self.brayton.coolingStates.solve_points()
    
    def adjust_brayton_flow_rate(self, deltaT, showGraph=False):
        if showGraph:
            plt.figure()
            plt.xlabel("iteration")
            plt.ylabel("Brayton Relative Mass Flow Rate")
            plt.title("Brayton Flow Rate Convergence for Combined Brayton-Rankine System")
        
        transfers = self.rankine.get_q_in_and_temps()
        m1 = 0
        m2 = 1
        iteration = 0
        expanded = False
        while not expanded:
            self.braytonFlowRate = m2
            self.add_brayton_cooling_states()
            expanded = True
            for i in range(0, len(transfers)):
                startTemp = self.brayton.postTurbStates.T[-1]
                if i > 0:
                    startTemp = self.brayton.coolingStates.T[i-1]
                
                if self.brayton.coolingStates.T[i] == None:
                    # failed to converge because cryogenic temperature, mass rate is obviously way too low
                    expanded = False
                    m1 = m2
                    m2 = m2 * 1.5
                
                elif transfers[i][0] + deltaT > startTemp:
                    if i == 0:
                        print("Brayton Turbine exit temperature is lower than rankine boiler temperature")
                        return False
                    expanded = False
                    m1 = m2
                    m2 = m2 * 1.5
                
                elif transfers[i][1] + deltaT > self.brayton.coolingStates.T[i]:
                    expanded = False
                    m1 = m2
                    m2 = m2 * 1.5

            if showGraph:
                iteration += 1
                plt.plot([iteration, iteration], [m1, m2])
                plt.show(block=False)
                #plt.pause(0.5)

        while m2 - m1 > 0.0001:
            m = (m1 + m2) / 2
            self.braytonFlowRate = m
            self.add_brayton_cooling_states()
            consistent = True
            for i in range(0, len(transfers)):
                startTemp = self.brayton.postTurbStates.T[-1]
                if i > 0:
                    startTemp = self.brayton.coolingStates.T[i-1]

                if transfers[i][0] + deltaT > startTemp:
                    consistent = False

                elif transfers[i][1] + deltaT > self.brayton.coolingStates.T[i]:
                    consistent = False

            if consistent:
                m2 = m
            else:
                m1 = m
            
            if showGraph:
                iteration += 1
                plt.plot([iteration, iteration], [m1, m2])
                plt.show(block=False)
                #plt.pause(0.5)

        self.transferredHeat = 0
        for part in transfers:
            self.transferredHeat += part[2]
        return m2

    def ts_diagram(self):
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
        fig.suptitle("Combined Cycle Ts Diagram")
        ax1.set_title("Brayton")
        ax2.set_title("Rankine")
        ax1.set_ylabel("Temperature - $K$")
        ax1.set_xlabel("Specific Entropy - $J kg^{-1} K^{-1}$")
        ax2.set_xlabel("Specific Entropy - $J kg^{-1} K^{-1}$")
        self.brayton.ts_diagram(ax=ax1)
        self.rankine.ts_diagram(ax=ax2, isRankine=True)
        
    def print_energy_calcs(self):
        print()
        print("---------------------------------------------------------------")
        print()
        print("Energy Calculations:")
        print("All units in J/kg noramalized using the mass of the rankine cycle")
        print()
        print()
        print("Brayton relative mass flow rate:", self.braytonFlowRate)
        print()
        print()
        print("Brayton q_in:", self.brayton.q_in() * self.braytonFlowRate)
        print("Brayton recycled heat:", self.brayton.recycledHeat * self.braytonFlowRate)
        print("Brayton external q_in:", (self.brayton.q_in() - self.brayton.recycledHeat) * self.braytonFlowRate)
        print()
        print("Rankine q_in:", self.rankine.q_in())
        print("Rankine recycled heat:", self.rankine.recycledHeat)
        print("Heat transeferred from brayton:", self.transferredHeat)
        print("Rankine external q_in:", self.rankine.q_in() - self.rankine.recycledHeat - self.transferredHeat)
        print()
        print("Total external q_in:", self.external_q_in())
        print()
        print()
        print("Brayton q_out:", self.brayton.q_out() * self.braytonFlowRate)
        print("Brayton recycled heat:", self.brayton.recycledHeat * self.braytonFlowRate)
        print("Heat transferred to rankine:", self.transferredHeat)
        print("Brayton external q_out:", (self.brayton.q_out() - self.brayton.recycledHeat) * self.braytonFlowRate - self.transferredHeat)
        print()
        print("Rankine q_out:", self.rankine.q_out())
        print("Rankine recycled heat:", self.rankine.recycledHeat)
        print("Rankine external q_out:", self.rankine.q_out() - self.rankine.recycledHeat)
        print()
        print("Total external q_out:", self.external_q_out())
        print()
        print()
        print("Brayton w_in:", self.brayton.w_in() * self.braytonFlowRate)
        print("Rankine w_in:", self.rankine.w_in())
        print("Total w_in;", self.w_in())
        print()
        print("Brayton w_out", self.brayton.w_out() * self.braytonFlowRate)
        print("Rankine w_out", self.rankine.w_out())
        print("Total w_out", self.w_out())
        print()
        print()
        print("Combined system properties:")
        print()
        print("Net Work:", self.w_net())
        print("Back-Work Ratio:", self.back_work())
        print("Efficiency:", self.efficiency())
        print("2nd Law Efficiency:", self.efficiency_2nd_law())
        print("Fuel Efficiency:", self.fuel_efficiency())
        print()
        print("---------------------------------------------------------------")
        print()
        
        
        

    def copy(self):
        return copy.deepcopy(self)
    


def graph_efficiency_vs_fuel_mass_frac(parentCycle, fuel, fuelCompressorEfficiency, minFrac, maxFrac, numPoints):
    fracs = []
    efficiencies = []
    for i in range(0, numPoints):
        frac = minFrac + (maxFrac - minFrac) * i / numPoints
        cycle = copy.deepcopy(parentCycle)
        cycle.inject_fuel_after_compressor(Mixture([["CH4", 1]]), frac, fuelCompressorEfficiency)
        cycle.combust()
        efficiencies.append(cycle.efficiency())
        fracs.append(frac)
    
    plt.figure()
    plt.xlabel("Fuel Mass fraction")
    plt.ylabel("Cycle Efficiency")
    plt.scatter(fracs, efficiencies)
    plt.show(block=False)

def graph_fuel_efficiency_vs_fuel_mass_frac(parentCycle, fuel, fuelCompressorEfficiency, minFrac, maxFrac, numPoints):
    fracs = []
    efficiencies = []
    for i in range(0, numPoints):
        frac = minFrac + (maxFrac - minFrac) * i / numPoints
        cycle = copy.deepcopy(parentCycle)
        cycle.inject_fuel_after_compressor(Mixture([["CH4", 1]]), frac, fuelCompressorEfficiency)
        cycle.combust()
        efficiencies.append(cycle.fuel_efficiency())
        fracs.append(frac)
    
    plt.figure()
    plt.xlabel("Fuel Mass fraction")
    plt.ylabel("Fuel Efficiency ($\\dot{w}_{out} / \\dot{m}_{fuel}$)")
    plt.scatter(fracs, efficiencies)
    plt.show(block=False)


def fuel_mass_fuel_efficiency_gradient(parentCycle, fuel, fuelCompressorEfficiency, frac):
    interval = 0.00001
    cycle1 = copy.deepcopy(parentCycle)
    cycle2 = copy.deepcopy(parentCycle)
    cycle1.inject_fuel_after_compressor(Mixture([["CH4", 1]]), frac - interval / 2, fuelCompressorEfficiency)
    cycle2.inject_fuel_after_compressor(Mixture([["CH4", 1]]), frac + interval / 2, fuelCompressorEfficiency)
    cycle1.combust()
    cycle2.combust()
    return (cycle2.fuel_efficiency() - cycle1.fuel_efficiency()) / interval

def optimize_fuel_mass_frac(parentCycle, fuel, fuelCompressorEfficiency, frac1, frac2, showGraph=False):
    if showGraph:
        plt.figure()
        plt.xlabel("iteration")
        plt.ylabel("Injected Fuel Mass Fraction")
        plt.title("Fuel Mass Fraction Convergence")
    
    tolerance = 0.0001
    i = 0
    while frac2 - frac1 > tolerance:
        frac = (frac1 + frac2) / 2
        grad = fuel_mass_fuel_efficiency_gradient(parentCycle, fuel, fuelCompressorEfficiency, frac)
        if grad > 0:
            frac1 = frac
        elif grad < 0:
            frac2 = frac
        else:
            return frac
        
        if showGraph:
            i += 1
            plt.plot([i, i], [frac1, frac2])
            plt.show(block=False)

    return (frac1 + frac2) / 2
    

def graph_efficiency_vs_feedwater_pressure(parentCycle, tempDiff, turbineIndex, numPoints, removeReheater=True):
    maxP = parentCycle.preTurbStates.P[0] - 10
    if turbineIndex > 0:
        maxP = parentCycle.postTurbStates.P[turbineIndex - 1] - 10
    
    minP = parentCycle.postTurbStates.P[turbineIndex + 1] + 10
    
    pressures = []
    efficiencies = []
    for i in range(0, numPoints):
        P = minP + (maxP - minP) * i / numPoints
        cycle = copy.deepcopy(parentCycle)
        cycle.postTurbStates.P[turbineIndex] = P
        cycle.set_closed_feedwater_heater_from_temp_diff(tempDiff, turbineIndex, removeReheater=removeReheater)
        cycle.solve_cp()

        if cycle.closedFeedwaterMassRatio < 1:
            pressures.append(P)
            efficiencies.append(cycle.efficiency())
    
    plt.figure()
    plt.xlabel("Feedwater Heater Pressure - $Pa$")
    plt.ylabel("Cycle Efficiency")
    plt.scatter(pressures, efficiencies)
    plt.show(block=False)

def feedwater_pressure_efficiency_gradient(parentCycle, P, tempDiff, turbineIndex, removeReheater=True):
    interval = 1
    cycle1 = copy.deepcopy(parentCycle)
    cycle2 = copy.deepcopy(parentCycle)
    cycle1.postTurbStates.P[turbineIndex] = P - interval / 2
    cycle2.postTurbStates.P[turbineIndex] = P + interval / 2
    cycle1.set_closed_feedwater_heater_from_temp_diff(tempDiff, turbineIndex, removeReheater=removeReheater)
    cycle2.set_closed_feedwater_heater_from_temp_diff(tempDiff, turbineIndex, removeReheater=removeReheater)
    cycle1.solve_cp()
    cycle2.solve_cp()
    
    if cycle2.efficiency() == 0:
        return -1 # stop trying to set feedwater above boiler temperature
    return (cycle2.efficiency() - cycle1.efficiency()) / interval

def optimize_feedwater_pressure(parentCycle, tempDiff, turbineIndex, removeReheater=True, showGraph=False):
    P2 = parentCycle.preTurbStates.P[0]
    if turbineIndex > 0:
        P2 = parentCycle.postTurbStates.P[turbineIndex - 1]
    P2 -= 10

    P1 = parentCycle.postTurbStates.P[turbineIndex + 1] + 10
    
    if showGraph:
        plt.figure()
        plt.xlabel("iteration")
        plt.ylabel("Feedwater Heater Pressure - $Pa$")
        plt.title("Feedwater Pressure Optimization Convergence")
    
    tolerance = 1
    i = 0
    while P2 - P1 > tolerance:
        P = (P1 + P2) / 2
        grad = feedwater_pressure_efficiency_gradient(parentCycle, P, tempDiff, turbineIndex, removeReheater=removeReheater)
        if grad > 0:
            P1 = P
        elif grad < 0:
            P2 = P
        else:
            return P
        
        if showGraph:
            i += 1
            plt.plot([i, i], [P1, P2])
            plt.show(block=False)
    
    return (P1 + P2) / 2


    

def optimize_turbine_exit_pressure_for_quality(parentCycle, turbineIndex, targetQuality, showGraph=False):
    if showGraph:
        plt.figure()
        plt.xlabel("iteration")
        plt.ylabel("Boiler Pressure - $Pa$")
        plt.title("Turbine " + str(turbineIndex) + " Pressure Optimization for Quality Convergence")
    
    P2 = parentCycle.preTurbStates.P[0]
    if turbineIndex > 0:
        P2 = parentCycle.postTurbStates.P[turbineIndex - 1]
    P2 -= 10

    i = 1
    P1 = None
    while P1 == None and i < len(parentCycle.postTurbStates.P):
        P1 = parentCycle.postTurbStates.P[turbineIndex + 1]
        i += 1
    if P1 == None:
        cycle = parentCycle.copy()
        cycle.solve_cp()
        cycle.propagate_assuming_no_pressure_losses()
        P1 = cycle.postTurbStates.P[-1]
    if P1 == None:
        P1 = 10
        print("ERROR No minimum pressure for cycle")

    P1 += 10
    
    tolerance = 1
    i = 0
    while P2 - P1 > tolerance:
        P = (P1 + P2) / 2
        cycle = copy.deepcopy(parentCycle)
        cycle.postTurbStates.P[turbineIndex] = P
        cycle.solve_cp()
        
        quality = cycle.postTurbStates.X[turbineIndex]
        if quality == -1: # superheated vapour
            quality = 1.1

        if quality == None:
            P1 = P
        
        elif quality < targetQuality:
            P1 = P
        elif quality > targetQuality:
            P2 = P
        else:
            return P
        
        if showGraph:
            i += 1
            plt.plot([i, i], [P1, P2])
            plt.show(block=False)

    return (P1 + P2) / 2



def graph_turbine_inlet_temp_vs_optimized_rankine_efficiency(parentCycle, turbineIndex, minT, maxT, turbineExitQualities, numPoints, closedFeedwaterTempDiff=10, closedFeedwaterIndex=None):
    # only pass an unsolved cycle to this, otherwise will create very broken cycles
    temperatures = []
    efficiencies = []
    for i in range(0, numPoints):
        T = minT + (maxT - minT) * i / numPoints
        cycle = parentCycle.copy()
        cycle.preTurbStates.T[turbineIndex] = T
        cycle.rankine_optimize_no_temp_changes(turbineExitQualities, closedFeedwaterTempDiff=closedFeedwaterTempDiff, closedFeedwaterIndex=closedFeedwaterIndex)
        temperatures.append(T)
        efficiencies.append(cycle.efficiency())

    plt.figure()
    plt.xlabel("Turbine " + str(turbineIndex) + " Inlet Temperature - $K$")
    plt.ylabel("Cycle Efficiency")
    plt.scatter(temperatures, efficiencies)
    plt.show(block=False)


def rankine_turbine_inlet_temp_efficiency_gradient(parentCycle, temp, turbineIndex, turbineExitQualities, feedwaterArgs=None):
    interval = 1
    cycle1 = copy.deepcopy(parentCycle)
    cycle2 = copy.deepcopy(parentCycle)
    cycle1.preTurbStates.T[turbineIndex] = temp - interval / 2
    cycle2.preTurbStates.T[turbineIndex] = temp + interval / 2
    cycle1.rankine_optimize_no_temp_changes(turbineExitQualities, feedwaterArgs=feedwaterArgs)
    cycle2.rankine_optimize_no_temp_changes(turbineExitQualities, feedwaterArgs=feedwaterArgs)

    X1 = cycle1.postTurbStates.X[-1]
    X2 = cycle2.postTurbStates.X[-1]

    if X1 == -1:
        X1 = 1.1
    if X2 == -1:
        X2 = 1.1

    # prevent decreasing temperature too much to a point it decreases quality past limits,
    # which would make other parts of the code kill efficiency:
    if X1 < turbineExitQualities or X2 < turbineExitQualities:
        return 1
    
    return (cycle2.efficiency() - cycle1.efficiency()) / interval

def rankine_optimize_turbine_inlet_temperature(parentCycle, turbineIndex, turbineExitQualities, feedwaterArgs=None, showGraph=False):
    # will only decrease from previous value, assuming increasing temperature is not an option
    if showGraph:
        plt.figure()
        plt.xlabel("iteration")
        plt.ylabel("Temperature - $Pa$")
        plt.title("Turbine " + str(turbineIndex) + " Inlet Temperature Optimization Convergence")
    
    solvedBase = parentCycle.copy()
    solvedBase.rankine_optimize_no_temp_changes(turbineExitQualities, feedwaterArgs=feedwaterArgs)
    T1 = max([solvedBase.postTurbStates.T[turbineIndex - 1], solvedBase.postTurbStates.T[turbineIndex]]) + 5
    T2 = solvedBase.preTurbStates.T[turbineIndex] - 5
    
    tolerance = 1
    i = 0
    while T2 - T1 > tolerance:
        T = (T1 + T2) / 2
        grad = rankine_turbine_inlet_temp_efficiency_gradient(parentCycle, T, turbineIndex, turbineExitQualities, feedwaterArgs=feedwaterArgs)
        if grad > 0:
            T1 = T
        elif grad < 0:
            T2 = T
        else:
            return T
        
        if showGraph:
            i += 1
            plt.plot([i, i], [T1, T2])
            plt.show(block=False)

    return T2 # use upper bound, not average to avoid risk of ending up just under required quality


def graph_turbine_exit_quality_vs_boiler_pressure(parentCycle, turbineExitQualities, numPoints, feedwaterArgs=None):
    P2 = cp.PropsSI(parentCycle.postCompStates.fluid[-1].to_coolprop(), "pcrit") - 20
    solvedParent = parentCycle.copy()
    solvedParent.solve_cp()
    solvedParent.propagate_assuming_no_pressure_losses()
    solvedParent.solve_cp()
    P1 = solvedParent.preCompStates.P[-1] + 10

    pressures = []
    qualities = []

    for i in range(0, numPoints):
        P = P1 + (P2 - P1) * i / numPoints
        cycle = parentCycle.copy()
        cycle.postCompStates.P[-1] = P
        cycle.propagate_assuming_no_pressure_losses()
        cycle.rankine_optimize_fixed_compression(turbineExitQualities, feedwaterArgs=feedwaterArgs)
        quality = cycle.postTurbStates.X[-1]
        
        if quality == -1:
            quality = 1

        pressures.append(P)
        qualities.append(quality)

    plt.figure()
    plt.scatter(pressures, qualities)
    plt.xlabel("Boiler Pressure - $Pa$")
    plt.ylabel("Turbine exit quality")
    plt.show(block=False)
    

def rankine_optimize_boiler_pressure(parentCycle, turbineExitQualities, feedwaterArgs=None, showGraph=False):
    if showGraph:
        plt.figure()
        plt.xlabel("iteration")
        plt.ylabel("Pressure - $Pa$")
        plt.title("Compressor Pressure Optimization Convergence")
    
    P2 = cp.PropsSI(parentCycle.postCompStates.fluid[-1].to_coolprop(), "pcrit") - 20
    solvedParent = parentCycle.copy()
    solvedParent.solve_cp()
    solvedParent.propagate_assuming_no_pressure_losses()
    solvedParent.solve_cp()
    P1 = solvedParent.preCompStates.P[-1] + 10
    
    if showGraph:
        plt.plot([-1, -1], [P1, P2])
    
    searchDelta = (P2 - P1) / 20
    alreadyDecreased = False
    
    while P2 > P1:
        cycle = parentCycle.copy()
        cycle.postCompStates.P[-1] = P2
        cycle.propagate_assuming_no_pressure_losses()
        cycle.rankine_optimize_fixed_compression(turbineExitQualities, feedwaterArgs=feedwaterArgs)
        quality = cycle.postTurbStates.X[-1]
        
        if quality == -1:
            quality = 1.1
        
        if quality >= turbineExitQualities:
            P1 = P2
            if alreadyDecreased:
                P2 += searchDelta
            break
        
        P2 -= searchDelta
        alreadyDecreased = True

    tolerance = 1
    i = 0
    while P2 - P1 > tolerance:
        P = (P1 + P2) / 2
        cycle = parentCycle.copy()
        cycle.postCompStates.P[-1] = P
        cycle.propagate_assuming_no_pressure_losses()
        cycle.rankine_optimize_fixed_compression(turbineExitQualities, feedwaterArgs=feedwaterArgs)
        quality = cycle.postTurbStates.X[-1]
        if quality == -1:
            quality = 1.1

        if quality > turbineExitQualities:
            P1 = P
        elif quality < turbineExitQualities:
            P2 = P
        else:
            return P

        
        if showGraph:
            i += 1
            plt.plot([i, i], [P1, P2])
            plt.show(block=False)

    return P1  # use lower bound, not average to avoid risk of ending up just under required quality


def rankine_fixed_quality_efficiency_worker(cycle, quality, feedwaterArgs):
    try:
        cycle.rankine_optimize_fixed_target_quality(quality, feedwaterArgs=feedwaterArgs)
        return cycle.efficiency()
    except:
        return 0


def graph_efficiency_vs_target_quality(parentCycle, minQuality, numPoints, feedwaterArgs=None):
    tasks = []
    qualities = []
    for i in range(0, numPoints):
        X = minQuality + (1 - minQuality) * i / numPoints
        tasks.append([parentCycle.copy(), X, feedwaterArgs])
        qualities.append(X)
    with multiprocessing.Pool(processes=3) as pool:
        efficiencies = pool.starmap(rankine_fixed_quality_efficiency_worker, tasks)

    plt.figure()
    plt.xlabel("Target Turbine Exit Quality")
    plt.ylabel("Efficiency")
    plt.scatter(qualities, efficiencies)
    plt.show(block=False)
    



# will have to look more into this, I don't think the newton-raphson is working
# might have local peaks or something:
def rankine_target_quality_efficiency_gradient(parentCycle, quality, feedwaterArgs=None):
    interval = 0.001
    cycle1 = parentCycle.copy()
    cycle2 = parentCycle.copy()

    if quality + interval/2 > 1:
        excess = 1 - quality - interval/2
        quality -= excess / 2
        interval -= excess

    tasks = [[cycle1, quality - interval / 2, feedwaterArgs], [cycle2, quality + interval / 2, feedwaterArgs]]
    
    with multiprocessing.Pool(processes=2) as pool:
        efficiencies = pool.starmap(rankine_fixed_quality_efficiency_worker, tasks)
    
    return (efficiencies[1] - efficiencies[0]) / 2
    

def rankine_optimize_target_quality(parentCycle, turbineExitQualities, feedwaterArgs=None):
    # starts with minimum quality for avoiding damaging blades, and tries to see if increasing quality will increase efficiency


    # checking minimum quality and saturated first, faster in most cases
    cycle1 = parentCycle.copy()
    cycle2 = parentCycle.copy()
    cycle3 = parentCycle.copy()
    cycle4 = parentCycle.copy()

    tasks = [[cycle1, turbineExitQualities, feedwaterArgs],
             [cycle2, turbineExitQualities + 0.001, feedwaterArgs],
             [cycle3, 0.999, feedwaterArgs],
             [cycle4, 1, feedwaterArgs]]

    with multiprocessing.Pool(processes=4) as pool:
        efficiencies = pool.starmap(rankine_fixed_quality_efficiency_worker, tasks)

    if efficiencies[0] == max(efficiencies):
        return turbineExitQualities
    if efficiencies[3] == max(efficiencies):
        return 1

    
    X1 = turbineExitQualities
    X2 = 1
    while X2 - X1 > 0.005:
        X = (X1 + X2) / 2
        grad = rankine_target_quality_efficiency_gradient(parentCycle, X, feedwaterArgs=feedwaterArgs)
        
        if grad > 0:
            X1 = X
        elif grad < 0:
            X2 = X
        else:
            return X

    return (X1 + X2) / 2



def init():
    multiprocessing.set_start_method('spawn')
    import warnings
    warnings.filterwarnings("ignore", message=".*Please use \"calc_sat_range\".*")
    warnings.filterwarnings("ignore", message=".*Detected an incomplete phase envelope, fixing it numerically.*")
