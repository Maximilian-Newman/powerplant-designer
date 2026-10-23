import cea
from CoolProp import CoolProp as cp # haha very funny
from CoolProp.Plots import PropertyPlot
import CoolProp
import matplotlib.pyplot as plt
import math

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
                self.T[i] = cp.PropsSI("T", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])
            if self.P[i] == None:
                self.P[i] = cp.PropsSI("P", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])
            if self.v[i] == None:
                self.v[i] = 1 / cp.PropsSI("D", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])
            if self.h[i] == None:
                self.h[i] = cp.PropsSI("Hmass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])
            if self.s[i] == None:
                self.s[i] = cp.PropsSI("Smass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])
            if self.X[i] == None:
                self.X[i] = cp.PropsSI("Q", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])

    
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
        
        intermediates.T.append(cp.PropsSI("T", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
        intermediates.P.append(cp.PropsSI("P", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
        intermediates.v.append(1 / cp.PropsSI("D", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
        intermediates.h.append(cp.PropsSI("Hmass", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
        intermediates.s.append(cp.PropsSI("Smass", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
        intermediates.X.append(cp.PropsSI("Q", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
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

    h_isen = cp.PropsSI("Hmass", knowns[0][0], knowns[0][1], "S", fromSet.s[fromIndex], toSet.fluid[i])
    toSet.h[i] = fromSet.h[fromIndex] + enthalpyMultiplier * (h_isen - fromSet.h[fromIndex])
    
    knowns.append(["Hmass", toSet.h[i]])

    if toSet.T[i] == None:
        toSet.T[i] = cp.PropsSI("T", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])
    if toSet.P[i] == None:
        toSet.P[i] = cp.PropsSI("P", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])
    if toSet.v[i] == None:
        toSet.v[i] = 1 / cp.PropsSI("D", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])
    if toSet.s[i] == None:
        toSet.s[i] = cp.PropsSI("Smass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])
    if toSet.X[i] == None:
        toSet.X[i] = cp.PropsSI("Q", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])










class Cycle: # should be valid for both rankine and brayton
    def __init__(self, numComp=1, numTurb=1, fluid = "Water"):
        self.preCompStates = SetOfStates(numComp, fluid)
        self.postCompStates = SetOfStates(numComp, fluid)
        self.preTurbStates = SetOfStates(numTurb, fluid)
        self.postTurbStates = SetOfStates(numTurb, fluid)

        self.allStateSets = [self.preCompStates, self.postCompStates, self.preTurbStates, self.postTurbStates]
        
        self.turbineEfficiencies = [1] * numTurb
        self.compressorEfficiencies = [1] * numComp

    def combined_states(self):
        states = SetOfStates(0)
        for i in range(0, len(self.preCompStates.T)):
            states.extend(self.preCompStates.get_point(i))
            states.extend(self.postCompStates.get_point(i))
            
        for i in range(0, len(self.preTurbStates.T)):
            states.extend(self.preTurbStates.get_point(i))
            states.extend(self.postTurbStates.get_point(i))
        return states # returned object is READ ONLY


    def solve(self):
        for subset in self.allStateSets:
            subset.solve_points()

        for i in range(0, len(self.preCompStates.T)):
            if self.preCompStates.point_is_solved(i) and self.postCompStates.point_is_solved(i):
                h_isen = cp.PropsSI("Hmass", "P", self.postCompStates.P[i], "Smass", self.preCompStates.s[i])
                self.compressorEfficiencies[turbNum] = (self.postCompStates.h[i] - self.preCompStates.h[i]) / (h_isen - self.preCompStates.h[i])
            elif self.preCompStates.point_is_solved(i):
                solve_point_inefficient(self.postCompStates, i, self.preCompStates, i, 1 / self.compressorEfficiencies[i])
            elif self.postCompStates.point_is_solved(i):
                solve_point_inefficient(self.preCompStates, i, self.postCompStates, i, self.compressorEfficiencies[i])

        for i in range(0, len(self.preTurbStates.T)):
            if self.preTurbStates.point_is_solved(i) and self.postTurbStates.point_is_solved(i):
                h_isen = cp.PropsSI("Hmass", "P", self.postTurbStates.P[i], "Smass", self.preTurbStates.s[i])
                self.turbineEfficiencies[turbNum] = (self.postTurbStates.h[i] - self.preTurbStates.h[i]) / (h_isen - self.preTurbStates.h[i])
            elif self.preTurbStates.point_is_solved(i):
                solve_point_inefficient(self.postTurbStates, i, self.preTurbStates, i, self.turbineEfficiencies[i])
            elif self.postTurbStates.point_is_solved(i):
                solve_point_inefficient(self.preTurbStates, i, self.postTurbStates, i, 1 / self.turbineEfficiencies[i])


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

        x = []
        y = []
        

        for i in range(0, len(self.preCompStates.T)): # compressor
            inter = get_intermediate_points(self.preCompStates, i, self.postCompStates, i, ["h", "s"])
            x.extend(inter.s)
            y.extend(inter.T)

            if i+1 < len(self.preCompStates.T): # intercooling
                inter = get_intermediate_points(self.postCompStates, i, self.preCompStates, i+1, ["T", "P"])
                x.extend(inter.s)
                y.extend(inter.T)

        # boiler / combustion chamber
        inter = get_intermediate_points(self.postCompStates, -1, self.preTurbStates, 0, ["s", "P"])
        x.extend(inter.s)
        y.extend(inter.T)
        
        for i in range(0, len(self.preTurbStates.T)): # turbine
            inter = get_intermediate_points(self.preTurbStates, i, self.postTurbStates, i, ["h", "s"])
            x.extend(inter.s)
            y.extend(inter.T)

            if i+1 < len(self.preTurbStates.T): # reheating
                inter = get_intermediate_points(self.postTurbStates, i, self.preTurbStates, i+1, ["T", "P"])
                x.extend(inter.s)
                y.extend(inter.T)

        # condenser
        inter = get_intermediate_points(self.postTurbStates, -1, self.preCompStates, 0, ["s", "P"])
        x.extend(inter.s)
        y.extend(inter.T)
        

        ax.plot(x, y)
        #ax.scatter(x, y, marker=".") # for debug to ensure no large jumps

        if isRankine:
            original_xlim = ax.get_xlim()
            original_ylim = ax.get_ylim()
            propertyPlot = PropertyPlot(self.postCompStates.fluid[0], 'TS', axis=ax, unit_system='SI')
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




testRankine = Cycle(1, 1)
testRankine.preCompStates.P = [8000]
testRankine.preCompStates.X = [0]
testRankine.postCompStates.P = [8000000]
testRankine.preTurbStates.T = [600]
testRankine.preTurbStates.P = [8000000]
testRankine.postTurbStates.P = [8000]
testRankine.turbineEfficiencies[0] = 1
testRankine.compressorEfficiencies[0] = 0.7

testRankine.solve()

states = testRankine.combined_states()
rankineTS = testRankine.ts_diagram(label="unmodified", isRankine=True)



testRankine = Cycle(2, 2)
testRankine.preCompStates.T = [None, 400]
testRankine.preCompStates.P = [8000, 4000000]
testRankine.preCompStates.X = [0, None]
testRankine.postCompStates.P = [4000000, 8000000]
testRankine.preTurbStates.T = [600, 600]
testRankine.preTurbStates.P = [8000000, 4000000]
testRankine.postTurbStates.P = [4000000, 8000]
testRankine.turbineEfficiencies = [0.7, 0.7]

testRankine.solve()

states = testRankine.combined_states()
testRankine.ts_diagram(ax=rankineTS, label = "reheater and useless intercooler")


testBrayton = Cycle(1, 1, "Air")
testBrayton.preCompStates.P = [101325]
testBrayton.postCompStates.P = [500000]
testBrayton.preCompStates.T = [298]
testBrayton.preTurbStates.T = [600]
testBrayton.preTurbStates.P = [500000]
testBrayton.postTurbStates.P = [101325]

testBrayton.solve()
braytonTS = testBrayton.ts_diagram(label="isentropic")


testBrayton = Cycle(1, 1, "Air")
testBrayton.preCompStates.P = [101325]
testBrayton.postCompStates.P = [500000]
testBrayton.preCompStates.T = [298]
testBrayton.preTurbStates.T = [600]
testBrayton.preTurbStates.P = [500000]
testBrayton.postTurbStates.P = [101325]
testBrayton.turbineEfficiencies = [0.7]
testBrayton.compressorEfficiencies = [0.7]

testBrayton.solve()
testBrayton.ts_diagram(ax=braytonTS, label="70% efficient turbine&compressor")





