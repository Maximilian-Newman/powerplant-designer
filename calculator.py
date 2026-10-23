import cea
from CoolProp import CoolProp as cp # haha very funny
import matplotlib.pyplot as plt
import math

class SetOfStates():
    def __init__(self, numStates):
        self.T = [None] * numStates
        self.P = [None] * numStates
        self.v = [None] * numStates
        self.h = [None] * numStates
        self.s = [None] * numStates
        self.X = [None] * numStates
        self.fluid = ["Water"] * numStates

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

    def get_intermediate_points(self, a, b, numPoints = 20):
        #if self.v[a] != self.v[b]:
        #    return self.polytropic_interpolation(a, b, numPoints)

        
        # I strongly suspect that my linear interpolation of entropy approach is very inaccurate. to revisit later.
        # TS looks better with P interpolation, but still wrong
        T = []
        P = []
        v = []
        h = []
        s = []
        X = []
        for h_loc in range(0, numPoints):
            #s_loc = self.s[a] + h_loc * (self.s[b] - self.s[a]) / numPoints
            P_loc = self.P[a] + h_loc * (self.P[b] - self.P[a]) / numPoints
            h_loc = self.h[a] + h_loc * (self.h[b] - self.h[a]) / numPoints
            #T.append(cp.PropsSI("T", "Hmass", h_loc, "Smass", s_loc, self.fluid[a]))
            #P.append(cp.PropsSI("P", "Hmass", h_loc, "Smass", s_loc, self.fluid[a]))
            h.append(h_loc)
            #s.append(s_loc)
            #x.append(cp.PropsSI("Q", "Hmass", h_loc, "Smass", s_loc, self.fluid[a]))
            T.append(cp.PropsSI("T", "Hmass", h_loc, "P", P_loc, self.fluid[a]))
            P.append(P_loc)
            v.append(1 / cp.PropsSI("D", "Hmass", h_loc, "P", P_loc, self.fluid[a]))
            s.append(cp.PropsSI("Smass", "Hmass", h_loc, "P", P_loc, self.fluid[a]))
            X.append(cp.PropsSI("Q", "Hmass", h_loc, "P", P_loc, self.fluid[a]))

        return T, P, h, s, X

    def polytropic_interpolation(self, a, b, numPoints = 20):
        intermediates = SetOfStates(0)
        k = math.log(self.P[b] / self.P[a]) / math.log(self.v[a] / self.v[b])
        c = self.P[a] + self.v[a] ** k

        print()
        print(k, c)

        for i in range(0, numPoints):
            v = self.v[a] + (self.v[b] - self.v[a]) * i / numPoints
            P = c / v**k
            intermediates.T.append(None)
            intermadiates.P.append(P)
            intermediates.v.append(v)
            intermediates.h.append(None)
            intermediates.s.append(None)
            intermediates.X.append(None)
        intermediates.solve_points()
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





def get_intermediate_points(fromSet, fromI, toSet, toI, linearProperties, numPoints=20):
    intermediates = SetOfStates(0)
    starts = [fromSet.get_property(fromI, symbol) for symbol in linearProperties]
    ends = [toSet.get_property(toI, symbol) for symbol in linearProperties]
    for i in range(0, numPoints):
        linValues = []
        for j in range(0, len(linearProperties)):
            val = starts[j] + (ends[j] - starts[j]) * i / numPoints
            linValues.append(val)
            if linearProperties[i] == "T":
                intermediates.T.append(val)
        
        if "T" not in linearProperties:
            intermediates.T.append(cp.PropsSI("T", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
        intermediates.P.append(cp.PropsSI("P", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
        intermediates.v.append(1 / cp.PropsSI("D", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
        intermediates.h.append(cp.PropsSI("h", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
        intermediates.s.append(cp.PropsSI("s", linearProperties[0], linValues[0], linearProperties[1], linValues[1], fromSet.fluid[fromI]))
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
    #if toSet.h[i] != None:
    #    knowns.append(["Hmass", toSet.h[i]])
    #if toSet.s[i] != None:
    #    knowns.append(["Smass", toSet.s[i]])
    if toSet.X[i] != None:
        knowns.append(["Q", toSet.X[i]])


    if len(knowns) == 0:
        print("Couldn't solve point", i+1, "in solve_point_inefficient")
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
    #if toSet.h[i] == None:
    #    toSet.h[i] = cp.PropsSI("Hmass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])
    if toSet.s[i] == None:
        toSet.s[i] = cp.PropsSI("Smass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])
    if toSet.X[i] == None:
        toSet.X[i] = cp.PropsSI("Q", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])










class Cycle: # should be valid for both rankine and brayton
    def __init__(self, numComp=1, numTurb=1):
        self.preCompStates = SetOfStates(numComp)
        self.postCompStates = SetOfStates(numComp)
        self.preTurbStates = SetOfStates(numTurb)
        self.postTurbStates = SetOfStates(numTurb)

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
                solve_point_inefficient(self.postCompStates, i, self.preCompStates, i, self.compressorEfficiencies[i])
            elif self.postCompStates.point_is_solved(i):
                solve_point_inefficient(self.preCompStates, i, self.postCompStates, i, 1 / self.compressorEfficiencies[i])

        for i in range(0, len(self.preTurbStates.T)):
            if self.preTurbStates.point_is_solved(i) and self.postTurbStates.point_is_solved(i):
                h_isen = cp.PropsSI("Hmass", "P", self.postTurbStates.P[i], "Smass", self.preTurbStates.s[i])
                self.turbineEfficiencies[turbNum] = (self.postTurbStates.h[i] - self.preTurbStates.h[i]) / (h_isen - self.preTurbStates.h[i])
            elif self.preTurbStates.point_is_solved(i):
                solve_point_inefficient(self.postTurbStates, i, self.preTurbStates, i, self.turbineEfficiencies[i])
            elif self.postTurbStates.point_is_solved(i):
                solve_point_inefficient(self.preTurbStates, i, self.postTurbStates, i, 1 / self.turbineEfficiencies[i])


    def ts_diagram(self, ax=None, label = None):
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

        ax.scatter(x, y, label)

        x = []
        y = []

        #for i in range(0, len(allPoints.T) - 1):
        #    inter = allPoints.get_intermediate_points(i, i+1)
        #    x.extend(inter.s)
        #    y.extend(inter.T)
        #inter = allPoints.get_intermediate_points(-1, 0)
        #x.extend(inter.s)
        #y.extend(inter.T)

        for i in range(0, len(self.preCompStates.T)):
            inter = get_intermediate_points(self.preCompStates, i, self.postCompStates, i, ["h", "T"])
            x.extend(inter.s)
            y.extend(inter.T)

            if i+1 < len(self.preCompStates.T): # intercooling
                inter = get_intermediate_points(self.postCompStates, i, self.preCompStates, i+1, ["h", "T"])
                x.extend(inter.s)
                y.extend(inter.T)

        inter = get_intermediate_points(self.postCompStates, -1, self.preTurbStates, 1, ["h", "s"])
        x.extend(inter.s)
        y.extend(inter.T)
        
        for i in range(0, len(self.preTurbStates.T)):
            inter = get_intermediate_points(self.preTurbStates, i, self.postTurbStates, i, ["h", "T"])
            x.extend(inter.s)
            y.extend(inter.T)

            if i+1 < len(self.preTurbStates.T): # reheating
                inter = get_intermediate_points(self.postTurbStates, i, self.preTurbStates, i+1, ["h", "T"])
                x.extend(inter.s)
                y.extend(inter.T)

        ax.plot(x, y)
        ax.scatter(x, y, marker=".")
        
        if label != None:
            ax.legend()
        plt.show(block=False)




testRankine = Cycle(1, 1)
testRankine.preCompStates.T = [298]
testRankine.preCompStates.P = [8000]
testRankine.postCompStates.P = [8000000]
testRankine.preTurbStates.T = [600]
testRankine.preTurbStates.P = [8000000]
testRankine.postTurbStates.P = [8000]
testRankine.turbineEfficiencies[0] = 0.7
testRankine.compressorEfficiencies[0] = 0.7

testRankine.solve()

states = testRankine.combined_states()
print(states.T)
print(states.P)
print(states.v)
print(states.h)
print(states.s)
print(states.X)
testRankine.ts_diagram()



testRankine = Cycle(2, 2)
testRankine.preCompStates.T = [298, 298]
testRankine.preCompStates.P = [8000, 4000000]
testRankine.postCompStates.P = [4000000, 8000000]
testRankine.preTurbStates.T = [600, 600]
testRankine.preTurbStates.P = [8000000, 4000000]
testRankine.postTurbStates.P = [4000000, 8000]
testRankine.turbineEfficiencies = [0.7, 0.7]
testRankine.compressorEfficiencies = [0.7, 0.7]

testRankine.solve()

states = testRankine.combined_states()
print(states.T)
print(states.P)
print(states.h)
print(states.s)
print(states.X)
testRankine.ts_diagram()



