import cea
from CoolProp import CoolProp as cp # haha very funny
import matplotlib.pyplot as plt

class SetOfStates():
    def __init__(self, numStates):
        self.T = [None] * numStates
        self.P = [None] * numStates
        self.h = [None] * numStates
        self.s = [None] * numStates
        self.X = [None] * numStates
        self.fluid = ["Water"] * numStates

    def extend(self, otherState):
        self.T.extend(otherState.T)
        self.P.extend(otherState.P)
        self.h.extend(otherState.h)
        self.s.extend(otherState.s)
        self.X.extend(otherState.X)
        self.fluid.extend(otherState.fluid)

    def solve_points(self):
        for i in range(0, len(self.T)):
            knowns = []
            if self.T[i] != None:
                knowns.append(["T", self.T[i]])
            if self.P[i] != None:
                knowns.append(["P", self.P[i]])
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
        # I strongly suspect that my linear interpolation of entropy approach is very inaccurate. to revisit later.
        # TS looks better with P interpolation, but still wrong
        T = []
        P = []
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
            s.append(cp.PropsSI("Smass", "Hmass", h_loc, "P", P_loc, self.fluid[a]))
            X.append(cp.PropsSI("Q", "Hmass", h_loc, "P", P_loc, self.fluid[a]))

        return T, P, h, s, X







def solve_point_inefficient(toSet, i, fromSet, fromIndex, enthalpyMultiplier):
    knowns = []
    if toSet.T[i] != None:
        knowns.append(["T", toSet.T[i]])
    if toSet.P[i] != None:
        knowns.append(["P", toSet.P[i]])
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
    #if toSet.h[i] == None:
    #    toSet.h[i] = cp.PropsSI("Hmass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])
    if toSet.s[i] == None:
        toSet.s[i] = cp.PropsSI("Smass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])
    if toSet.X[i] == None:
        toSet.X[i] = cp.PropsSI("Q", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], toSet.fluid[i])










class Cycle: # should be valid for both rankine and brayton
    def __init__(self, numComp=1, numTurb=1):
        #self.T = [None] * (numComp + numTurb + 2)
        #self.P = [None] * (numComp + numTurb + 2)
        #self.h = [None] * (numComp + numTurb + 2)
        #self.s = [None] * (numComp + numTurb + 2)
        #self.x = [None] * (numComp + numTurb + 2)

        self.preCompStates = SetOfStates(numComp)
        self.postCompStates = SetOfStates(numComp)
        self.preTurbStates = SetOfStates(numTurb)
        self.postTurbStates = SetOfStates(numTurb)

        self.allStateSets = [self.preCompStates, self.postCompStates, self.preTurbStates, self.postTurbStates]
        
        self.turbineEfficiencies = [1] * numTurb
        self.compressorEfficiencies = [1] * numComp

    def combined_states(self):
        # READ ONLY
        states = SetOfStates(0)
        for subset in self.allStateSets:
            states.extend(subset)
        return states

    """
    def post_compressor_indices(self):
        return [i+1 for i in range(0, len(self.compressorEfficiencies))]

    def post_turbine_indices(self):
        return [i + len(self.compressorEfficiencies) + 2 for i in range(0, len(self.turbineEfficiencies))]


    def solve_points(self):
        for i in range(0, len(self.T)):
            knowns = []
            if self.T[i] != None:
                knowns.append(["T", self.T[i]])
            if self.P[i] != None:
                knowns.append(["P", self.P[i]])
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
            if self.h[i] == None:
                self.h[i] = cp.PropsSI("Hmass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])
            if self.s[i] == None:
                self.s[i] = cp.PropsSI("Smass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])
            if self.X[i] == None:
                self.X[i] = cp.PropsSI("Q", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])


    def solve_point_inefficient(self, i, fromIndex, enthalpyMultiplier):
        knowns = []
        if self.T[i] != None:
            knowns.append(["T", self.T[i]])
        if self.P[i] != None:
            knowns.append(["P", self.P[i]])
        #if self.h[i] != None:
        #    knowns.append(["Hmass", self.h[i]])
        #if self.s[i] != None:
        #    knowns.append(["Smass", self.s[i]])
        if self.X[i] != None:
            knowns.append(["Q", self.X[i]])


        if len(knowns) == 0:
            print("Couldn't solve point", i+1, "in solve_point_inefficient")
            return

        h_isen = cp.PropsSI("Hmass", knowns[0][0], knowns[0][1], "S", self.s[fromIndex], self.fluid[i])
        self.h[i] = self.h[fromIndex] + enthalpyMultiplier * (h_isen - self.h[fromIndex])
        
        knowns.append(["Hmass", self.h[i]])

        if self.T[i] == None:
            self.T[i] = cp.PropsSI("T", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])
        if self.P[i] == None:
            self.P[i] = cp.PropsSI("P", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])
        #if self.h[i] == None:
        #    self.h[i] = cp.PropsSI("Hmass", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])
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


    def solve(self):
        self.solve_points()

        compNum = 0
        for i in self.post_compressor_indices():
            if self.point_is_solved(i-1) and self.point_is_solved(i):
                h_isen = cp.PropsSI("Hmass", "P", self.P[i], "Smass", self.s[i-1])
                self.compressorEfficiencies[compNum] = (self.h[i] - self.h[i-1]) / (h_isen - self.h[i-1])
            elif self.point_is_solved(i-1):
                self.solve_point_inefficient(i, i-1, self.compressorEfficiencies[compNum])
            elif self.point_is_solved(i):
                self.solve_point_inefficient(i-1, i, 1 / self.compressorEfficiencies[compNum])
            compNum += 1

        turbNum = 0
        for i in self.post_turbine_indices():
            if self.point_is_solved(i-1) and self.point_is_solved(i):
                h_isen = cp.PropsSI("Hmass", "P", self.P[i], "Smass", self.s[i-1])
                self.turbineEfficiencies[turbNum] = (self.h[i] - self.h[i-1]) / (h_isen - self.h[i-1])
            elif self.point_is_solved(i-1):
                self.solve_point_inefficient(i, i-1, self.turbineEfficiencies[turbNum])
            elif self.point_is_solved(i):
                self.solve_point_inefficient(i-1, i, 1 / self.turbineEfficiencies[turbNum])
            turbNum += 1
    """

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

    """
    def get_intermediate_points(self, a, b, numPoints = 20):
        # I strongly suspect that my linear interpolation of entropy approach is very inaccurate. to revisit later.
        # TS looks better with P interpolation, but still wrong
        T = []
        P = []
        h = []
        s = []
        x = []
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
            s.append(cp.PropsSI("Smass", "Hmass", h_loc, "P", P_loc, self.fluid[a]))
            x.append(cp.PropsSI("Q", "Hmass", h_loc, "P", P_loc, self.fluid[a]))

        return T, P, h, s, x
        
    
    def ts_diagram(self, ax=None, label = None):
        if ax == None:
            ax = plt.figure().add_subplot()
            ax.set_xlabel("s - $J kg^{-1} K^{-1}$")
            ax.set_ylabel("T - $K$")

        x = []
        y = []
        for i in range(0, len(self.T)):
            x.append(self.s[i])
            y.append(self.T[i])

        ax.scatter(x, y, label)

        x = []
        y = []

        for i in range(0, len(self.T) - 1):
            T, P, h, s, qual = self.get_intermediate_points(i, i+1)
            x.extend(s)
            y.extend(T)
        T, P, h, s, qual = self.get_intermediate_points(-1, 0)
        x.extend(s)
        y.extend(T)

        ax.plot(x, y)
        ax.scatter(x, y, marker=".")
        
        if label != None:
            ax.legend()
        plt.show(block=False)
    """

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

        for i in range(0, len(allPoints.T) - 1):
            T, P, h, s, qual = allPoints.get_intermediate_points(i, i+1)
            x.extend(s)
            y.extend(T)
        T, P, h, s, qual = allPoints.get_intermediate_points(-1, 0)
        x.extend(s)
        y.extend(T)

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
print(states.h)
print(states.s)
print(states.X)
testRankine.ts_diagram()



"""
testRankine = Cycle()
testRankine.T = [298, None, 600, None]
testRankine.P = [8000, 8000000, 8000000, 8000]
testRankine.turbineEfficiencies[0] = 0.7
testRankine.compressorEfficiencies[0] = 0.7

testRankine.solve()

print(testRankine.T)
print(testRankine.P)
print(testRankine.h)
print(testRankine.s)
print(testRankine.x)

testRankine.ts_diagram()




testRankine = Cycle(2,2)
testRankine.T = [298, None, None, 600, None, None]
testRankine.P = [8000, 4000000, 8000000, 8000000, 4000000, 8000]
#testRankine.X[4] = 1
testRankine.turbineEfficiencies[0] = 0.7
testRankine.compressorEfficiencies[0] = 0.7

testRankine.solve()

print(testRankine.T)
print(testRankine.P)
print(testRankine.h)
print(testRankine.s)
print(testRankine.x)

testRankine.ts_diagram()
"""
