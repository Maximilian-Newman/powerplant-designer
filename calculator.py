import cea
from CoolProp import CoolProp as cp # haha very funny
import matplotlib.pyplot as plt

class Cycle: # should be valid for both rankine and brayton
    def __init__(self):
        self.T = [None, None, None, None]
        self.P = [None, None, None, None]
        self.h = [None, None, None, None]
        self.s = [None, None, None, None]
        self.x = [None, None, None, None]
        self.turbineEfficiencies = [1]
        self.compressorEfficiencies = [1]
        self.fluid = ["Water", "Water", "Water", "Water"]

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
            if self.x[i] != None:
                knowns.append(["Q", self.x[i]])

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
            if self.x[i] == None:
                self.x[i] = cp.PropsSI("Q", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])

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
        if self.x[i] != None:
            knowns.append(["Q", self.x[i]])


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
        if self.x[i] == None:
            self.x[i] = cp.PropsSI("Q", knowns[0][0], knowns[0][1], knowns[1][0], knowns[1][1], self.fluid[i])

    def point_is_solved(self, i):
        if self.T[i] == None: return False
        if self.P[i] == None: return False
        if self.h[i] == None: return False
        if self.s[i] == None: return False
        if self.x[i] == None: return False
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
