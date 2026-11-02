if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    rankine = pd.Cycle(1, 3, pd.Mixture([["H2O", 1]]))
    rankine.preCompStates.P = [8000]
    rankine.preCompStates.X = [0]
    rankine.postCompStates.P = [7000000]
    rankine.preTurbStates.T = [600, 600, None]
    rankine.preTurbStates.P = [7000000, None, None]
    rankine.postTurbStates.P = [4000000, None, None]
    
    rankine.turbineEfficiencies = [0.85, 0.85, 0.85]
    rankine.compressorEfficiencies = [0.8]
    rankine.propagate_assuming_no_pressure_losses()

    pd.graph_efficiency_vs_feedwater_pressure(rankine, 10000, 2500000, 10, 1, 30)
    pd.plt.pause(1)

    rankine.postTurbStates.P[1] = 200000
    rankine.set_closed_feedwater_heater_from_temp_diff(10, 1)

    rankine.solve_cp()
    rankine.print_energy_calcs()
    rankineTS = rankine.ts_diagram(isRankine=True)
    rankineTS.set_title("Rankine Test Plot Ts")
