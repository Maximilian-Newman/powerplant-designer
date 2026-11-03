if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    rankine = pd.Cycle(1, 3, pd.Mixture([["H2O", 1]]))
    rankine.preCompStates.P = [8000]
    rankine.preCompStates.X = [0]
    rankine.preTurbStates.T = [570, 570, None]
    rankine.postTurbStates.P = [4000000, None, None]
    
    rankine.turbineEfficiencies = [0.85, 0.85, 0.85]
    rankine.compressorEfficiencies = [0.8]

    pd.graph_efficiency_vs_compressor_pressure(rankine, 100000, 8000000, 0, 50, feedwaterArgs=[10, 1, True])
    pd.plt.pause(1)
    
    boilerPressure = pd.optimize_compressor_pressure(rankine, 100000, 8000000, 0, feedwaterArgs=[10, 1, True], showGraph=True)
    rankine.postCompStates.P = [boilerPressure]
    rankine.propagate_assuming_no_pressure_losses()

    pd.graph_efficiency_vs_feedwater_pressure(rankine, 10, 1, 200)
    pd.plt.pause(1)

    P = pd.optimize_feedwater_pressure(rankine, 10, 1, showGraph=True)
    rankine.postTurbStates.P[1] = P
    pd.plt.pause(1)
    rankine.set_closed_feedwater_heater_from_temp_diff(10, 1)

    rankine.solve_cp()
    rankine.print_energy_calcs()
    rankineTS = rankine.ts_diagram(isRankine=True)
    rankineTS.set_title("Rankine Test Plot Ts")
