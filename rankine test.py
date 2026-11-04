if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    turbineExitQualities = 0.95

    rankine = pd.Cycle(1, 3, pd.Mixture([["H2O", 1]]))
    rankine.preCompStates.T = [273.15 + 40]
    rankine.preCompStates.X = [0]
    rankine.preTurbStates.T = [570, 570, None]
    
    rankine.turbineEfficiencies = [0.85, 0.85, 0.85]
    rankine.compressorEfficiencies = [0.8]
    
    noFeedwater = rankine.copy()
    noFeedwater.remove_reheat_step(1)
    
    pd.graph_turbine_quality_vs_compressor_pressure(noFeedwater, 100000, 8000000, 0, -1, 30, reheaterQuality=turbineExitQualities)
    
    boilerPressure = pd.optimize_compressor_pressure_for_turbine_quality(noFeedwater, 100000, 8000000, 0, -1, turbineExitQualities, showGraph=True)
    rankine.postCompStates.P = [boilerPressure]
    rankine.propagate_assuming_no_pressure_losses()

    rankine.postTurbStates.P[0] = pd.optimize_turbine_exit_pressure_for_quality(rankine, 0, turbineExitQualities, showGraph=True)
    rankine.solve_cp()
    rankine.propagate_assuming_no_pressure_losses()
    rankine.solve_cp()

    pd.graph_efficiency_vs_feedwater_pressure(rankine, 10, 1, 200)

    P = pd.optimize_feedwater_pressure(rankine, 10, 1, showGraph=True)
    rankine.postTurbStates.P[1] = P
    rankine.set_closed_feedwater_heater_from_temp_diff(10, 1)

    rankine.solve_cp()
    rankine.print_energy_calcs()
    pd.plt.pause(1)
    rankineTS = rankine.ts_diagram(isRankine=True)
    rankineTS.set_title("Rankine Test Plot Ts")
    pd.plt.show()
