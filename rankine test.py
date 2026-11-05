if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    turbineExitQualities = 0.98
    feedwaterIndex = 1

    rankine = pd.Cycle(1, 3, pd.Mixture([["H2O", 1]]))
    rankine.preCompStates.T = [273.15 + 40]
    rankine.preCompStates.X = [0]
    rankine.preTurbStates.T = [570 + 273.15, 570 + 273.15, None]
    
    rankine.turbineEfficiencies = [0.85, 0.85, 0.85]
    rankine.compressorEfficiencies = [0.8]
    
    rankine.propagate_assuming_no_pressure_losses()

    #pd.graph_turbine_inlet_temp_vs_optimized_rankine_efficiency(rankine, 1, 370, 570 + 273.15, turbineExitQualities, 30, closedFeedwaterIndex=feedwaterIndex)
    #print(pd.rankine_optimize_turbine_inlet_temperature(rankine, 1, turbineExitQualities, closedFeedwaterIndex=feedwaterIndex, showGraph=True))

    #rankine.rankine_optimize_no_temp_changes(turbineExitQualities, closedFeedwaterIndex=feedwaterIndex)
    
    rankine.rankine_optimize(turbineExitQualities, closedFeedwaterIndex=feedwaterIndex)

    rankine.print_design_info()
    rankine.print_energy_calcs()
    pd.graph_reheater_pressure_vs_efficiency(rankine, 0, 30)
    
    pd.plt.pause(3)
    rankineTS = rankine.ts_diagram(isRankine=True)
    rankineTS.set_title("Rankine Test Plot Ts")
    pd.plt.show(block=False)
