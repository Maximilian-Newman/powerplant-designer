if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    turbineExitQualities = 0.9
    feedwaterIndex = 1

    rankine = pd.Cycle(1, 3, pd.Mixture([["H2O", 1]]))
    rankine.preCompStates.T = [273.15 + 40]
    rankine.preCompStates.X = [0]
    rankine.preTurbStates.T = [570 + 273.15, 570 + 273.15, None]
    #rankine.preTurbStates.T = [570, 570, None]
    
    rankine.turbineEfficiencies = [0.87, 0.87, 0.87]
    rankine.compressorEfficiencies = [0.85]
    
    rankine.propagate_assuming_no_pressure_losses()

    #pd.graph_efficiency_vs_target_quality(rankine, turbineExitQualities-0.05, 10, feedwaterArgs=[10, feedwaterIndex, True])
    #pd.plt.pause(5)
    
    rankine.rankine_optimize(turbineExitQualities, feedwaterArgs=[10, feedwaterIndex, True])
    
    rankine.print_design_info()
    rankine.print_energy_calcs()
    rankine.save("rankine.txt")
    rankineTS = rankine.ts_diagram(isRankine=True)
    rankineTS.set_title("Rankine Test Plot Ts")
    pd.plt.show(block=False)
