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

    rankine.rankine_optimize(turbineExitQualities, closedFeedwaterIndex=feedwaterIndex)

    rankine.print_design_info()
    rankine.print_energy_calcs()
    pd.graph_reheater_pressure_vs_efficiency(rankine, 0, 30)
    pd.plt.pause(5)
    rankineTS = rankine.ts_diagram(isRankine=True)
    rankineTS.set_title("Rankine Test Plot Ts")
    pd.plt.show(block=False)
