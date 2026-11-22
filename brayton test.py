if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()
    
    brayton = pd.Cycle(2, 1, pd.AIR())
    brayton.preCompStates.P = [101325, None]
    brayton.postCompStates.P = [None, 500000]
    brayton.preCompStates.T = [298, 310]
    brayton.propagate_assuming_no_pressure_losses()

    brayton.turbineEfficiencies = [0.87]
    brayton.compressorEfficiencies = [0.85, 0.85]
    
    #pd.graph_fuel_efficiency_vs_intercooler_pressure(brayton, 0, pd.Mixture([["CH4", 1]]), 0.85, 0.08, regeneratorDeltaT=10)
    #import sys
    #sys.exit()
    
    brayton.brayton_optimize(pd.Mixture([["CH4", 1]]), 0.85, 0.05, regeneratorDeltaT=10)
    
    brayton.save("cycles/regenerated brayton with intercooler.ppdc")
    brayton.print_design_info()
    brayton.print_all_states()
    brayton.print_energy_calcs()
    
    brayton.mixture_pie_charts("Brayton Cycle")
    brayton.ts_diagram()
    brayton.pv_diagram()
