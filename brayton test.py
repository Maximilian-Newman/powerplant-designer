if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    brayton = pd.Cycle(1, 1, pd.AIR())
    brayton.preCompStates.P = [101325]
    brayton.postCompStates.P = [500000]
    brayton.preCompStates.T = [298]

    brayton.turbineEfficiencies = [0.87]
    brayton.compressorEfficiencies = [0.85]
    
    brayton.propagate_assuming_no_pressure_losses()
    
    #brayton.add_pre_combustion_heat_exchange(100000)

    frac = pd.optimize_fuel_mass_frac(brayton, pd.Mixture([["CH4", 1]]), 0.85, 0.005, 0.08, showGraph=True)
    
    #pd.graph_fuel_efficiency_vs_fuel_mass_frac(brayton, pd.Mixture([["CH4", 1]]), 0.85, 0.005, 0.08, 20)
    #pd.graph_efficiency_vs_fuel_mass_frac(brayton, pd.Mixture([["CH4", 1]]), 0.85, 0.005, 0.08, 20)

    print("Optimal fuel mass fraction:", frac)
    brayton.inject_fuel_after_compressor(pd.Mixture([["CH4", 1]]), frac, 0.85)
    brayton.combust()

    brayton.save("brayton.ppdc")
    brayton.print_design_info()
    brayton.print_all_states()
    brayton.print_energy_calcs()
    brayton.mixture_pie_charts("Test Brayton Cycle")
    brayton.ts_diagram()
    brayton.pv_diagram()
