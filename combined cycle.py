if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    brayton = pd.load_cycle("cycles/basic brayton.ppdc")
    rankine = pd.load_cycle("cycles/500 600 rankine.ppdc")
    brayton.recycledHeat = 100000
    system = pd.CombinedBraytonRankineCycle(brayton, rankine)
    #system.adjust_brayton_flow_rate(10, showGraph=True)
    #system.braytonFlowRate = 4.590877532958984 # 600 600 rankine
    #system.braytonFlowRate = 4.166891098022461 # 500 600 rankine
    #system.braytonFlowRate = 4.868659973144531 # 500 600 rankine with 100000 brayton recuperation
    #system.add_brayton_cooling_states()

    print("BRAYTON CYCLE:")
    system.brayton.print_all_states()
    print()
    print()
    print("RANKINE CYCLE:")
    system.rankine.print_all_states()
    print()
    print()
    print("COMBINED SYSTEM:")
    system.print_energy_calcs()

    system.transfer_th_diagram()
    #pd.plt.pause(5)
    #system.ts_diagram()
