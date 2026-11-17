if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    brayton = pd.load_cycle("cycles/basic brayton.ppdc")
    rankine = pd.load_cycle("cycles/nuclear rankine.ppdc")
    system = pd.CombinedBraytonRankineCycle(brayton, rankine)
    #system.adjust_brayton_flow_rate(20)
    system.braytonFlowRate = 4.346672058105469
    system.add_brayton_cooling_states()

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
    pd.plt.pause(5)
    system.ts_diagram()

