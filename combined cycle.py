if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    brayton = pd.load_cycle("basic brayton.ppdc")
    rankine = pd.load_cycle("nuclear rankine.ppdc")
    system = pd.CombinedBraytonRankineCycle(brayton, rankine)
    system.adjust_brayton_flow_rate(10, showGraph=True)

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
    system.ts_diagram()
