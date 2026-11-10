if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    rankine = pd.load_cycle("rankine.txt")
    
    rankine.print_design_info()
    rankine.print_energy_calcs()
    rankineTS = rankine.ts_diagram(isRankine=True)
    rankineTS.set_title("Rankine From File Test Plot Ts")
    pd.plt.show(block=False)
