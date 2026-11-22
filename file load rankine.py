if __name__ == "__main__":
    import powerplant_designer as pd
    from tkinter import filedialog
    pd.init()

    name = filedialog.askopenfilename(initialdir = __file__, filetypes=(("Power-Plant Designer Cycle", "*.ppdc"),))
    rankine = pd.load_cycle(name)
    
    rankine.print_design_info()
    rankine.print_all_states()
    rankine.print_energy_calcs()
    rankineTS = rankine.ts_diagram(isRankine=True)
    rankineTS.set_title("Rankine From File Test Plot Ts")
    pd.plt.show(block=False)
