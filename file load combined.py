if __name__ == "__main__":
    import powerplant_designer as pd
    from tkinter import filedialog
    pd.init()

    path = filedialog.askopenfilename(initialdir = __file__, filetypes=(("Power-Plant Designer Cycle", "*.ppds"),))
    system = pd.load_combined(path)
    
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
    system.brayton.mixture_pie_charts()
    pd.plt.pause(5)
    system.ts_diagram(True)
    system.ts_diagram()
