if __name__ == "__main__":
    import powerplant_designer as pd
    from tkinter import filedialog
    pd.init()

    name = filedialog.askopenfilename(initialdir = __file__, filetypes=(("Power-Plant Designer Cycle", "*.ppdc"),))
    brayton = pd.load_cycle(name)
    
    brayton.print_design_info()
    brayton.print_all_states()
    brayton.print_energy_calcs()
    brayton.mixture_pie_charts()
    TS = brayton.ts_diagram()
    TS.set_title("Brayton From File Ts Diagram")
    PV = brayton.pv_diagram()
    PV.set_title("Brayton From File Pv Diagram")
    pd.plt.show(block=False)
