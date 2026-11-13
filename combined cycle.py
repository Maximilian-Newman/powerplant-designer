if __name__ == "__main__":
    import powerplant_designer as pd
    pd.init()

    brayton = pd.load_cycle("brayton.ppdc")
    rankine = pd.load_cycle("nuclear rankine.ppdc")
    system = pd.CombinedBraytonRankineCycle(brayton, rankine)
    system.adjust_brayton_flow_rate(10, showGraph=True)

    system.ts_diagram()
