def writePlumed(dihedral, METAD=False, RATES=False,RESTRAIN=False, PRINT=True, temp=None):
    if temp == None:
        temp = 310.0

    out = "UNITS LENGTH=A TIME=0.001\n"
    out += f"t: TORSION ATOMS={dihedral[0]},{dihedral[1]},{dihedral[2]},{dihedral[3]}\n"
    out += f"a: ALPHABETA ATOMS1={dihedral[0]},{dihedral[1]},{dihedral[2]},{dihedral[3]} REFERENCE=3.14\n\n"

    if RESTRAIN:
        out += "RESTRAINT ARG=a AT=0 KAPPA=2000.0 LABEL=restraint\n\n"

    if RATES:
        out += "COMMITTOR ...\n"
        lines = [
        "ARG=a",
        "STRIDE=10 ",
        "BASIN_LL1=0.9",
        "BASIN_UL1=1.0"
        ]
        for word in lines:
           out += f"\t{word}\n"
        out += "...\n\n"

    if METAD:
        pace = 200
        if RATES:
            pace = 500000

        lines = [
        "LABEL=metad",
        "ARG=t",
        f"PACE={pace}",
        "HEIGHT=1.0",
        "SIGMA=0.1",
        "GRID_MIN=-pi",
        "GRID_MAX=pi",
        "GRID_BIN=100",
        "CALC_RCT",
        "FILE=HILLS",
        "BIASFACTOR=60",
        f"TEMP={temp}"
        ]
        out += "METAD ...\n"
        for word in lines:
           out += f"\t{word}\n"

        if RATES:
            out += f"\tACCELERATION\n"
        
        out += "... METAD\n\n"

    if PRINT:
        out += "PRINT FILE=COLVAR ARG=t,a,metad.bias,metad.rbias,metad.rct STRIDE=10\n"
        if RATES:
            out += "PRINT FILE=ACCEL ARG=metad.acc\n"

        out += "FLUSH STRIDE=10\n"
        
    
    return out



if __name__ == "__main__":
    with open("plumedTest", "w") as file:
        file.write(writePlumed([1,2,3,4],RESTRAIN=True, RATES=True, METAD=True))