import numpy as np
import math

coords = np.empty((0))
with open("unitcell", "r") as file:
    for line in file.readlines()[2:]:
        line_s = line.split()
        if line_s[0] != "HETATM":
            break
        #print(line_s[5:8])line_s[5:8]
        coords = np.append(coords, np.array(line_s[5:8],dtype=float))

coords = coords.reshape(-1,3)

def CalcUnitCell(coords):
    def ShiftOrigin(coords):
        coordsT = coords.T
        x_min, y_min, z_min = coordsT[0].min() , coordsT[1].min(), coordsT[2].min()
        coordsT[0], coordsT[1], coordsT[2] = coordsT[0] + abs(x_min), coordsT[1] + abs(y_min), coordsT[2] + abs(z_min)
        return coordsT.T
    
    shiftedCoords = ShiftOrigin(coords)
    
    x_max, y_max, z_max = shiftedCoords.T[0].max() , shiftedCoords.T[1].max(), shiftedCoords.T[2].max()
    return shiftedCoords, (x_max, y_max, z_max)

coords , cell  = CalcUnitCell(coords)
x_offset = 15
atomsWritten = 0
with open("unitcell", "r") as file:
    with open("NewCell.pdb", "w") as writer:
        for i,line in enumerate(file.readlines()):
            if i == 2:
                writer.write(f"CRYST1   {math.ceil(cell[0])+x_offset}.000   {math.ceil(cell[1])}.000   {math.ceil(cell[2])}.000  90.00  90.00  90.00 P 1           1\n")
            if line.startswith(("HETATM", "ATM")):
                lineS = line.split()  
                writer.write(f"{lineS[0]}   {lineS[1]:>2}  {lineS[2]}   {lineS[3]}     {lineS[4]}      {coords[atomsWritten][0]:.3f}  {coords[atomsWritten][1]:.3f}  {coords[atomsWritten][2]:.3f}  {lineS[8]}  {lineS[9]}           {lineS[10]}\n")  
                atomsWritten += 1
            else:
                writer.write(line)