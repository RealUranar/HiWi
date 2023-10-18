import numpy as np
import math


def makeUnitcell(inFile :str, z_ySideLengh = None):
    """Encases a molecule inside a Unitcell

    Args:
        inName (str, optional): A Molecule in PDB format as input. Defaults to "NEWPDB.PDB".
        outName (str, optional): A Molecule in PDB format as output. Defaults to "SHIFTED.PDB".
        cornerAtom (str, optional): Places the selected atom at the 0,0,0 position in the unitcell, If None is passed, places the molecule just inside the cell. Defaults to "S". 
        z_ySideLengh (int, optional): Side lenght of the unitcell. Defaults to 6.

    Returns:
        _type_: Nothing
    """
    coords = np.empty((0))
    for line in inFile.split("\n"):
        if not line.startswith(("HETATM", "ATOM")):
            continue
        coords = np.append(coords, np.array(line[32:55].split(),dtype=float))

    coords = coords.reshape(-1,3)

    def CalcUnitCell(coords):
        def ShiftOrigin(coords):
            coordsT = coords.T
            x_min, y_min, z_min = coordsT[0].min() , coordsT[1].min(), coordsT[2].min()

            coordsT[0], coordsT[1], coordsT[2] = coordsT[0] + abs(x_min), coordsT[1] + abs(y_min), coordsT[2] + abs(z_min)
            return coordsT.T
        
        shiftedCoords = ShiftOrigin(coords)
        
        x_max, y_max, z_max = shiftedCoords.T[0].max() , shiftedCoords.T[1].max(), shiftedCoords.T[2].max()
        return shiftedCoords, [x_max, y_max, z_max]

    coords , cell  = CalcUnitCell(coords)
    x_offset = 15
    if z_ySideLengh != None:
        cell[1], cell[2] = z_ySideLengh, z_ySideLengh  #Make the cell 6x6 in y and z direction
    
    atomIndex = 0
    outFile = f"CRYST1   {math.ceil(cell[0])+x_offset}.000   {math.ceil(cell[1])}.000   {math.ceil(cell[2])}.000  90.00  90.00  90.00 P 1           1\n"
    for line in inFile.split("\n"):
        if "HETATM" in line or "ATOM" in line:
            line = f"{line[:30]}{coords[atomIndex][0]:>8.3f}{coords[atomIndex][1]:>8.3f}{coords[atomIndex][2]:>8.3f}{line[54:]}"
            atomIndex += 1
        elif "CRYST1" in line:
            continue
        outFile += f"{line}\n"

    return outFile


if __name__ == "__main__":
    with open("NEWPDB.PDB", "r") as file:
        outFile = makeUnitcell(file.read(), z_ySideLengh=6)
    with open("SHIFTED.PDB", "w") as file:
        file.write(outFile)