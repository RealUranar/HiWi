import xml.etree.ElementTree as ET
import numpy as np

import CreateFiles as CF


class VASPReader():
    def __init__(self, PathToData):
        self._pathToData = PathToData
        self._AUTOA = 0.529177249
        self._HARTREETOEV = 27.211396641308
        self._nKpoints = 0


        tree = ET.parse(f"{PathToData}/vasprun.xml")
        self._root = tree.getroot()

        

    def _FindIndexOfKeyword(self,content, findString) -> "list[int]":
        """Returns the index of every occurenc of a specific keyword in a textfile"""
        stringLen = len(findString)

        iLast = 0
        indexList = []
        while True:
            iNew = content[iLast:].find(findString)
            if iNew == -1:
                break
            iNew += stringLen + iLast
            indexList.append(iNew)
            iLast = iNew

        return indexList

    def _FindInText(self, content, findString, endString) -> "list[str]":
        """Returns the text between two given keywords in a file.
        If a keyword is used multiple times, it returns every occurence in a list"""
        indexList = self._FindIndexOfKeyword(content, findString)

        strings = []
        for i in indexList:
            string = ""
            for c in content[i:]:
                if c == endString:
                    break
                string += c
            strings.append(string)

        return strings



    def getUsedProgram(self):
        return "VASP"

    def getSym(self) -> "tuple[int, int]":
        """Returns the Symmetry of the System, only used to determine if the calculation was done correcly"""
        noSym = int(self._root.find("parameters//*[@name='symmetry']/*[@name='ISYM']").text)
        if noSym != 0 or noSym != -1:
            noSym = 1

        noInv = 0 #Nur für andere Programme verwendet
        return noSym, noInv

    def getEnergyForAnalyses(self) -> "tuple[float, float, int]":
        """Returns the:
            Energy minimum  (EMIN),
            Energy maximum  (EMAX) and
            Energy step     (NEDOS)
            """
        energyPath = self._root.find("parameters//*[@name='dos']")
        EMIN = float(energyPath.find("*[@name='EMIN']").text)
        EMAX = float(energyPath.find("*[@name='EMAX']").text)
        NEDOS = int(energyPath.find("*[@name='NEDOS']").text)
        return EMIN, EMAX, NEDOS

    def getECut(self) -> float:
        """get Cutoff Energie for the bands"""
        try:
            ECut = float(self._root.find("incar/*[@name='ENCUT']").text)
        except:
            ECut = float(self._root.find("parameters//*[@name='electronic']/*[@name='ENMAX']").text)
        return ECut

    def getGVector(self) -> "list[float]":
        gVec = self._root.findall("parameters//*[@name='grids']/")
        sizeOfGVectorGrid = [float(x.text) for x in gVec[0:3]] #Only take NGX, NGY, NGZ
        return sizeOfGVectorGrid

    def getFermi(self) -> float:
        """get Fermi-energie"""
        fermiEnergy = float(self._root.find("calculation/dos/*[@name='efermi']").text)
        return fermiEnergy

    def getAtoms(self)-> "tuple['list[str]','list[float]']":  #Ich speichere Echte Coordinaten und nicht fractional (diff. atom pos und ion pos?)
        """Get atomic information\n
            returns:
                [Element names] [string]\n
                [Element Positions] [float]
        """
        posTree = self._root.findall("structure[@name='finalpos']/*[@name='positions']/")
        atomTree = self._root.findall("atominfo/array[@name='atoms']/set/")

        atomPos = [[float(x) for x in pos.text.split()] for pos in posTree]
        atomElement = [atom[0].text for atom in atomTree]

        return atomElement, atomPos

    def getLattice(self) -> "tuple[list[float], list[float], float]":
        """Get lattice information\n
            returns:\n
                Celllattice real\n
                Celllattice rec\n
                Cell Volume
        """
        CellLatticeReal = np.zeros((3,3))
        CellLatticeRecFrac = np.zeros((3,3))
        for i, vec in enumerate(self._root.findall("structure[@name='finalpos']/crystal/*[@name='basis']/")):
            coord = vec.text.split()
            CellLatticeReal[0][i] = float(coord[0])
            CellLatticeReal[1][i] = float(coord[1])
            CellLatticeReal[2][i] = float(coord[2])
        
        for i, vec in enumerate(list(self._root.findall("structure[@name='finalpos']/crystal/*[@name='rec_basis']/"))):
            coord = vec.text.split()
            CellLatticeRecFrac[0][i] = float(coord[0])
            CellLatticeRecFrac[1][i] = float(coord[1])
            CellLatticeRecFrac[2][i] = float(coord[2])
        
        CellVolume = float(self._root.find("structure[@name='finalpos']/crystal/*[@name='volume']").text)

        return CellLatticeReal, CellLatticeRecFrac, CellVolume

    def getNBands(self) -> int:
        """Get number of bands of Calculation"""
        nbnd = int(self._root.find("parameters//*[@name='electronic']/*[@name='NBANDS']").text)
        return nbnd

    def getNSpins(self) -> int:
        """Get Spinchannels of Calculation"""
        nSpins = int(self._root.find("parameters//*[@name='electronic spin']/*[@name='ISPIN']").text)
        return nSpins

    def getNKpoints(self) -> int:
        if self._nKpoints == 0:
            self._nKpoints = len(self._root.findall("kpoints/*[@name='kpointlist']/"))
        
        return self._nKpoints

    def getKPoints(self, kpoint) -> "tuple[int, float, 'list[float]', 'list[float]', 'list[float]']":
        """
        BISHER NUR FÜR EINEN SPIN!!!!\n
        get Information about K-Points\n
            returns:\n
            Number of PWs at KPoint (int)\n
            Weight of Kpoint (float)\n
            Coordinates of KPoint (list(float))\n
            Eigenvalues of KPoint (list(float))\n
            Occupations of Kpoint (list(float))\n
        """
        with open(f"{self._pathToData}/OUTCAR") as file:
            content = file.read()
            nPW = [int(x) for x in self._FindInText(content, "plane waves:", "\n")]
                
        kPos = self._root.findall("kpoints/*[@name='kpointlist']/")
        kweight = self._root.findall("kpoints/*[@name='weights']/")
        kEigenAndOccs = self._root.findall("calculation/eigenvalues/array/set/set/")

        coords = [float(x) for x in kPos[kpoint].text.split()]
        weight = float(kweight[kpoint].text)
        eigen = [float(x.text.split()[0]) for x in kEigenAndOccs[kpoint]]
        occs =  [float(x.text.split()[1]) for x in kEigenAndOccs[kpoint]]
        return nPW[kpoint], weight, coords,  eigen, occs

    def getPPKeys(self):
        with open(f"{self._pathToData}/POTCAR") as potcar:
            content = potcar.read()
        PPs = {}
        iStart = 0
        for i in self._FindIndexOfKeyword(content, "End of Dataset"):
            elementIndex = content[iStart:].find("VRHFIN =") + 8 + iStart
            element = content[elementIndex:elementIndex+3].replace(":","").replace(" ","")
            PPs[element] = content[iStart:i]
            iStart=i+1

        self._PPs = PPs
        return PPs.keys()

    def getPPInfo(self, element):
        #Valence fehlt, da nicht in daten gegeben
        PPContent = self._PPs[element]

        PPSplit = PPContent.split("\n")
        title = PPSplit[0][1:].replace("  ", "")
        nElec = float(PPSplit[1])

        orbitals = [] #TODO!!!!!!!!
        return element ,title, nElec, orbitals

    def getPPGrids(self, element):
        PPContent = self._PPs[element]

        indexNonlocal = self._FindIndexOfKeyword(PPContent, "Non local Part")
        maximumGridValueOfProjectorRec = float(PPContent[indexNonlocal[0]-100:indexNonlocal[0]-20].split("\n")[-1])

        numberOfComingProj = []
        angularQuantumNumberOfProjector = []
        for nonLocalPart in indexNonlocal:
            firstRow = PPContent[nonLocalPart:].split("\n")[1].split()
            print(firstRow)
            angularQuantumNumberOfProjector.append(int(firstRow[0]))
            numberOfComingProj.append(int(firstRow[1]))
            maximumGridValueOfProjectorReal  = float(firstRow[2])

        indexReciprocal = self._FindIndexOfKeyword(PPContent, "Reciprocal Space Part")

        recReal = {}
        i = 0
        for angmom, nproj in zip(angularQuantumNumberOfProjector, numberOfComingProj):
            rec = []
            real = []
            for _ in range(nproj):
                index = indexReciprocal[i]
                rec.append([float(x) for x in PPContent[index:].split()[:100]])
                real.append([float(x) for x in PPContent[index:].split()[103:203]])
                i += 1

            recReal[angmom] = [rec,real]

        return recReal, maximumGridValueOfProjectorRec, maximumGridValueOfProjectorReal

    def getPPProjectors(self,element):
        PPContent = self._PPs[element]
        indexRadialSets = self._FindIndexOfKeyword(PPContent, "PAW radial sets")[0]
        
        radialSets = PPContent[indexRadialSets:].split()

        augmentationCharge = np.array([])
        for val in radialSets[7:]:
            try:
                val = float(val)
            except:
                break
            augmentationCharge = np.append(augmentationCharge, val)

        nProj = int(np.sqrt(len(augmentationCharge)))

        augmentationCharge = augmentationCharge.reshape(nProj,nProj)
        augmentationChargeList = [list(x) for x in augmentationCharge]
       
        gridSize = int(radialSets[0])
        

        indexGrids = self._FindIndexOfKeyword(PPContent, "grid\n")[0]
        grid = np.array([float(x) for x in PPContent[indexGrids:].split()[0:gridSize]])

        indexPSWave = self._FindIndexOfKeyword(PPContent, "pseudo wavefunction")
        psWave = []
        aeWave = []
        for i in indexPSWave:
            psTemp = np.array([float(x) for x in PPContent[i:].split()[0:gridSize]]) / grid
            aeTemp = np.array([float(x) for x in PPContent[i:].split()[gridSize+2:gridSize*2 +2]]) / grid
            psWave.append(list(psTemp))
            aeWave.append(list(aeTemp))

        return augmentationChargeList, nProj, psWave, aeWave, list(grid)
        
    def getPWCoeffs(self, kpoint, spin): 
        """Returns the Plane Wave coefficients for a specific kpoint and spin from the WAVECAR"""

        def whereRec(nKpoints, nBands, iSpin, ikpt, iband):
            rec = 2 + (iSpin -1) * nKpoints * (nBands +1) + (ikpt -1) * (nBands +1) + iband
            return int(rec) -1
            
        with open(f'{self._pathToData}/WAVECAR', 'rb') as f:
            f.seek(0,0)
            recl, nSpin, rtag = map(int, np.fromfile(f, dtype="float64", count=3))
            
            #recl = record lengh
            #rtag = coefficint Precision 
            #   45200 single precision : np.complex64
            #   45210 double precision : np.complex128
            if rtag == 45200.0:
                prec = np.complex64
            else:
                prec = np.complex128

            f.seek(recl)
            nKpoints, nBands = map(int, np.fromfile(f, dtype=np.float, count=2))
            ECut, *LatVec, EFermi = np.fromfile(f, dtype=np.float, count=11)
            

            rec = whereRec(nKpoints, nBands, spin+1, kpoint+1, 1)
            f.seek(recl * rec,0)
            nPW, *kvecs = np.fromfile(f, dtype=np.float,count=4)
            nPW = int(nPW)

            EigenVals, _ , Occs = np.fromfile(f, dtype=np.float,count=3*nBands).reshape(nBands,3).transpose()


            f.seek(recl * rec+ recl,0)
            PWCoeffs = np.empty((0,nPW))
            for _ in range(nBands):
                PWCoeffsBand = np.fromfile(f, dtype=prec, count=nPW)
                PWCoeffs = np.append(PWCoeffs, PWCoeffsBand)
                f.seek(recl - nPW*8 ,1)
        
            PWCoeffs = PWCoeffs.reshape(nBands, nPW)

        return PWCoeffs, None

if __name__ == "__main__":
    PathToData = "Data_H2O_VASP"
    #PathToData = "Data_CO_VASP"

    writer = CF.Writer(VASPReader, PathToData)
    writer.writeData()
    
    # Test = VASPReader(PathToData)
    # PPKeys = Test.getPPKeys()
    # for element in PPKeys:
    #     None
    
    # PwCoeffs, _ = Test.getPWCoeffs(0,0)
    # print(PwCoeffs)

    # for k in range(Test.getNKpoints()):
    #     PwCoeffs, _ = Test.getPWCoeffs(k,0)
        # print(len(PwCoeffs))
        # count = 0
        # for val in PwCoeffs[0]:
        #     if val == 0:
        #         count += 1
        # print(count)
        # print(len(PwCoeffs[0][-1]))
        # #PwCoeffs[8] = PwCoeffs[8][:count]