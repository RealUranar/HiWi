
import xml.etree.ElementTree as ET
import numpy as np
import os

import CreateFiles as CF

class QEReader():
    def __init__(self, PathToData):
        self._pathToData = PathToData
        self._AUTOA = 0.529177249
        self._HARTREETOEV = 27.211396641308
        self._nKpoints = 0

        self._kWeights = np.array([])
        self._atoms = []

        tree = ET.parse(f"{PathToData}/data-file-schema.xml")
        self._root = tree.getroot()

        self._PPFiles = []
        for f in os.listdir(self._pathToData):
            if "UPF" in f:
                self._PPFiles.append(f)


    def getUsedProgram(self):
        return "QuantumEspresso"

    def _findStringInList(self,list, stringToFind):
        """Searches for a string in a list of items\n
        return the first item with a matching string
        """
        for item in list:
            if item.find(stringToFind) != -1:
                return item

    def _getNPArrayFromString(self,listAsString:str):
        """Converts a long string of numbers to a numpy array"""
        return np.array([float(x) for x in listAsString.split()])

    def _StringToMatrix(self,matAsString, x, y):
        matNP = np.zeros((x,y))
        matString = []
        tempNum = [float(x) for x in matAsString.split()]

        for x_i in range(x):
            temp = []
            for y_i in range(y):
                matNP[x_i][y_i] = tempNum[x_i*x + y_i]
                temp.append(tempNum[x_i*x + y_i])
            matString.append(temp)

        return matNP, matString

    def getSym(self):
        noSym = int(bool(self._root.find("./input/symmetry_flags/nosym").text))
        noInv = int(bool(self._root.find("./input/symmetry_flags/noinv").text))
        return noSym, noInv

    def getEnergyForAnalyses(self) -> "tuple[float, float, int]":
        """Returns the:
            Energy minimum  (EMIN),
            Energy maximum  (EMAX) and
            Energy step     (NEDOS)
            """
        #For QE these Values are Fixed
        EMIN = -20.00
        EMAX =  20.00
        NEDOS = 100
        return EMIN, EMAX, NEDOS

    def getAtoms(self):  #Ich speichere Echte Coordinaten und nicht fractional (diff. atom pos und ion pos?)
        """Get atomic information\n
            returns:
            {
                element: str
                pos: [double]
            }
        """
        celldm1 = float(self._root.find("./output/atomic_structure").attrib["alat"])
        
        atoms = []
        coords = []
        for i, atom in enumerate(list((self._root.find("./output/atomic_structure/atomic_positions")))):
            atomPos = np.array([float(x) for x in atom.text.split(" ")])
            atomPos /= celldm1
            _, LatFrac ,_= self.getLattice()
            ionPos =  np.dot(np.transpose(LatFrac), atomPos)

            self._atoms.append(atom.attrib["name"])
            coords.append(list(ionPos))

        #_ionPositions.push_back(recLatticeFrac.transpose()*atomPosition);
        #Was macht das???? 342
        return self._atoms, coords

    def getLattice(self):
        """Get lattice information\n
            returns:\n
                Celllattice real\n
                Celllattice rec\n
                Cell Volume
        """
        CellLatticeReal = np.zeros((3,3))
        CellLatticeRecFrac = np.zeros((3,3))
        for i, vec in enumerate(list(self._root.find("./output/atomic_structure/cell"))):
            coord = vec.text.split(" ")
            CellLatticeReal[i][0] = float(coord[0]) * self._AUTOA
            CellLatticeReal[i][1] = float(coord[1]) * self._AUTOA
            CellLatticeReal[i][2] = float(coord[2]) * self._AUTOA
        
        for i, vec in enumerate(list(self._root.find("./output/basis_set/reciprocal_lattice"))):
            coord = vec.text.split()
            CellLatticeRecFrac[i][0] = float(coord[0])
            CellLatticeRecFrac[i][1] = float(coord[1])
            CellLatticeRecFrac[i][2] = float(coord[2])
        
        
        CellVolume = np.abs(np.dot(np.cross(CellLatticeReal[:,1],CellLatticeReal[:,2]),CellLatticeReal[:,0]))

        return CellLatticeReal, CellLatticeRecFrac, CellVolume

    def getECut(self):
        """get Cutoff Energie"""
        ECut = float(self._root.find("./output/basis_set/ecutwfc").text) * self._HARTREETOEV
        return ECut

    def getFermi(self):
        """get Fermi-energie"""
        fermiEnergy = float(self._root.find("./output/band_structure/fermi_energy").text) * self._HARTREETOEV
        return fermiEnergy

    def getNSpins(self):
        """Get Spinchannels of Calculation"""
        nSpins = (self._root.find("./output/magnetization/lsda").text == "true") +1
        return nSpins

    def getNBands(self):
        """Get number of bands of Calculation"""
        nbnd = int(self._root.find("./output/band_structure/nbnd").text)
        return nbnd

    def getNKpoints(self):
        if self._nKpoints == 0:
            self._nKpoints = int(self._root.find("./output/band_structure/nks").text)

        return self._nKpoints
    
    def getGVector(self):
        gVec = self._root.find("./output/basis_set/fft_smooth")
        sizeOfGVectorGrid = [float(gVec.attrib[x]) for x in gVec.keys()]
        return sizeOfGVectorGrid

    def getKPoints(self, k):
        """
        get Information about K-Points\n
            returns:\n
            for all KPoints:[
                number of PW\n
                Occupation \n
                Eigenvalues \n
                Coords \n
                Weights]
        """
        if len(self._kWeights) == 0:
            for kPoint in self._root.findall("./output/band_structure/ks_energies"):
                self._kWeights = np.append(self._kWeights, float(kPoint[0].attrib["weight"]))
            self._kWeights /= self._kWeights.sum()

        kPoint = self._root.findall("./output/band_structure/ks_energies")[k]

        nPW =int(kPoint[1].text)
        coords =[float(x) for x in kPoint[0].text.split()]
        eigenValues =  list(np.array([float(x) for x in kPoint[2].text.split()]) * self._HARTREETOEV)
        occupations =  [float(x) for x in kPoint[3].text.split()]
        
        return nPW, self._kWeights[k], coords, eigenValues, occupations


    def getPPFileNames(self):
        return self._PPFiles

    def getPPInfo(self, filename):
        PPtree = ET.parse(f"{self._pathToData}/{filename}")
        PProot = PPtree.getroot()

        element = PProot.find("PP_HEADER").attrib["element"].replace(" ", "")
        title = f'{element.lower()}-{PProot.find("PP_HEADER").attrib["functional"]}'
        
        nProj = int(PProot.find("PP_HEADER").attrib["number_of_proj"])
        
        nElec = float(PProot.find("PP_HEADER").attrib["z_valence"])
        
        infoText = PProot.find("PP_INFO").text
        infoText = infoText[infoText.find("Valence configuration:")+80:].split()[:-3]
        orbitals = []
        i = 0
        while True:
            if infoText[i*7].lower() == "generation":
                break
            orbitals.append(infoText[i*7].lower())
            i += 1

        return element, title, nElec, nProj, orbitals

    def getPPGrid(self, filename):
        PPtree = ET.parse(f"{self._pathToData}/{filename}")
        PProot = PPtree.getroot()

        nProj = int(PProot.find("PP_HEADER").attrib["number_of_proj"])

        rGrid = self._getNPArrayFromString(PProot.find("PP_MESH/PP_R").text)
        drGrid = self._getNPArrayFromString(PProot.find("PP_MESH/PP_RAB").text)
        _, auChargeString = self._StringToMatrix(PProot.find("PP_NONLOCAL/PP_AUGMENTATION/PP_Q").text, nProj, nProj)

        return list(rGrid), list(drGrid), auChargeString
        
    def getPPProjektors(self, filename):
        PPtree = ET.parse(f"{self._pathToData}/{filename}")
        PProot = PPtree.getroot()

        rGrid = self._getNPArrayFromString(PProot.find("PP_MESH/PP_R").text)
        nProj = int(PProot.find("PP_HEADER").attrib["number_of_proj"])
        # Kann man auch so machen, falls nProj != anz. BETA 
        # iter = 1  #Kann zu fehlern führen, wenn PP_BETA.1 nicht an zweiter Stelle steht
        # while True:
        #     PP_BETA = PProot.find(f"PP_NONLOCAL/PP_BETA.{iter}")
        #     if PP_BETA == None:
        #         break
        #     angMom = int(PP_BETA.attrib["angular_momentum"])
        #     iter += 1

        projektor = []
        AEpartialwave =[]
        PSpartialwave = []
        angMom = []

        wfConversion = np.power(self._AUTOA, -3.0/2.0)
        for i in range(nProj):
            angMom.append(int(PProot.find(f"PP_NONLOCAL/PP_BETA.{i+1}").attrib["angular_momentum"]))
            projektor.append(list((self._getNPArrayFromString(PProot.find(f"PP_NONLOCAL/PP_BETA.{i+1}").text) / rGrid) * wfConversion))
            AEpartialwave.append(list((self._getNPArrayFromString(PProot.find(f"PP_FULL_WFC/PP_AEWFC.{i+1}").text) / rGrid) * wfConversion))
            PSpartialwave.append(list((self._getNPArrayFromString(PProot.find(f"PP_FULL_WFC/PP_PSWFC.{i+1}").text) / rGrid) * wfConversion))


        return angMom , projektor, AEpartialwave, PSpartialwave


    def getPWCoeffs(self,kpoint, spin):  #https://gitlab.com/QEF/q-e/-/wikis/Developers/Format-of-wfc-files
        with open(f'{self._pathToData}/wfc{kpoint+1}.dat', 'rb') as f:
            # Moves the cursor 60 bytes to the right
            f.seek(60,1)
            igwx = np.fromfile(f, dtype='int32', count=1)[0]  #max number of PW (may be larger than ngw, not sure why)
            npol = np.fromfile(f, dtype='int32', count=1)[0]  #number of spin states for PWs: 2 for non-colinear case, 1 otherwise
            nbnd = np.fromfile(f, dtype='int32', count=1)[0]  #number of wavefunctions
            f.seek(88, 1)

            mill = np.fromfile(f, dtype='int32', count=3*igwx)
            mill = mill.reshape( (igwx, 3) ) #     miller indices: h=mill(1,i), k=mill(2,i), l=mill(3,i)  the i-th PW has wave vector (k+G)(:)=xk(:)+h*b1(:)+k*b2(:)+ l*b3(:)
            evc = np.zeros( (nbnd, npol*igwx), dtype="complex128")

            f.seek(8,1)
            for i in range(nbnd):
                evc[i,:] = np.fromfile(f, dtype='complex128', count=npol*igwx)
                f.seek(8, 1)
            
            evc.flatten(order="C")
            #COMPLEX(8) :: evc(npol*igwx,nbnd)
            #wave functions in the PW basis set, The first index runs on PW components, the second index runs on band states. 
            #For non-colinear case, each PW has a spin component, first  igwx components have PW with   up spin, second igwx components have PW with down spin
            
            return evc , mill


if __name__ == "__main__":
    PathToData = "Data_NaCl_QE/NaCl.save"
    PathToData = "Data_CO/CO.save" 

    writer = CF.Writer(QEReader, PathToData)
    writer.writeData()

    # Test  = QEReader(PathToData)
    # filenames = Test.getPPFileNames()
    # for file in filenames:
    #     print(len(Test.getPPGrid(file)[2]))



















# def getPWCoeffs_old(nFile):  #https://gitlab.com/QEF/q-e/-/wikis/Developers/Format-of-wfc-files
#     with open(f'{PathToData}/wfc{nFile}.dat', 'rb') as f:
#         # Moves the cursor 4 bytes to the right
#         f.seek(4)

#         ik = np.fromfile(f, dtype='int32', count=1)[0]   #k-point index (1 to number of k-points)
#         xk = np.fromfile(f, dtype='float64', count=3)   #k-point coordinates
#         ispin = np.fromfile(f, dtype='int32', count=1)[0]  #spin index for LSDA case: ispin=1 for spin-up, ispin=2 for spin-down for unpolarized or non-colinear cases, ispin=1 always
#         gamma_only = bool(np.fromfile(f, dtype='int32', count=1)[0])   #if .true. write or read only half of the plane waves
#         scalef = np.fromfile(f, dtype='float64', count=1)[0]    # scale factor applied to wavefunctions

#         # Move the cursor 8 byte to the right
#         f.seek(8, 1)

#         ngw = np.fromfile(f, dtype='int32', count=1)[0]   #number of plane waves (PW)
#         igwx = np.fromfile(f, dtype='int32', count=1)[0]  #max number of PW (may be larger than ngw, not sure why)
#         npol = np.fromfile(f, dtype='int32', count=1)[0] #number of spin states for PWs: 2 for non-colinear case, 1 otherwise
#         nbnd = np.fromfile(f, dtype='int32', count=1)[0]  #number of wavefunctions
#         # Move the cursor 8 byte to the right
#         f.seek(8, 1)

#         b1 = np.fromfile(f, dtype='float64', count=3)  #primitive reciprocal lattice vectors
#         b2 = np.fromfile(f, dtype='float64', count=3)
#         b3 = np.fromfile(f, dtype='float64', count=3)

#         f.seek(8,1)
        
#         mill = np.fromfile(f, dtype='int32', count=3*igwx)
#         mill = mill.reshape( (igwx, 3) ) #     miller indices: h=mill(1,i), k=mill(2,i), l=mill(3,i)  the i-th PW has wave vector (k+G)(:)=xk(:)+h*b1(:)+k*b2(:)+ l*b3(:)


#         evc = np.zeros( (nbnd, npol*igwx), dtype="complex128")

#         f.seek(8,1)
#         for i in range(nbnd):
#             evc[i,:] = np.fromfile(f, dtype='complex128', count=npol*igwx)
#             f.seek(8, 1)
        
#         evc.flatten(order="C")
#         #COMPLEX(8) :: evc(npol*igwx,nbnd)
#         #wave functions in the PW basis set, The first index runs on PW components, the second index runs on band states. 
#         #For non-colinear case, each PW has a spin component, first  igwx components have PW with   up spin, second igwx components have PW with down spin
#         with h5py.File(f"kPoint{nFile}.hdf5", "w") as hf:
#             hf.create_dataset('Miller', data=mill)
#             hf.create_dataset('PWCoeffs', data=evc)
#         return evc


