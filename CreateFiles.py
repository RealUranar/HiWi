import tarfile
import os
import h5py
from shutil import rmtree
import json
from multiprocessing import Pool
from itertools import repeat
from subprocess import call

class Writer():
    def __init__(self, program, pathToFile):
        self._reader = program(pathToFile)
        self._pathToFile = pathToFile

        self._usedProgram = self._reader.getUsedProgram()
        self._nSpins = self._reader.getNSpins()
        self._nKpoints = self._reader.getNKpoints()
        self._latReal, _ ,self._cellVol = self._reader.getLattice()


    def _setupKPointData(self):
        kPointData = []
        for k in range(self._nKpoints):
            nPW, weight, coords, eigen, occs = self._reader.getKPoints(k)
            kPointData.append(
                {
                    "coordinates": coords,
                    "weight": weight,
                    "occupations": occs,
                    "eingen_values": eigen,
                    "n_plane_waves": nPW
                }
            )
        return kPointData

    def _setupAtomData(self):
        atomData = []
        atoms, positions = self._reader.getAtoms()
        for atom, pos in zip(atoms, positions):
            atomData.append(
                {
                    "element": atom,
                    "coordinates": pos
                }
            )
        return atomData

    def _PPQE(self):
        PPFiles = self._reader.getPPFileNames()
        Data = []
        for fileName in PPFiles:
            element, title, nElec, nProj, orbitals = self._reader.getPPInfo(fileName)
            rGrid, drGrid, augCharge = self._reader.getPPGrid(fileName)
            angMom , projektor, AEpartialwave, PSpartialwave = self._reader.getPPProjektors(fileName)
            
            PP_Data = {
                "element": element,
                "title" : title,
                "n_electrons": nElec,
                "n_projections": nProj,
                "angularMomentum" : angMom,
                "valence": orbitals,
                "augmentationCharges" : augCharge,
                "grid" : rGrid,
                "drgrid" : drGrid,
                "projectors": projektor,
                "all_electron_wfc": AEpartialwave,
                "PSpartialwave": PSpartialwave
            }
            Data.append(PP_Data)
        return Data

    def _PPVASP(self):
        PPKeys = self._reader.getPPKeys()
        Data = []
        for element in PPKeys:
            augCharge, nProj, PSpartialwave, AEpartialwave, grid = self._reader.getPPProjectors(element)
            recReal, maxGridValueProjRec, maxGridValueProjReal = self._reader.getPPGrids(element)
            element, title, nElec, orbitals = self._reader.getPPInfo(element)

            angularMom = []
            for angmom in recReal.keys():
                angularMom.append(
                    {
                        "angularMomentum":angmom,
                        "rec": recReal[angmom][0],
                        "real": recReal[angmom][1],
                    }
                )

            PP_Data = {
                "element": element,
                "title" : title,
                "n_electrons": nElec,
                "n_projections": nProj,
                "angularQuantumProjectors" : angularMom,
                "valence":  orbitals,
                "maximumGridValueOfProjectorRec" : maxGridValueProjRec,
                "maximumGridValueOfProjectorReal": maxGridValueProjReal,
                "grid": grid,
                "augmentationCharges" : augCharge,
                "all_electron_wfc": AEpartialwave,
                "PSpartialwave": PSpartialwave
            }
            Data.append(PP_Data)

        return Data

    def _setupPP(self):
        if self._usedProgram == "QuantumEspresso":
            Data = self._PPQE()

        if self._usedProgram == "VASP":
            Data = self._PPVASP()
        return Data

    def _setupJSONFormat(self):
        DATA = {
            "file":{
                "usedProgram": self._usedProgram,
                "noSym" : self._reader.getSym()[0],
                "noInv": self._reader.getSym()[1],
            },
            "cell": {
                "lattice_vectors_real": [list(x) for x in self._latReal],
                "cell_volume": self._cellVol, 
                "atomic_structure" : self._setupAtomData()
            },
            "wave_function": {
                "EnergyMinimumForAnalysis" : self._reader.getEnergyForAnalyses()[0],
                "EnergyMaximumForAnalysis" : self._reader.getEnergyForAnalyses()[1],
                "EnergyStepForAnalysis" : self._reader.getEnergyForAnalyses()[2],
                "sizeOfGvectorGrid": self._reader.getGVector(),
                "cutoff_Energy": self._reader.getECut(),
                "fermi_Energy": self._reader.getFermi(),
                "n_Bands" : self._reader.getNBands(),
                "n_Spins" : self._nSpins,
                "spin_channels": [
                    {
                        "k_points": self._setupKPointData()
                    }
                ]
            },
            "paw": self._setupPP()
        }
        return DATA
    


    def _writeHDF5(self, work):
        k,s = work

        fileName = f"LOBSTER_Kpoints/kPoint{k+1}.hdf5"

        PWCoeffs, Miller = self._reader.getPWCoeffs(k, s)
        with h5py.File(fileName, "w") as hf:
            hf.create_dataset('PWCoeffs', data= PWCoeffs)
            if type(Miller) != type(None):
                hf.create_dataset('Miller', data=Miller)

            return
            

    def writeData(self):
        Data = self._setupJSONFormat()
        print("Writing 'LobsterInput.json'....")
        with open("LobsterInput.json", "w") as outfile:
            json.dump(Data, outfile, indent=4)

        try:
            print("Trying to delete leftover LOBSTER_Kpoints folder")
            rmtree("LOBSTER_Kpoints")
        except:
            print("Not found\n")
        
        print("Writing Kpoints to hdf5 files....")
        os.mkdir("LOBSTER_Kpoints")


        for s in range(self._nSpins):
            print(f"Spin: {s}")
            print(f"Writing Kpoints: {0} of {self._nKpoints}")

            with Pool() as p:
                work = [(k,s) for k,s in zip(list(range(self._nKpoints)) ,repeat(s))]
                res = p.imap_unordered(self._writeHDF5, work)
                p.close()

                completed = 0
                while(True):
                    if completed == res._index:
                        continue
                    completed = res._index
                    if (completed == self._nKpoints):
                        break
                    print ("\033[A                             \033[A")
                    print(f"Writing Kpoints: {completed} of {self._nKpoints}")
                
                        
                    
            print ("\033[A                             \033[A")
            print(f"Writing Kpoints: {self._nKpoints} of {self._nKpoints}")

            print("Adding Kpoints to tar-File (this may take some time)...")


            script = b"#!/bin/bash\n \
                    tar -zcf LobsterInput.tar.gz LOBSTER_Kpoints LobsterInput.json\n"
            try:
                rc = call(script, shell= True)
            except:
                with tarfile.open("LobsterInput.tar.gz", "w:gz") as tar_handle:
                    tar_handle.add("LobsterInput.json")
                    tar_handle.add(f"LOBSTER_Kpoints")

        os.remove("LobsterInput.json")
        rmtree("LOBSTER_Kpoints")
        
        print("Done!")






