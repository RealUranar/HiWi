import os, shutil, glob, sys

sys.path.append("Modules/Misc")
from InputFileReader import Reader
from Sbatch import JobScripts
from task import Task

from rdkit import Chem
from rdkit.Chem.rdDetermineBonds import DetermineBonds
from rdkit.Chem import AllChem

class MoleculeActions():
    def combineMolecules(mol1, mol2, index):  #Takes the rdkit mols and combines them
        """Function to combine two rdkit-mol objects at a specific point

        Args:
            mol1 (_type_): rdkit structure of molecule Nr.1
            mol2 (_type_): rdkit Structure of molecule Nr.2
            index (_type_): A tuple with the index of the atoms to connect (eg. (12, 24)). IMPORTANT: The order has to match the order of the given Molecules

        Returns:
            _type_: A New Molecule
        """
        combo = Chem.CombineMols(mol1, mol2)  #Combine
        edcombo = Chem.EditableMol(combo)     #Make editabel
        a1, a2 = index
        
        edcombo.AddBond(a1, mol1.GetNumAtoms()+ a2, order=Chem.rdchem.BondType.SINGLE)  #Combine the molecules with a bond
        combinedMol = edcombo.GetMol()  #Get the finished structure
        AllChem.EmbedMolecule(combinedMol)  #Somewhat relax the structure to make a belivable Molecule
        return Molecule(combinedMol)
    
    def Mol2COM(mol): #Takes the rdkit mol and give .com file
        from openbabel import openbabel
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "com")

        molFile =Chem.MolToMolBlock(mol)
        comFile = openbabel.OBMol()
        obConversion.ReadString(comFile, molFile)
        return obConversion.WriteString(comFile)
    

class Molecule():
    def __init__(self, molecule):
        if type(molecule) == str:
            self.mol = Chem.MolFromXYZFile(molecule)
            DetermineBonds(self.mol,charge = 0)
        else:
            self.mol = molecule
    
    def getMol(self):
        return self.mol
    
    def removeAtom(self, atomNr):
        edFrag = Chem.EditableMol(self.mol)
        edFrag.RemoveAtom(atomNr)
        self.mol = edFrag.GetMol()
    


class Gaussian_opt(Task,Reader):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gaussian"
        self.inputFile = ""
        self.executionOrder = [self.moveFiles,
                               self.writeInputFile,
                               self.generateJobScript,
                               self.submit]

    def moveFiles(self):
        os.mkdir(f"{self.job.location}/Gaussian") # Create new SubFolders
        optiFilePath = glob.glob(f"{self.job.location}/Orca_Opt/orca.xyz")[0]  #Get Path to the optimized structure
        shutil.copy(optiFilePath, f"{self.newPath}/orca_opt.xyz") #Copy xyz File to new directory


    def writeInputFile(self):
        comFile = self._setupMolecule()
        with open(f"{self.newPath}/combined.com", "w") as file:
                    file.writelines([
                        "%Chk=F1a_ortho.chk\n",
                        "#P RHF/6-31G* Opt",
                        "\n\nCOMMENT\n\n"])

                    for line in list(comFile.split("\n"))[5:]:
                        file.write(f"{line}\n")

                    file.writelines([
                        "--Link1--\n",
                        "%Chk=F1a_ortho.chk\n",
                        "#P HF/6-31G* SCF=Tight Geom=AllCheck Guess=Read\n",
                        "Pop=MK IOp(6/33=2, 6/41=10, 6/42=17)"
                    ])

    def generateJobScript(self):
        JobScripts().writeGausianJob(name = self.job.id, location=self.newPath)

    def submit(self):
        self.job.updateJob(Gaussian = 2)
        print(f"Submitted Gaussian job {self.name}")
        return super().submit(self.newPath)
        
    def isFinished(self):
        hasFinished, succesfull = False, True
        tail = self._readTail(self.newPath, gauss=True)
        hasFinished = "Normal termination of Gaussian" in tail
        #succesfull = "****ORCA TERMINATED NORMALLY****" in str(tail)
            
        if hasFinished:
            if succesfull:
                self.job.updateJob(Gaussian = 1, Gromacs = 3)
                print(f"Gaussian Job {self.name} has finished succesfull")
            else:
                self.job.updateJob(Gaussian = -1)
                print(f"Gaussian Job {self.name} run into a problem")
        print(f"Gaussian Job {self.name} is still running")

    def _setupMolecule(self):
        inputVars = self.readInputFile(f"{self.job.location}Input")
        molecule = Molecule(f"{self.newPath}/orca_opt.xyz")
        molecule.removeAtom(inputVars["removeAtomNr"])

        fragment = Molecule(f"Modules/fragment.xyz")
        fragment.removeAtom(inputVars["removeAtomFragmentNr"])

        combinedMolecules = MoleculeActions.combineMolecules(molecule.getMol(), fragment.getMol(), (inputVars["combineAtomAt"],inputVars["combineFragmentAt"]))
        comFile = MoleculeActions.Mol2COM(combinedMolecules.getMol())
        return comFile

if __name__ == "__main__":
    from excel import Excel
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Gaussian_opt(jobs[0])
    task.writeInputFile()
    #task.moveFiles()
    #task.generateJobScript()
    #task.submit()
