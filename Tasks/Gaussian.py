import sys, os, shutil, glob
sys.path.append("../HPC_Jobs/")

from Sbatch import JobScripts
from task import Task
from excel import Excel

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
    


class Gaussian_opt(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gaussian"
        #self.inputFile = "" #Looks dead no idea what it was for...

    def moveFiles(self):
        os.mkdir(f"{self.job.location}/Gaussian") # Create new SubFolders
        optiFilePath = glob.glob(f"{self.job.location}/Orca_Opt/orca.xyz")[0]  #Get Path to the optimized structure
        shutil.copy(optiFilePath, f"{self.newPath}/orca_opt.xyz") #Copy xyz File to new directory

    def _setupMolecule(self, removeAtomNr = 22, removeAtomFragmentNr = 30, combineAt = (19,29)):
        molecule = Molecule(f"{self.newPath}/orca_opt.xyz")
        molecule.removeAtom( removeAtomNr)

        fragment = Molecule(f"/home/seifert/src/HiWi/Tasks/fragment.xyz")
        fragment.removeAtom(removeAtomFragmentNr)

        combinedMolecules = MoleculeActions.combineMolecules(molecule.getMol(), fragment.getMol(), combineAt)
        comFile = MoleculeActions.Mol2COM(combinedMolecules.getMol())
        return comFile


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
        JobScripts().writeGausianJob(name = self.job.name, location=self.newPath)

    def submit(self):
        return super().submit(self.newPath)
        

if __name__ == "__main__":
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Gaussian_opt(jobs[0])
    #task.moveFiles()
    #task.writeInputFile()
    task.generateJobScript()
    #task.submit()
