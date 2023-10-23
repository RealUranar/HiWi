import glob, os, sys
import subprocess
sys.path.append("Modules/Misc")
from convertFile import convertFile
from rdkit.Chem.rdchem import Mol

class Task():
    def __init__(self, job):
        self.job : Job = job
        self.newPath = ""
        self.executionOrder = []

    def moveFiles(self):
        pass

    def writeInputFile(self):
        pass

    def submit(self, taskFolder, joScriptName = "run_job.sh"):
        shell = subprocess.run(f"sbatch {joScriptName}", cwd=taskFolder, shell=True)
        print(shell.stdout)

    def isFinished(self):
        pass

    def _findSubstring(self, smilesString: str, inStructure : str|Mol, inFormat = "xyz") -> list:
        from rdkit import Chem
        from rdkit.Chem.rdmolfiles import MolFromXYZBlock
        if inFormat != "xyz":
            inStructure = convertFile(inStructure, inFormat=inFormat, outFormat="xyz")

        if type(inStructure) == Mol:
            import copy
            mol = copy.deepcopy(inStructure)
        else:
            mol = MolFromXYZBlock(inStructure)  #Convert File to rdkit Structure

        def setBondsManual(mol, newBondbetween):
            for bond in mol.GetBonds():  #Find the correct Bond between two N=N
                if bond.GetBeginAtom().GetSymbol() == newBondbetween[0] and bond.GetEndAtom().GetSymbol() == newBondbetween[1]:
                    bond.SetBondType(Chem.rdchem.BondType.DOUBLE)

        try:  #Here try/except is nessecary, because Determine bonds does not like my 4x4 cell
            if mol.GetNumAtoms() > 100:
                raise TimeoutError("Too many atoms!")
            from rdkit.Chem.rdDetermineBonds import DetermineBonds
            DetermineBonds(mol,charge = 0)
            if smilesString == "CN=NC":
                smilesString = 'cN=Nc'
        except:
            from rdkit.Chem.rdDetermineBonds import DetermineConnectivity
            DetermineConnectivity(mol,charge = 0)  #This is fine for small molecules, but does not work for big ones
            if "=" in smilesString or "#" in smilesString or "$" in smilesString:
                atoms = smilesString.split("=")
                setBondsManual(mol, newBondbetween=(atoms[0][-1], atoms[1][0]))

        subStringIdx = list(mol.GetSubstructMatches(Chem.MolFromSmarts(smilesString)))
        if len(subStringIdx) > 16:
            print("Something went very wrong when looking up the dihedral angle!!\nTHIS IS VERY BAD LOOK AT PLUMED AND table_fourier!!!!!")
            subStringIdx = [subStringIdx[4], "NO", "NO", "NO", "NO", "NO",subStringIdx[4*6]]  #THIS IS VERY BAD 

        return subStringIdx

    def _readTail(self,folder, file = "output.*.txt"):
            try:
                folder = glob.glob(f"{folder}/{file}")[0]

            except IndexError:
                print(f"File {folder}/{file} not found!")
                return ""
            
            try:
                with open(folder, "r") as file:
                    tail = str(file.readlines()[-15:])
            except:
                with open(folder, "r", encoding = "ISO-8859-1") as file:
                    tail = str(file.readlines()[-15:])
            return tail
    
    def _runCondaScript(self,script,taskFolder, environmentName = 'AmberTools23'):
        activate_script = os.path.join("~","miniconda3", 'etc', 'profile.d', 'conda.sh')
        command = "source " + activate_script + ' && conda activate ' + environmentName + f"&& {script}"
        ret = subprocess.run(command, executable='/bin/bash', shell=True, capture_output=True, cwd=taskFolder)
        if ret.returncode != 0:
            print(ret.stdout)
            raise RuntimeError("Some AmberTools Command Failed!")

if __name__ == "__main__":
    import sys
    sys.path.append("Modules/Misc")
    from job import Job
    job = Job(name = "Test", id = 666, location="Calculations/TESTING/", tasks={"Amber":1})

    task = Task(job)
    with open(f"{task.job.location}Gromacs/System.gro", "r") as file:
        dihed = task._findSubstring(smilesString="CN=NC", inStructure=file.read(), inFormat="gro")
    print(dihed)
    # tail = task._readTail(f"{task.job.location}Gromacs", "prod.log")
    # print("Constraint error in algorithm" not in str(tail))




