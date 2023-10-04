import sys, glob
sys.path.append("Modules/")
from job import Job

class Task():
    def __init__(self, job):
        self.job : Job = job

    def moveFiles(self):
        pass

    def writeInputFile(self):
        pass

    def submit(self, taskFolder):
        import subprocess
        shell = subprocess.run("sbatch run_job.sh", cwd=taskFolder, shell=True)
        print(shell.stdout)

    def isFinished(self):
        pass

    def _findNNDihedral(self):
        from rdkit import Chem
        from rdkit.Chem.rdDetermineBonds import DetermineBonds
        mol = Chem.MolFromXYZFile(f"{self.job.location}Orca_Opt/startMolecule.xyz")
        DetermineBonds(mol,charge = 0)

        for bond in mol.GetBonds():  #Find the correct Bond between two N=N
            if bond.GetBeginAtom().GetSymbol() == "N" and bond.GetEndAtom().GetSymbol() == "N":
                break
            
        N1 = bond.GetBeginAtom()
        N2 = bond.GetEndAtom()
        atomsIndex = [N1.GetIdx()]  #Save the Index of the first N
        
        for atom in N1.GetNeighbors():  #Check the neighbors of the first N and save the Index of the bonded atom.
            if atom.GetSymbol() != "N":
                atomsIndex.insert(0, atom.GetIdx())
            if atom.GetSymbol() == "N":
                atomsIndex.insert(2, atom.GetIdx())

        for atom in N2.GetNeighbors():  # Find the index of the last non N-atom
            if atom.GetSymbol() != "N":
                atomsIndex.insert(3, atom.GetIdx())
                
        return atomsIndex

    def _readTail(self,folder, gauss=False):
            tail = ""
            try:
                if gauss:
                    folder = glob.glob(f"{folder}*.log")[0]
                else:
                    folder = glob.glob(f"{folder}output.*.txt")[0]

            except IndexError:
                if gauss:
                    print(f"File {folder}*.log not found!")
                else:
                    print(f"File {folder}output.*.txt not found!")
                
                return tail
            
            try:
                with open(folder, "r") as file:
                    tail = str(file.readlines()[-15:])
            except:
                with open(folder, "r", encoding = "ISO-8859-1") as file:
                    tail = str(file.readlines()[-15:])
            return tail