import sys
sys.path.append("Modules/")
from job import Job
from excel import Excel

class Task():
    def __init__(self, job):
        self.job : Job = job

    def moveFiles(self):
        pass

    def writeInputFile(self):
        pass

    def submit(self, taskFolder):
        import subprocess
        shell = subprocess.run(["sbatch run_job.sh"], cwd=taskFolder, shell=True)
        print(shell.stdout)

    def _findNNDihedral(self):
        from rdkit import Chem
        from rdkit.Chem.rdDetermineBonds import DetermineBonds
        mol = Chem.MolFromXYZFile(f"{self.job.location}/Orca_Opt/start_molecule.xyz")
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
    
if __name__ == "__main__":
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    ts = Task(jobs[0])
    print(ts._findNNDihedral())