import sys, glob
sys.path.append("Modules/")
#from job import Job

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
        import subprocess
        shell = subprocess.run(f"sbatch {joScriptName}", cwd=taskFolder, shell=True)
        print(shell.stdout)

    def isFinished(self):
        pass

    def _findNNDihedral(self, inFile = "Orca_Opt/startMolecule.xyz") -> list:
        from rdkit import Chem
        def convertFile(fileName, inFormat = "gro", outFormat = "xyz"): ##Function to convert a given file to xyz-Format
            with open(fileName, "r") as file:
                inFile = file.read()

            from openbabel import openbabel
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats(inFormat, outFormat)
            outFile = openbabel.OBMol()
            obConversion.ReadString(outFile, inFile)
            return obConversion.WriteString(outFile)

        def determineBondsManual(mol):
            from rdkit.Chem.rdDetermineBonds import DetermineConnectivity
            DetermineConnectivity(mol,charge = 0)
            NNBonds = []
            for bond in mol.GetBonds():  #Find the correct Bond between two N=N
                if bond.GetBeginAtom().GetSymbol() == "N" and bond.GetEndAtom().GetSymbol() == "N":
                    NNBonds.append(bond)
                
            for NNBond in NNBonds:
                NNBond.SetBondType(Chem.rdchem.BondType.DOUBLE)
                a1 = NNBond.GetBeginAtom()     ####The code works even without this section... i have no idea why
                for otherBond in a1.GetBonds():
                    if otherBond.GetBondType() == NNBond.GetBondType():
                        continue
                    otherBond.SetBondType(Chem.rdchem.BondType.SINGLE) ####This scares me...

        if inFile.endswith(".xyz"):
            mol = Chem.MolFromXYZFile(f"{self.job.location}{inFile}")
        else:
            end = inFile.split(".")[-1]
            data = convertFile(f"{self.job.location}{inFile}", inFormat=end)
            mol = Chem.MolFromXYZBlock(data)  ##Convert the input file to xyz then open it


        try:
            from rdkit.Chem.rdDetermineBonds import DetermineBonds
            if mol.GetNumAtoms() > 100:
                raise Exception("Too many Atoms")
            DetermineBonds(mol,charge = 0)
            dihedralIdx = list(mol.GetSubstructMatches(Chem.MolFromSmarts('cN=Nc'))) #Use RDkits inbuilt function to determine the CNNC dihedrals
        except:
            determineBondsManual(mol)
            dihedralIdx = list(mol.GetSubstructMatches(Chem.MolFromSmarts('CN=NC')))
            
        return dihedralIdx


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
    
if __name__ == "__main__":
    sys.path.append("Modules/Misc")
    from job import Job
    job = Job(name = "Test", id = 666, location="Modules/TESTING/", tasks={"Amber":1})

    task = Task(job)
    hi = task._findNNDihedral(inFile="Amber/System.gro")


