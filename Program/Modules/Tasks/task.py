import glob, os
import subprocess

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

    def _findSubstring(self, smilesString, inFile = "Orca_Opt/startMolecule.xyz") -> list:
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

        def setBondsManual(mol, newBondbetween):
            for bond in mol.GetBonds():  #Find the correct Bond between two N=N
                if bond.GetBeginAtom().GetSymbol() == newBondbetween[0] and bond.GetEndAtom().GetSymbol() == newBondbetween[1]:
                    bond.SetBondType(Chem.rdchem.BondType.DOUBLE)

        if inFile.endswith(".xyz"):
            mol = Chem.MolFromXYZFile(f"{self.job.location}{inFile}")
        else:
            end = inFile.split(".")[-1]
            data = convertFile(f"{self.job.location}{inFile}", inFormat=end)
            mol = Chem.MolFromXYZBlock(data)  ##Convert the input file to xyz then open it


        try:
            if mol.GetNumAtoms() > 100:
                raise TimeoutError("Too many atoms!")
            from rdkit.Chem.rdDetermineBonds import DetermineBonds
            DetermineBonds(mol,charge = 0)
            if smilesString == "CNNC":
                smilesString = 'cN=Nc'
        except:
            from rdkit.Chem.rdDetermineBonds import DetermineConnectivity
            DetermineConnectivity(mol,charge = 0)
            if "=" in smilesString or "#" in smilesString or "$" in smilesString:
                atoms = smilesString.split("=")
                setBondsManual(mol, newBondbetween=(atoms[0][-1], atoms[1][0]))

        subStringIdx = list(mol.GetSubstructMatches(Chem.MolFromSmarts(smilesString)))
            
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
    job = Job(name = "Test", id = 666, location="Modules/TESTING/", tasks={"Amber":1})

    task = Task(job)
    hi = task._findSubstring(inFile="Amber/System.gro")


