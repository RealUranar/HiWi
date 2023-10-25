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

    def _findSubstring(self, smilesString: str, inStructure : str, inFormat = "xyz") -> list:
        """_summary_

        Args:
            smilesString (str): For CSC use "csc", for CNNC use "*N=N*"
            inStructure (str): string of the stucrure you want to input
            inFormat (str, optional): Format of the file. Defaults to "xyz".

        Returns:
            list: _description_
        """        
        from openbabel import pybel
        mol = pybel.readstring(inFormat, inStructure)
        smarts = pybel.Smarts(smilesString)
        subStringIdx = smarts.findall(mol)
        if len(subStringIdx) == 0:
            subStringIdx = pybel.Smarts("*NN*").findall(mol)
            print(f"Could not find {smilesString} string in Structure, using '*NN*'. Make sure this still works!")

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
        dihed = task._findSubstring(smilesString="*N=N*", inStructure=file.read(), inFormat="gro")
    print(dihed, len(dihed))
    # tail = task._readTail(f"{task.job.location}Gromacs", "prod.log")
    # print("Constraint error in algorithm" not in str(tail))




