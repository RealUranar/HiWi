import glob, os, sys, shutil

sys.path.append("Modules/Misc")
from InputFileReader import Reader
from excel import Excel
from job import Job

sys.path.append("Modules/Tasks")
from Orca_Opt import Orca_opt
from Orca_Dihedral import Orca_Dihedral
from Gaussian import Gaussian_opt
from GromacsEnergy import GromacsEnergy
from GromacsEquillibration import GromacsEquill
from GromacsProduction import GromacsProd
from Amber import Amber

# Setting up the new Calculation in Excel spread sheet
def setupNewCalculation(NewFolder):
    keyDict = Reader().readInputFile(f"{NewFolder}Input")
    name = keyDict["Name"]
    with Excel() as schedule:
        schedule.createJob(name, location = f"Calculations/{name}/") #Setup New Job in the Excel Spreadsheet

    os.makedirs(f"Calculations/{name}/") #Create new Folder
    files = glob.glob(f"{NewFolder}/*")
    for file in files:
        shutil.copy(file, f"Calculations/{name}/")  #Copy files temporarily
    shutil.rmtree(NewFolder, ignore_errors = False) #Remove input folder

#Handle new Jobs
for NewFolder in glob.glob("Input/*/"):
    setupNewCalculation(NewFolder)

#Get new Jobs to set up
with Excel() as scheduler:
    jobs = scheduler.readJobs()

for job in jobs:
    for runningTask in job.getRunningTask():
        runningTask.isFinished() #Check if a job has finished

    for nextTask in job.getNextTask():  #Execute every function defined in the execute order
        for func in nextTask.executionOrder:
            func()


    with Excel() as scheduler:
        scheduler.updateRow(job)
    
