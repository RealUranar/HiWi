import glob, os, sys, shutil
sys.path.append("Modules/Jobs")
sys.path.append("Modules/Tasks")
from excel import Excel
from job import Job
from Orca_Opt import Orca_opt
from Orca_Dihedral import Orca_Dihedral
from Gaussian import Gaussian_opt
from GromacsEnergy import GromacsEnergy
from GromacsEquillibration import GromacsEquill
from GromacsProduction import GromacsProd
from Amber import Amber

# Setting up the new Calculation in Excel spread sheet
def setupNewCalculation(NewFolder):
    with open(f"{NewFolder}Input", "r") as InputFile:
        Name = InputFile.readline().split("=")[1].strip()
    with Excel() as schedule:
        schedule.createJob(Name, location = f"Calculations/{Name}/") #Setup New Job in the Excel Spreadsheet

    os.makedirs(f"Calculations/{Name}/") #Create new Folder
    files = glob.glob(f"{NewFolder}/*")
    for file in files:
        shutil.copy(file, f"Calculations/{Name}/")  #Copy files temporarily
    shutil.rmtree(NewFolder, ignore_errors = False) #Remove input folder

#Handle new Jobs
for NewFolder in glob.glob("Input/*/"):
    setupNewCalculation(NewFolder)

#Get new Jobs to set up
with Excel() as scheduler:
    jobs = scheduler.readJobs()

for job in jobs:
    for runningTask in job.getRunningTask():
        runningTask.isFinished()

    for nextTask in job.getNextTask():  #Execute every function defined in the execute order
        for func in nextTask.executionOrder:
            func()


    with Excel() as scheduler:
        scheduler.updateRow(job)
    
