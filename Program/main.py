import glob, os, sys, shutil

sys.path.append("Modules/Misc")
from InputFileReader import Reader
from job import Job
from database import JobDatabase

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
    inputVars = Reader(f"{NewFolder}Input")
    name = inputVars.getKeyword("Name")
    id = JobDatabase.getLastID("Database.json") + 1
    taks = inputVars.getKeyword("Tasks")[1]
    newLocation = f"/hpcwork/om962181/Calculations/{name}/"

    JobDatabase.saveJob({"id": id,
                        'name': name,
                        'location': newLocation,
                        "tasks": {
                            'waitingtasks': taks,
                            "runningtasks": [],
                            "finnishedtasks": [],
                            "failedtasks" : []}},
                        filepath= 'Database.json')

    os.makedirs(newLocation) #Create new Folder
    files = glob.glob(f"{NewFolder}/*")
    for file in files:
        shutil.copy(file, newLocation)  #Copy files temporarily
    shutil.rmtree(NewFolder, ignore_errors = False) #Remove input folder

#Handle new Jobs
for NewFolder in glob.glob("Input/*/"):
    setupNewCalculation(NewFolder)

#Get new Jobs to set up
jobs = []
for job in JobDatabase.loadJobs("Database.json"):
    jobs.append(Job(job))

for job in jobs:
    for runningTask in job.getRunningTask():
        runningTask.isFinished() #Check if a job has finished

    if job.getRunningTask() == [] and job.getNextTask() != []: #If no job is running and there is a job to run
        for func in job.getNextTask()[0].executionOrder:
            func()  #Execute every function defined in the execute order


    JobDatabase.saveJob(job.toDict(), "Database.json") #Save the job to the database
