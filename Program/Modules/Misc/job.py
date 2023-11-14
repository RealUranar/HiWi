import sys
sys.path.append("Modules/Tasks")
from task import Task
from Orca_Opt import Orca_opt
from Orca_Dihedral import Orca_Dihedral
from Gaussian import Gaussian_opt
from GromacsEnergy import GromacsEnergy
from GromacsEquillibration import GromacsEquill
from GromacsProduction import GromacsProd
from Amber import Amber


class Job():
    def __init__(self,dictIn):
        for key in dictIn.keys():
            setattr(self, key, dictIn[key])

        self._convertTasks()

    def toDict(self):
        dictOut = {}
        for key in self.__dict__.keys():
            if key != "tasks":
                dictOut[key] = self.__dict__[key]
            else:
                dictOut[key] = {}
                for task in self.tasks:
                    dictOut[key][task] = []
                    for program in self.tasks[task]:
                        dictOut[key][task].append(type(program).__name__)

        return dictOut


    def _convertTasks(self):
        taskToProgram = {
            "Orca_opt": Orca_opt,
            "Orca_Dihedral": Orca_Dihedral,
            "Gaussian_opt" : Gaussian_opt,
            "Amber" : Amber,
            "GromacsEnergy": GromacsEnergy,
            "GromacsEquill" : GromacsEquill,
            "GromacsProd": GromacsProd}
        
        for task in self.tasks:
            convertedPrograms = []
            for program in self.tasks[task]:
                convertedPrograms.append(taskToProgram[program](self))
            self.tasks[task] = convertedPrograms

    def getNextTask(self) -> list[Task]:
        return self.tasks["waitingtasks"]
    
    def getFinishedTask(self) -> list[Task]:
        return self.tasks["finnishedtasks"]
    
    def getRunningTask(self) -> list[Task]:
        return self.tasks["runningtasks"]

    def getLocation(self):
        return self.location

    def updateJob(self, *args, **kwargs):
        task = list(kwargs.keys())[0]
        value = list(kwargs.values())[0]
        if task not in list(self.tasks.keys()):
            print("Invalid task!")
            return None
        self.tasks[task].append(value)
        if task == "finnishedtasks":
            self.tasks["runningtasks"].remove(value)
        elif task == "runningtasks":
            self.tasks["waitingtasks"].remove(value)
        elif task == "failedtasks":
            self.tasks["failedtasks"].remove(value)

if __name__ == "__main__":
    from database import JobDatabase
    # JobDatabase.saveJob({"id": 1,'name': 'test', 'location': 'Calculations/TESTING/', "tasks": {'waitingtasks': ["Orca_Opt", "Orca_Dihedral", "Gromacs_Prod]"], "finnishedtasks": [], "failedtasks" : [], "runningtasks": []}}, 'test.json')
    thing = Job(JobDatabase.loadJobs('test.json')[0])
    print(thing.getNextTask()[0])
    thing.updateJob(runningtasks = thing.getNextTask()[0])
    # thing.id = 5


    JobDatabase.saveJob(thing.toDict(), 'test.json')
    # print(thing.__dict__)