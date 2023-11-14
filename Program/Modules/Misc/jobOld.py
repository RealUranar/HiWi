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
    def __init__(self,name:str, id:int, location:str, tasks:dict):
        self.name = name
        self.id  = id
        self.location = location
        self.tasksPandas = tasks
        self._sortTasks(tasks)

    def _sortTasks(self, tasksPandas):
        taskToProgram = {
            "Orca_Opt": Orca_opt,
            "Orca_Dihedral": Orca_Dihedral,
            "Gaussian" : Gaussian_opt,
            "Amber" : Amber,
            "GromacsEnergy": GromacsEnergy,
            "GromacsEquil" : GromacsEquill,
            "GromacsProduction": GromacsProd}
        
        self.tasks = {
            "Done" : [],
            "Ready" : [],
            "Running" : [],
            "Error" : [],
        }
        for task in tasksPandas.keys():
            taskClass = taskToProgram[task](self)

            match tasksPandas[task]:
                case 1:
                    self.tasks["Done"].append(taskClass)
                case 2:
                    self.tasks["Running"].append(taskClass)
                case 3:
                    self.tasks["Ready"].append(taskClass)
                case -1:
                    self.tasks["Error"].append(taskClass)

    def getNextTask(self) -> list[Task]:
        return self.tasks["Ready"]
    
    def getFinishedTask(self) -> list[Task]:
        return self.tasks["Done"]
    
    def getRunningTask(self) -> list[Task]:
        return self.tasks["Running"]

    def getPanda(self) -> dict:
        return list(self.tasksPandas.values())

    def getLocation(self):
        return self.location

    def updateJob(self, *args,**kwargs):
        """
        
        """
        for newTask in kwargs.keys():
            if newTask not in list(self.tasksPandas.keys()):
                print("Invaid task!")
                continue
            self.tasksPandas[newTask] = kwargs[newTask]

        self._sortTasks(self.tasksPandas)

if __name__ == "__main__":
    thing = Job(name = "TEST", id=2, location="Here/" ,tasks= dict({
                                "Orca_Opt": 3,
                                "Orca_Dihedral": 0,
                                "Gaussian" : 0,
                                "Amber" : 0,
                                "GromacsEnergy":0,
                                "GromacsEquil" : 0,
                                "GromacsProduction": 0}))
    

    newThing = thing.getNextJob()[0]
    print(newThing.newPath)