class Job():
    def __init__(self,name:str, id:int, location:str, tasks:dict):
        self.name = name
        self.id  = id
        self.location = location
        self.tasksPandas = tasks
        self._sortTasks(tasks)

    def _sortTasks(self, tasksPandas):
        self.tasks = {
            "Done" : [],
            "Ready" : [],
            "Running" : [],
            "Error" : [],
        }
        for task in tasksPandas.keys():
            match tasksPandas[task]:
                case 1:
                    self.tasks["Done"].append(task)
                case 2:
                    self.tasks["Running"].append(task)
                case 3:
                    self.tasks["Ready"].append(task)
                case -1:
                    self.tasks["Error"].append(task)

    def getNextJob(self):
        return self.tasks["Ready"]
    
    def getFinishedJobs(self):
        return self.tasks["Done"]
    
    def getRunningJobs(self):
        return self.tasks["Running"]

    def getPanda(self) -> dict:
        return list(self.tasksPandas.values())

    def getLocation(self):
        return self.location

    def updateJob(self, *args,**kwargs):
        """
        
        """
        print(kwargs)
        for newTask in kwargs.keys():
            if newTask not in list(self.tasksPandas.keys()):
                print("Invaid task!")
                continue
            self.tasksPandas[newTask] = kwargs[newTask]

        self._sortTasks(self.tasksPandas)
