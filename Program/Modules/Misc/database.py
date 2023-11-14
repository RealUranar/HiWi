import json
import os
from job import Job

class JobDatabase:
    @classmethod
    def saveJob(cls, job_new, filepath):
        existing_jobs = cls.loadJobs(filepath)
        for job_old in existing_jobs:
            if job_old["id"] == job_new["id"]:
                existing_jobs.remove(job_old)
                continue

        existing_jobs.insert(job_new["id"]-1, job_new)

        with open(filepath, 'w') as f:
            json.dump(existing_jobs, f, indent=4)

    @classmethod
    def loadJobs(cls, filepath):
        if not os.path.exists(filepath):
            return []
        
        with open(filepath, 'r') as f:
            if os.stat(filepath).st_size == 0:
                return []
            jobs = json.load(f)

        return jobs
    
    @classmethod
    def getLastID(cls, filepath):
        jobs = cls.loadJobs(filepath)
        if len(jobs) == 0:
            return 0
        return jobs[-1]["id"]
    

if __name__ == "__main__":
    # JobDatabase.saveJob({"id": 1,'name': 'test', 'location': 'Calculations/TESTING/', "tasks": {'waitingtasks': ["Orca_Opt", "Orca_Dihedral", "GromacsProduction"], "finnishedtasks": [], "failedtasks" : [], "runningtasks": []}}, 'test.json')
    # JobDatabase.saveJob({"id": 2,'name': 'test2', 'location': 'Calculations/HalloTEST/', "tasks": {'waitingtasks': ["Orca_Opt", "Orca_Dihedral", "GromacsProduction"], "finnishedtasks": [], "failedtasks" : [], "runningtasks": []}}, 'test.json')

    job = Job(JobDatabase.loadJobs('test.json')[0])
    # job.updateJob(runningtasks=[job.getNextTask()[0]])
    # job.id = 1

    JobDatabase.saveJob(job.toDict(), 'test.json')
