from dataclasses import dataclass

@dataclass
class SBatchKEywords():
    head = "#!/usr/bin/zsh"
    jobName = "#SBATCH --job-name="
    outputFormat = "#SBATCH --output="
    time = "#SBATCH --time="  #### Request the time you need for execution. The full format is D-HH:MM:SS
    cpusPerTask = "#SBATCH --cpus-per-task=" # This corresponds to the number of processors (no hyperthreading possible)
    memPerCPU = "#SBATCH --mem-per-cpu="  #Limit for maximum memory per slot (in MB)
    nodes = "#SBATCH --nodes="
    tasks = "#SBATCH --ntasks="
    constraint = "#SBATCH --constraint="
    partition = "#SBATCH --partition="
    account = "#SBATCH --account="


class JobScripts(SBatchKEywords):
        def writeGausianJob(self, name, location,
                        inputFile = "combined.com",
                        outputFormat="GAUSSIANJOB.%J",
                        time = "96:00:00",
                        cpusPerTask = 8,
                        memPerCpu = "2000M"):

                with open(f"{location}/run_job.sh", "w") as output:
                       output.writelines([
                              f"{self.head}\n\n",
                              "########## start of batch directives #######\n",
                              f"{self.jobName}ID {name}\n",
                              f"{self.outputFormat}{outputFormat}\n",
                              f"{self.time}{time}\n",
                              f"{self.cpusPerTask}{cpusPerTask}\n",
                              f"{self.memPerCPU}{memPerCpu}\n",
                              f"\n###### start of shell commands ######\n\n",
                              "module load Gaussian/16.C.01-AVX2\n\n",  # load the necessary module files
                              f"g16 < {inputFile} > gauss.log"  # execute the gaussian binary
                       ])

        def writeOrcaJob(self, name, location, 
                        inputFile = "orca.inp",
                        outputFormat="output.%J.txt",
                        time = "0-24:00:00",
                        nodes = 1,
                        tasks = 16,
                        constraint = "Rocky8",
                        memPerCpu = "3900",
                        partiton = "c18m",
                        account = "p0020506"):
        
                with open(f"{location}/run_job.sh", "w") as output:
                       output.writelines([
                                f"{self.head}\n\n",
                                "########## start of batch directives #######\n",
                                f"{self.jobName}ID {name}\n",
                                f"{self.outputFormat}{outputFormat}\n",
                                f"{self.time}{time}\n",
                                f"{self.nodes}{nodes}\n",
                                f"{self.tasks}{tasks}\n",
                                f"{self.constraint}{constraint}\n",
                                f"{self.memPerCPU}{memPerCpu}\n",
                                f"{self.partition}{partiton}\n",
                                f"{self.account}{account}",
                                f"\n###### start of shell commands ######\n\n",
                                "module load GCC/11.3.0 OpenMPI/4.1.4  ORCA/5.0.4\n\n",
                                f"/cvmfs/software.hpc.rwth.de/Linux/RH8/x86_64/intel/skylake_avx512/software/ORCA/5.0.4-gompi-2022a/bin/orca  {inputFile}\n"
                       ])


        def writeGromacsJob(self, name, location, 
                        inputFile = "prod.tpr",
                        outputFormat="output.%J.txt",
                        time = "0-24:00:00",
                        nodes = 1,
                        tasks = 16,
                        constraint = "Rocky8",
                        memPerCpu = "3900",
                        partiton = "c18m",
                        account = "p0020506"):
        
                with open(f"{location}/run_job.sh", "w") as output:
                       output.writelines([
                              f"{self.head}\n\n",
                              "########## start of batch directives #######\n",
                              f"{self.jobName}ID {name}\n",
                              f"{self.outputFormat}{outputFormat}\n",
                              f"{self.time}{time}\n",
                              f"{self.nodes}{nodes}\n",
                              f"{self.tasks}{tasks}\n",
                              f"{self.constraint}{constraint}\n",
                              f"{self.memPerCPU}{memPerCpu}\n",
                              f"{self.partition}{partiton}\n",
                              f"{self.account}{account}",
                              f"\n###### start of shell commands ######\n\n",
                              "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n\n",  # load the necessary module files
                              f"mpirun -np 1 gmx_mpi mdrun -s {inputFile} -v -deffnm prod -ntomp 1 -plumed plumed.dat -tableb table_d0.xvg"  # execute the gaussian binary
                       ])

if __name__ == "__main__":
    #writeGausianJob("TEST")
    script = JobScripts().writeGromacsJob(name = "Test", location="")
    


