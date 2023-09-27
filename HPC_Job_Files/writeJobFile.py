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


def writeGausianJob(name, location ,
                    inputFile = "orca.com",
                    outputFormat="GAUSSIANJOB.%J",
                    time = "96:00:00",
                    cpusPerTask = 8,
                    memPerCpu = "2000M"):
    
    with open(f"{location}run_job.sh", "w") as output:
        output.write(f"{SBatchKEywords.head}\n\n")
        output.write("########## start of batch directives #######\n")
        output.write(f"{SBatchKEywords.jobName}{name}\n\
{SBatchKEywords.outputFormat}{outputFormat}\n\
{SBatchKEywords.time}{time}\n\
{SBatchKEywords.cpusPerTask}{cpusPerTask}\n\
{SBatchKEywords.memPerCPU}{memPerCpu}\n")
        output.write(f"\n###### start of shell commands ######\n\n")
        output.write("module load Gaussian/16.C.01-AVX2\n\n")  # load the necessary module files
        output.write(f"g16 < {inputFile} > gauss.log")  # execute the gaussian binary

def writeOrcaJob(name, location, 
                inputFile = "orca.inp",
                outputFormat="output.%J.txt",
                time = "0-24:00:00",
                nodes = 1,
                tasks = 16,
                constraint = "Rocky8",
                memPerCpu = "3900",
                partiton = "c18m",
                account = "p0020506"):
    
    with open(f"{location}run_job.sh", "w") as output:
        output.write(f"{SBatchKEywords.head}\n\n")
        output.write("########## start of batch directives #######\n")
        output.write(f"{SBatchKEywords.jobName}{name}\n\
{SBatchKEywords.outputFormat}{outputFormat}\n\
{SBatchKEywords.time}{time}\n\
{SBatchKEywords.nodes}{nodes}\n\
{SBatchKEywords.tasks}{tasks}\n\
{SBatchKEywords.constraint}{constraint}\n\
{SBatchKEywords.memPerCPU}{memPerCpu}\n\
{SBatchKEywords.partition}{partiton}\n\
{SBatchKEywords.account}{account}")
    
        output.write(f"\n###### start of shell commands ######\n\n")
        output.write("module load GCC/11.3.0 OpenMPI/4.1.4  ORCA/5.0.4\n\n")  # load the necessary module files
        output.write(f"/cvmfs/software.hpc.rwth.de/Linux/RH8/x86_64/intel/skylake_avx512/software/ORCA/5.0.4-gompi-2022a/bin/orca  {inputFile}")  # execute the gaussian binary



if __name__ == "__main__":
    #writeGausianJob("TEST")
    writeOrcaJob("Test")
    


