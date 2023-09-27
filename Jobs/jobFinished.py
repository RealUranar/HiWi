import sys
import glob

class HasFinished():
    def orca(folder):
        try:
            with open(glob.glob(f"{folder}/output.*.txt")[0], "r") as file:
                tail = file.readlines()[-4:]
        except IndexError:
            print("File not found!")
            return

        hasFinished = "TOTAL RUN TIME:" in str(tail)
        succesfull = "****ORCA TERMINATED NORMALLY****" in str(tail)
        return hasFinished, succesfull
    

    def orcaDihedral(folder):
        subfolders = ["singlet_left", "singlet_right", "triplet_left", "triplet_right"]
        for subfolder in subfolders:
            if HasFinished.orca(f"{folder}/{subfolder}/") != (True,True):
                return False
        return True


    def gaussian(folder):
        try:
            with open(glob.glob(f"{folder}/*.log")[0], "r") as file:
                tail = file.readlines()[-4:]
        except IndexError:
            print("File not found!")
            return

        hasFinished = "Normal termination of Gaussian" in str(tail)
        #succesfull = "****ORCA TERMINATED NORMALLY****" in str(tail)
        return hasFinished#, succesfull


if __name__ == "__main__":
    #print(HasFinished.orcaDihedral("../Jobs"))
    print(HasFinished.gaussian("../Jobs"))