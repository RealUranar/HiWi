import sys
import glob

class HasFinished():
    def _readTail(folder, gauss=False):
        tail = ""
        try:
            if gauss:
                folder = glob.glob(f"{folder}*.log")[0]
            else:
                folder = glob.glob(f"{folder}output.*.txt")[0]

        except IndexError:
            print(f"File {folder}output.*.txt not found!")
            return tail
        
        with open(folder, "r") as file:
            tail = str(file.readlines()[-15:])
        return tail

    def orca(folder):
        hasFinished, succesfull = False, False
        tail = HasFinished._readTail(folder)
        hasFinished = "TOTAL RUN TIME:" in tail
        succesfull = "****ORCA TERMINATED NORMALLY****" in tail
        return hasFinished, succesfull
    

    def orcaDihedral(folder):
        subfolders = ["singlet_left", "singlet_right", "triplet_left", "triplet_right"]
        for subfolder in subfolders:
            if HasFinished.orca(f"{folder}/{subfolder}/") != (True,True):
                return False
        return True


    def gaussian(folder):
        hasFinished, succesfull = False, True
        tail = HasFinished._readTail(folder, gauss=True)
        hasFinished = "Normal termination of Gaussian" in tail
        #succesfull = "****ORCA TERMINATED NORMALLY****" in str(tail)
        return hasFinished, succesfull


    def gromacs(folder):
        hasFinished, succesfull = False, False

        tail = HasFinished._readTail(folder)
        if "Segmentation fault" in str(tail) or "DUE TO TIME LIMIT" in str(tail):
            hasFinished = True
        elif "Writing final coordinates." in str(tail): 
            hasFinished = True
            succesfull = True

        return hasFinished, succesfull


if __name__ == "__main__":
    #print(HasFinished.orcaDihedral("../Jobs"))
    print(HasFinished.orca("../Jobs/"))