import numpy as np
import warnings
warnings.filterwarnings("error")

with open("Data_CO_VASP/WAVECAR", "rb") as file:
    file.seek(0,0)
    zeroCount = 0
    numCount = 0
    bitesRead = 0
    while bitesFile > bitesRead:
        
        num = np.fromfile(file, dtype = np.float64, count = 1)

        if len(num) == 0:
            print(num)
            break

        if np.abs(num) < 0.0000001 and np.abs(num) != 0:
            file.seek(file.tell()-8,0)
            num = np.fromfile(file, dtype = np.complex64, count = 1)

        bitesRead = file.tell()


        if num == 0:
            if numCount != 0:
                print(f"Bite {file.tell()-(numCount*8+8)} = {numCount} : Nums")
                numCount = 0

            zeroCount += 1

        if num != 0:
            if zeroCount != 0:
                print(f"Bite {file.tell()-(zeroCount*8+8)} = {zeroCount} : Zeros")
                zeroCount = 0
            numCount += 1

        bitesRead += 8




#from subprocess import call

# Path ="Data_NaCl_QE"
# call(f"#!/bin/bash\n \
#     tar -zcvf ZIP.tar.gz {Path}/NaCl.save {Path}/NaCl.xml",shell=True)
