import h5py
import numpy as np
import json

if __name__ == '__main__':

    f1 = h5py.File("kPoint1QE.hdf5")
    #f2 = h5py.File("kPoint1Gut.hdf5")

    print(f1["PWCoeffs"])


    # for key in f1.keys():
    #     print(f1[key], f2[key])
    #     for S, G in zip(f1[key], f2[key]):
    #         if (S-G).sum() != 0:
    #             print("Hier stimmt was nicht!!")
    #             print(f"{S} - {G} != 0")


    