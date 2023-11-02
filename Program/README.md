# How to do it easier!

I have written a small program to perform most of the tasks described in on the main page. How to run, and use it is decribed here.

**DISCLAIMER:**
**THIS HAS ONLY BEEN TESTED UNDER LINUX PROCEED WITH CAUTION!!**

## Table of content
  1. [Installation](#installation)
  1. [Calculation Preparation](#preparing-a-calculation)
  1. [Running the script](#running-the-script)
    1. [Program flow](#program-flow)
    2. [Things to do manually](#things-you-still-have-to-do-manually)
  1. [Last words of advice](#last-words-of-advice)

### Installation
To begin with, copy the *Program*-Folder into a location you like.
Then we create a virtual environment and install the requiered packages by executing the `ìnit.sh` file.
You can activate the environment by typing `source .venv/bin/activate` into a console.

### Preparing a Calculation 
To conduct a calculation you need to create a new folder in the `Input`-folder.
This folder contains the following files:
- The azobenzene in xyz-file format
- The fragment which is added later on (named fragment.xyz)
- The input file for this script (named Input)

<details>
  <summary><b>molecule.xyz</b></summary>

```
24
Coordinates for Azobenzene
  C   0.31295725236628     -0.71263358933773      1.03429269133599
  C   1.12031582619380      0.35777213222635      0.69146369118406
  C   0.64983145836037      1.31312303197833     -0.20500370415785
  C   -0.63561587373163      1.21008740182182     -0.72373278187689
  C   -1.42203081807197      0.11608624580387     -0.40496009309985
  C   -0.95099136724937     -0.84741151475525      0.47605733097021
  H   2.10893800095033      0.46050684024495      1.12071218613814
  H   -0.99908014721191      1.99109457672673     -1.38131582976800
  C   3.27057377682094      1.30099621562246     -1.15631499521902
  C   2.79030242883628      0.35946523145230     -2.06205098060979
  C   4.55294366702228      1.17666058705542     -0.63448512120610
  C   3.58457824881301     -0.71875956473376     -2.41111871846065
  H   1.80452446603223      0.47889552531611     -2.49357345277709
  C   5.32574808771716      0.07489356639787     -0.95962875730840
  H   4.92445211924435      1.94728817845603      0.03080069966782
  C   4.84473072024485     -0.87502324151455     -1.85001176084028
  H   3.21454071135506     -1.44633262706876     -3.12408651295957
  H   6.31520516427904     -0.03254322686432     -0.53085813626680
  H   5.45741915920509     -1.72677111431232     -2.12018994560903
  N   1.37527254118309      2.50117790140879     -0.50324627022334
  N   2.55983483559053      2.49555240321097     -0.84893480618743
  H   0.67554022397524     -1.45100405531037      1.73998276516609
  H   -2.41425884479749      0.02548639635230     -0.83124505514160
  H   -1.57423163712762     -1.69310730017724      0.74114755724938
```
</details>     

<details>
  <summary><b>fragment.xyz</b></summary>

```
31
PEDOT fragment to add
 C         1.440396        4.801444        9.924212
 H         0.420934        4.924465        9.505257
 H         1.373304        4.366064       10.942320
 C         2.217904        3.834611        9.089238
 C         3.413521        4.037653        8.432923
 H         4.004199        4.941150        8.387568
 N         1.854891        2.551740        8.908151
 N         2.801747        1.978917        8.183927
 N         3.744844        2.855785        7.887933
 C         4.968839        2.528574        7.166572
 H         5.823638        2.772959        7.830162
 H         5.008622        1.438863        6.959251
 C         5.144086        3.288605        5.843910
 H         5.094391        4.386146        5.996163
 H         6.172697        3.046066        5.502711
 C         4.205446        2.839638        4.711996
 O         4.669889        3.350338        3.504406
 H         4.247698        1.731141        4.637550
 C         2.718527        3.244275        4.888015
 H         2.256388        2.798423        5.784165
 H         2.666479        4.349419        4.978505
 O         1.992593        2.858912        3.765562
 C         2.537858        2.985267        2.534218
 C         3.897087        3.248111        2.401686
 C         4.331914        3.378119        1.096605
 H         5.342398        3.579343        0.773771
 S         2.988533        3.131043        0.108561
 C         1.860201        2.872014        1.334916
 H         0.807930        2.667487        1.209033
 O         2.001264        5.886481        9.980112
 H         2.420493        6.697648       10.021880
```
</details>     


<details>
  <summary><b>Input</b></summary>

```
Name=Azobenzene    #Atom index starts at 1!!! (e.g 1,2,3,.. NOT 0,1,2,...)
removeAtomNr=24             #This is the H-atom at the bonding site
combineAtomAt = 6           #This is the atom which forms a new Bond with the fragment
removeAtomFragmentNr=31     #This is the H-atom at the bonding site
combineFragmentAt = 30      #This is the atom which forms a new bond with the starting molecule
```
</details>

|        Keyword       | Description |
|:--------------------:|------------:|
|         Name         |  Name of the calculation  |
|     removeAtomNr     |  Index of the atom from the starting molecule to replace |
|     combineAtomAt    |  Index of the atom at which a new bond to the fragment is formed|
| removeAtomFragmentNr |  Index of the atom from the fragment to replace |
|   combineFragmentAt  |  Index to the atom at wich a new bond to the starting molecule is formed|

### Running the script
To run the script just execute the file `main.py` by typing `python3.10 main.py` into your console.
This will create a new file `schedule.xlsx` which contains information on all of your run calculations.

#### Program flow
The program contains instructions to perform most of the steps detailed in the main document.
The following steps are implemented:
1. Performing a orca optimization
2. Perform the orca dihedral scans
3. Combine the molecule and the fragment and optimize the structure
4. Calculate the force-field using Amber
5. All the mentioned Gromacs steps

The current state of the calculation can be tracked using the `schedule.xlsx`-file.
The first row `pandas-index` can be ignored as it is used internally for bookkeeping. The different numbers below the tasks correspond to the following code:
| Number Code |     Meaning    |
|:-----------:|:--------------:|
|      0      |    Not ready   |
|      1      |    Finnished   |
|      2      |     Running    |
|      3      | Ready to start |
|      -1     |      Error     |

**ONLY EDIT THIS FILE IF YOU KNOW WHAT YOU ARE DOING**

### Things you still have to do manually
This program has the capability to run mostly autonomous but there where a few things i was not able to automate.

#### Check if the combined molecule is correct
When starting with the Gaussian calculation have a quick look at the `combined.com` file. Sometimes my combination algorithm combines the in a weird way. As such it is your responsibility that the azobenzene unit points upwards away from the sulfure backbone.

#### Change CHANGEME.PDB to SHIFTED.PDB
The program will warn you about this one. Here you will have to place the molecule inside the given unitcell. Just open the file using VMD, then type `pbc box` into the console and move the molecule by pressing `8`. You can rotate the molcule by pressing `shift`. The tutorial on the main paige also applies. To save the new structure, select it in the vmd-MAIN window, then selct `save coordinates`.
*This should also be possible by performing some vector math, but i can not be bothered right now.*

### Last words of advice
Sometimes calculations will fail, i have not jet encountered every error message there is. This means that sometimes a calculation crashes and my program does not realize. Keep an eye on your running jobs using `squeue -u abc123456` and cross referencing the IDs.
To restart a calculation, either just rerun it yourself unsing `sbatch run_job.sh`. Or you can set the corresponding task to `3` in the `schedule.xlsx` file and delete/rename the folder containing the old calculation.

**Good Luck**