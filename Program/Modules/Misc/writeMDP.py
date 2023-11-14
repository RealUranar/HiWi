JOB_PARAMETERS = {
        "energy_minimization": {
        "title": "Energy Minimization",
        "define": "-DPOSRES",
        "integrator": "steep",
        "emtol": 100.0,
        "emstep": 0.01,
        "nsteps": 500000,
        "nstlist": 15,
        "cutoff-scheme": "Verlet",
        "ns_type": "grid",
        "coulombtype": "PME",
        "rcoulomb": 1.0,
        "rvdw": 1.0,
        "pbc": "xyz",
        "nstxout": 100,
        "nstvout": 1000,
        "nstenergy": 1000,
        "nstlog": 1000
    },
    "equilibration":  {
        "title": "equilibration NVT",
        "define": "-DPOSRES",
        "integrator": "md",
        "nsteps": 400000,
        "dt": 0.001,
        "nstxout": 1000,
        "nstvout": 1000,
        "nstenergy": 1000,
        "nstlog": 1000,
        "cutoff-scheme": "Verlet",
        "nstlist": 20,
        "ns_type": "grid",
        "pbc": "xyz",
        "rlist": 0.5,
        "coulombtype": "PME",
        "pme_order": 4,
        "fourierspacing": 0.16,
        "ewald_rtol": 1e-05,
        "tcoupl": "V-rescale",
        "tc-grps": "system",
        "tau_t": 0.1,
        "ref_t": 310,
        "Pcoupl": "no",
        "gen_vel": "yes",
        "gen_temp": 310,
        "gen_seed": -1,
        "constraint_algorithm": "lincs",
        "constraints": "h-bonds",
        "lincs_iter": 1,
        "lincs_order": 4,
        "annealing": "single",
        "annealing_npoints": 2,
        "annealing_time": "0 400",
        "annealing_temp": "0 310"
    },
    "production": {
        "title": "Production",
        "define": "-DPOSRES  -DPOTENTIAL",

        "integrator": "md", #Run parameters
        "nsteps": 40000000,
        "dt": 0.002,

        "nstxout": 10000, # Output control
        "nstvout": 10000,
        "nstenergy": 2000,
        "nstlog": 10000,

        "nstxtcout" : 1000,
        "xtc-precision" : 1000,  #Output frequency and precision for xtc file
        "xtc-grps"  : "F1",
        "energygrps"  :"F1",

        "continuation": "no", #Bond parameters 
        "constraint_algorithm": "lincs",
        "constraints": "h-bonds",
        "lincs_iter": 1,
        "lincs_order": 4,

        "cutoff-scheme": "Verlet", #Neighborsearching parameters
        "ns_type": "grid",   #search neighboring grid cells
        "nstlist": 10,      #10 fs, largely irrelevant with Verlet scheme
        "rcoulomb": 0.9,    #short-range electrostatic cutoff (in nm)
        "rvdw": 0.9,        #short-range van der Waals cutoff (in nm)

        "coulombtype": "PME",  # Particle Mesh Ewald for long-range electrostatics
        "pme_order": 4,        #cubic interpolation
        "fourierspacing": 0.16, #grid spacing for FFT
        "rlist"           : 1.0, #cut-off distance for the short-range neighbor list

        #Temperature coupling is on
        "tcoupl": "v-rescale",  # modified Berendsen thermostat
        "tc-grps": "System",
        "tau_t": 0.1,  #time constant, in ps
        "ref_t": 310,  #reference temperature, one for each group, in K

        "pbc"             :"xyz",           # 3-D PBC

        "DispCorr"        : "EnerPres",      #account for cut-off vdW scheme

        "gen_vel": "yes", #Velocity generation is on
        "gen_seed"       : -1 
    }
}

def writeMdpFile(job_type, temp=None, nsteps=None):
    if job_type not in JOB_PARAMETERS:
        print("Invalid job type")
        return

    parameters = JOB_PARAMETERS[job_type]
    if nsteps is not None:
        parameters["nsteps"] = nsteps
    if temp is not None:
        if "ref_t" in parameters:
            parameters["ref_t"] = temp
        if "gen_temp" in parameters:
            parameters["gen_temp"] = temp

    out = ""
    for key, value in parameters.items():
        out += f"{key} = {value}\n"
    return out

