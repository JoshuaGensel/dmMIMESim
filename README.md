# dmMIMESim

First round MIME:
1. Simulate sequence variants in pool of a wild type sequence (random mutagenesis 1)
2. Randomly assign fitness values for each position-wise mutation and pairwise epistatic effects
3. Derive species frequencies in equilibrium
4. Add statistical noise
5. Output positionwise mutant counts, sequences and the groundtruth

Second round MIME:

6. Read in results from first round of MIME
7. Simulate sequence variants in unbound pool of first round experiment (random mutagenesis 2)
8. Derive species frequencies in equilibrium
9. Add statistical noise
10. Output positionwise mutant counts, sequences and the groundtruth

*Detailed explanation will be added*

## Install

Go into the MIMESim directory and create a build directory and cmake and make from there.

```
cd /path/to/MIMESim
mkdir build
cd build
cmake ..
make
```
This will create the binary file *.../MIMESim/build/bin/MIMESim_prog*

## Run


Running the program with an optional parameter, giving the path to the result directory. 
If no parameter is given, the output will be written to the newly created directory *../results*.

```
/path/to/MIMESim_prog [result_directory1] # first round simulation
/path/to/MIMESim_prog [result_directory2] [result_directory1] # second round simulation
```

If no parameter file is present in the given result directory (or no directory was given at all), the programm will run the simulation with default paraters *(tba)*.
After the run it will create the parameter file in the result directory together with the other output. 

If you want to set the simulation parameters, create the *parameters.txt* file in the desired result directory and add the following parameter names and values (tab separated):

| parameter          | type          | description  |
| :-------------     |:-------------:| :------------|
| kd_wt              | (float)       |   (absolute Kd for the wildtype sequence (default: 1))|
| p_effect           | (float)       |   probability per position to have an effect on function when mutated (default: 0.5) |
| p_epistasis        | (float)       |   probability for each pair of mutations to have epistatic effect (default: 0.3) |
| L                  | (int)         |   Length of the sequence (default: 50) |
| q                  | (int)         |   (number of possible symbols per position (default: 2))  |
| M                  | (int)         |   (size of initial population (default: 12000000)  |
| p_mut              | (float)       |   probability for each position of being mutated (default: 0.01) |
| p_error            | (float)       |   probability of introducing an error (default: 0.001) |
| seed               | (int)         |   seed for random number generator used |
| epistasis_mut_excl | (bool)        |   when true, each position can be in an epistatic interaction with at most one other position (default: false) |
| B_tot              | (float)       |   total amount of protein in competition experiment, given in relation to the amount of sequences (default: 2.0) |
| max_mut            | (int)         |   maximal number of mutation allowed per sequence (default=-1, no max) |

**Comments:**
1. Right now, kd_wt is always 1, but this can be easily adjusted in the code.
2. It is only working for q=2 atm.


## Output

* 1d counts (with noise)
* 2d counts (with noise)
* full-length sequences (with noise)
* species (sequence variants without sequencing noise)
* amount of free protein in equilibrium, per binding experiment
* ground truth Kds (position and pairwise) and epistasis values

```
result_directory1
├── 1d
│   ├── 1.txt
│   ├── 2.txt
│   ├── 3.txt
│   └── 4.txt
├── 2d
│   ├── 1.txt
│   ├── 2.txt
│   ├── 3.txt
│   └── 4.txt
├── free_protein_1_2.txt
├── free_protein_3_4.txt
├── pairwise_epistasis.txt
├── pairwise_kds.txt
├── parameters.txt
├── sequences
│   ├── 1.txt
│   ├── 2.txt
│   ├── 3.txt
│   └── 4.txt
├── single_kds.txt
└── species
    ├── 1.txt
    ├── 2.txt
    ├── 3.txt
    └── 4.txt
```

```
results_directory2
├── 1d
│   ├── 5.txt
│   ├── 6.txt
│   ├── 7.txt
│   └── 8.txt
├── 2d
│   ├── 5.txt
│   ├── 6.txt
│   ├── 7.txt
│   └── 8.txt
├── free_protein_5_6.txt
├── free_protein_7_8.txt
├── pairwise_epistasis.txt
├── pairwise_kds.txt
├── parameters.txt
├── sequences
│   ├── 5.txt
│   ├── 6.txt
│   ├── 7.txt
│   └── 8.txt
├── single_kds.txt
└── species
    ├── 5.txt
    ├── 6.txt
    ├── 7.txt
    └── 8.txt
```
### Sequence pool encoding
* 1 = wt_bound,
* 2 = wt_unbound,
* 3 = mut_bound,
* 4 = mut_unbound,
* 5 = mut_bound_bound,
* 6 = mut_bound_unbound,
* 7 = mut_unbound_bound,
* 8 = mut_unbound_unbound,
