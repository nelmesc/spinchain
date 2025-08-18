## Description

This code solves the Schrodinger equation for a set of interacting two-level quantum spins in a specified topology. More specifically, this branch contains modifications to allow the use of a genetic algorithm to optimise the couplings between such qubits to tailor the network for a specific purpose, such as fast Perfect State Transfer (PST) or performing rudimentary quantum gates.

See https://github.com/estaremp/spinchain/wiki for more info on the original code.

![Animation showcasing the method](https://github.com/lumorti/spinchain/raw/genetic/figures/combined.gif "Animation showcasing the method")

## Compilation

The program requires some MPI implementation, such as mpich, along with Python3, tk, matplotlib, PIL and numpy. 

On Debian based distributions these can be easily installed using:

```bash
sudo apt install mpich gfortran python3 liblapack-dev libblas-dev python3-tk libopenblas-dev python3-matplotlib python3-numpy python3-pil make
```

To compile this software use the make command in the directory containing the makefile, creating a "bin" folder containing binaries:

```bash
make
```

## Usage

For specific help with the resulting software type the following command to view the command line help:

```bash
bin/spinnet --help
```

An example which opens the visualiser for a network, allowing interactive modification of the genome:

```bash
bin/spinnet -v "<A|E>AB500BC500CD500DE500#0000"
```

An example which uses the genetic algorithm to optimise the couplings using the default settings, also generating an animation of the results:

```bash
bin/spinnet -o -a "<A|E>AB500BC500CD500DE500#0000"
```

An example which uses the genetic algorithm to optimise ONLY the on-site energies, keeping the coupling the same within the time window:

```bash
bin/spinnet -o --energy "<A|E>AB500BC500CD500DE500#0000"
```

An example which just uses the original code to evaluate the network, without optimising:

```bash
bin/spinnet "<A|E>AB500BC500CD500DE500#0000"
```

An example which adds a central uniformly coupled region of 10 spins with maximum nearest and next-nearest coupling:

```bash
bin/spinnet "<A|F>AB500BC500...10=9999,9999,...DE500EF500"
```
An example which demonstrates near-perfect state transfer for a 20-site NNN chain with J1 = 1 and J2 = 0.5:

```bash
bin/spinnet "<A|M>AA8889AB5347BC9999AC2999BD1763CD9999CE3712DE9989DF4097...20=9999,5000,...GJ4097HJ9989HK3712JK9999JL1763KM2999KL9999LM5347MM8889"
```
The structure for the NNN chain is similar to the NN optimisation, except the couplings at the edges are only optimised, and the genomic structure is as follows:

```bash
bin/spinnet "<A|F>AA9999AB9999AC9999...N=J1,J2,...DE9999EF9999FF9999"
```
where the repeated AA or FF terms represent the on-site energies (to be optimised or removed otherwise). 

Results are placed in a folder called "output-latest" as well as a backup folder labelled with the date/time of the run. The output includes graphs/data files depending on the run parameters, generally the main file needed is "genetic.out" which contains the fitness and evolution of the run.

## Interesting Genomes

| Genome                                                                          |   Fidelity    | Transfer Time * J<sub>max</sub> |  Description   |
| ------------------------------------------------------------------------------- | :-----------: | :-------------------: | ------------- |
| "<A\|G>AB6364BC8216CD9000DE9000EF8216FG6364#000000"                             | 99.9% | 5.6 | 7-site linear PST chain, theoretically 100% fidelity |
| "<G\|B>AB7487AC9569CD9952DE9898EF9525FG7464#088888"                             | 99.7% | 5.4 | slightly faster than the 7-site PST chain at the expense of a small amount of fidelity |
| "<A\|G>AB9924BC9363CD9820DE9895EF9309FG9972FH9251<br>DH9943DI9978BI9289#0E2E206262" | 99.7% | 3.8 | significantly faster than the 7-site linear PST chain at the expense of two extra non-linear nodes |
| "@12.4<R+A+RA\|S+F-SF>AB2416BC4591CD8793DE6017<br>EF2600BG0013CH0050GH0246HI9520DI0022IJ9646EJ0446GK8888KL0222<br>HL5897LM9999IM0517MN7651JN9918NO0129OP5618MP0043KQ1048QR2468<br>OS2543QT4939LT0104PT8801#00000CC00C0CC0C0C0CC8CC800C8" | 99.8% | 12.4 | a two-qubit controlled phase gate using a 4x4 grid network |


