## Description

This code solves the Schrodinger equation for a set of interacting two-level quantum spins in a specified topology. More specifically, this branch contains modifications to allow the use of a genetic algorithm to optimise the couplings between such qubits to tailor the network for a specific purpose, such as fast Perfect State Transfer (PST) or performing rudimentary quantum gates.

See https://github.com/estaremp/spinchain/wiki for more info on the original code.

## Compilation

The program requires some MPI implementation, such as mpich, along with Python3, tk, matplotlib, numpy. On Debian based distributions:
```bash
sudo apt install mpich python3 python3-matplotlib python3-numpy python3-tk make
```

To compile this software use the make command in the directory containing the makefile, creating a "bin" folder containing binaries:
```bash
make
```

## Usage

To use the resulting software type the following command to view the command line help:

```bash
bin/spinnet --help
```

An example which opens the visualiser for a network, allowing interactive modification of the genome:

```bash
bin/spinnet -v "<A|E>AB500BC500CD500DE500#0000"
```

An example which uses the genetic algorithm to optimise the couplings for speed:

```bash
bin/spinnet -o "<A|E>AB500BC500CD500DE500#0000"
```

An example which just uses the original code to evaluate the network:

```bash
bin/spinnet -o "<A|E>AB500BC500CD500DE500#0000"
```

Results are placed in a folder called "output-latest" as well as a backup folder labelled with the date/time of the run. The output includes graphs/data files depending on the run parameters, generally the main file needed is "genetic.out" which contains the fitness and evolution of the run.

## Interesting Genomes

| Genome                                                                          |   Fidelity    | Transfer Time * J<sub>max</sub> |  Description   |
| ------------------------------------------------------------------------------- | :-----------: | :-------------------: | ------------- |
| "<A\|G>AB6364BC8216CD9000DE9000EF8216FG6364#000000"                             | 99.9% | 5.6 | 7-site linear PST chain, theoretically 100% fidelity |
| "<G\|B>AB7487AC9569CD9952DE9898EF9525FG7464#088888"                             | 99.7% | 5.4 | slightly faster than the 7-site PST chain at the expense of a small amount of fidelity |
| "<A\|G>AB9924BC9363CD9820DE9895EF9309FG9972FH9251<br>DH9943DI9978BI9289#0E2E206262" | 99.7% | 3.8 | significantly faster than the 7-site linear PST chain at the expense of two extra non-linear nodes |
| "@12.4<R+A+RA\|S+F-SF>AB2446BC4482CD8587D<br>E6111EF2608BG1270CH0187GH0288HI9016<br>DI0174IJ9763EJ0573GK8936KL0052HL6187<br>LM9969IM0278MN7059JN9980NO0601OP5783MP0141KQ1592<br>QR2524OS2541QT4791LT0023<br>PT8683#00000CC00C0CC0C0C0CC8CC800C8" | 94.9% | 12.4 | a two-qubit controlled phase gate using a 4x4 grid network |


