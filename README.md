This code solves the Schrodinger equation of a set of interacting qubits in a specified topology. More specifically this branch contains tools to use a genetic algorithm to optimise the couplings between such qubits to tailor the network for a specifc purpose, such as fast Perfect State Transfer (PST) or performing rudimentary quantum gates.

See https://github.com/estaremp/spinchain/wiki for more info on the original code.

# Compilation

To compile this software simply use the make command in the directory containing the makefile:

```bash
make
```

# Usage

To use the resulting software type the following command to view the command line help:

```bash
spinnet --help
```

An example which opens the visualiser for a network:

```bash
spinnet -v "<A|D>AB500BC500CD500#000"
```

An example which uses the genetic algorithm to optimise the couplings for speed:

```bash
spinnet -o "<A|D>AB500BC500CD500#000"
```


