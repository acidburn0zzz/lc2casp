# lc2casp

A translator for logic programs with constraint atoms to CASP.

## Usage

The basic usage to ground, translate, and solve logic programs with constraint atoms is
```bash
gringo PROGRAM | lc2casp | clingcon --mode=clasp
```

## Installation

The translator depends on the clingo sources and has to be configured accordingly.
To generate a basic configuration, run
```bash
make FLAGS
```
and adjust the generated FLAGS file accordingly.

After configuration, the default make target builds the translator.
There is also a test target to run a small suite of acceptance tests.
