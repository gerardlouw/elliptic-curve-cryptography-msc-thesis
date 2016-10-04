# Elliptic curve cryptography (MSc thesis)

## Dependencies

For compilation of the LaTeX source, you need a full TeX Live installation, as well as SageMath (used to generate some plots). To run the Sage code, you just need SageMath.

These packages can be installed on Ubuntu 16.04 as follows.

```
sudo apt install texlive-full sagemath
```

## Running

To compile the LaTeX, just run `make -C latex`, which should generate `latex/ecc.pdf`. Run `make -C clean` to remove all generated artifacts. A precompiled `ecc.pdf` has been included in the repository root for convenience.

To do a simple test of the Sage code, run `make -C sage`. For your own tests, just load the relevant files in a Sage environment (see the `sage/test.sage` file for an example).
