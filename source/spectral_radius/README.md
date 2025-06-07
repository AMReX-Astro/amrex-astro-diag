# `spectral_radius`

This tool returns the thermodynamic state information about for
the zone where the reaction network Jacobian has the largest
spectral radius (maximum absolute value eigenvalue).  This can be
used to help assess the stiffness of the network for a particular
problems.

To build, do:

```
make DIM=2
```

changing the `DIM` line to match the dimension of your plotfile.

It is also important that the network you build with matches
the one used for generating the plotfile.  This is set via
the `NETWORK_DIR` parameter in the `GNUmakefile`.

Runtime parameters are managed by AMReX's ParmParse.  To run,
you specify the plotfile via `diag.plotfile`, either in an inputs
file or on the command line, e.g.:

```
./sprad2d.gnu.ex diag.plotfile=plt00000
```



