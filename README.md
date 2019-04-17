# lattice_expansion
Python script originally written by Vishikh Athavale (05/2017). It calculates the lattice expansion of any (noncubic) crystal at a given temperature. It necessitates at least 10 phonon calculations (corresponding to small displacements of the lattice parameters) as input.

usage: quasihNewApproach.py [-h] [-f [FOLDER [FOLDER ...]]] [-o OUTPUT] [-a A]
                            [-b B] [-c C] [-t [TEMPERATURE [TEMPERATURE ...]]]
                            [-at] [-plt]

Gives the best-fit unit-cell parameters at a given temperature, based on input data of unit-cell parameter from geometry relaxation, and corresponding free-energy of the unit-cell from phonopy.
   

Please type `python quasihNewApproach.py --help` for more detail.

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                         -----------
                         |IMPORTANT|
                         -----------
The script will print out the 3 optimized lattice parameters x0, y0 and z0.
It will also output the hessian and the corresponding eigenvalues.
These eigenvalues must be positive! If not, the solution found does not correspond to a minimum.
In general, if you observe instabilities, try changing the initial guesses for the lattice parameters. There is usually a large range of initial guesses that will yield the same (stable) solution.
