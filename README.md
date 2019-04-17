# lattice_expansion
Python script originally written by Vishikh Athavale (05/2017). It calculates the lattice expansion of any (noncubic) crystal at a given temperature. It necessitates at least 10 phonon calculations (corresponding to small displacements of the lattice parameters) as input. Please have a look at the README for additional information.

usage: quasihNewApproach.py [-h] [-f [FOLDER [FOLDER ...]]] [-o OUTPUT] [-a A]
                            [-b B] [-c C] [-t [TEMPERATURE [TEMPERATURE ...]]]
                            [-at] [-plt]

Gives the best-fit unit-cell parameters at a given temperature, based on input data of unit-cell parameter from geometry relaxation, and corresponding free-energy of the unit-cell from phonopy.
   
Input file:  The first three columns of the input file should be the lattice-parameters a, b, c, and fourth column should be the free energy for that configuration computed using phonopy.

optional arguments:
  -h, --help            show this help message and exit
  -f [FOLDER [FOLDER ...]], --folder [FOLDER [FOLDER ...]]
                        List of the phonopy data folders (default is '../paracetamol/f1')
  -o OUTPUT, --output OUTPUT
                        Name of the output folder (default is current directory)
  -a A                  Initial guess for first unit-cell parameter "a"
  -b B                  Initial guess for second unit-cell parameter "b"
  -c C                  Initial guess for third unit-cell parameter "c"
  -t [TEMPERATURE [TEMPERATURE ...]], --temperature [TEMPERATURE [TEMPERATURE ...]]
                        Specify the list of temperatures anywhere in  0-300 K in
                        multiples of 10 K. Default is 300 K
                        (use -at for all temperatures)
  -at, --allt           Specify whether calculation is to be done at all
                        temperatures (don't use -t if you're using -at)
  -plt, --plot          Specify whether plots of a,b,c and V have to be made


- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                         -----------
                         |IMPORTANT|
                         -----------
The script will print out the 3 optimized lattice parameters x0, y0 and z0.
It will also output the hessian and the corresponding eigenvalues.
These eigenvalues must be positive! If not, the solution found does not correspond to a minimum.
In general, if you observe instabilities, try changing the initial guesses for the lattice parameters. There is usually a large range of initial guesses that will yield the same (stable) solution.
