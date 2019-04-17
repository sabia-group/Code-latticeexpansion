#This python code finds the optimal lattice parameters for a (noncubic) crystal
#Vishikh Athavale, 29/05/2017

from __future__ import print_function
import numpy as np
from scipy.optimize import curve_fit
import argparse, os, re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats import linregress
from matplotlib.font_manager import FontProperties

def makeFunc(angles):
    def func(X, e0, x0, y0, z0, hxx, hxy, hxz, hyy, hyz, hzz):
        """
        Computes the free energy of the unit-cell 
        using a second order taylor polynomial 
        expanded about the unit-cell parameter x0, y0, z0.
        F = e0 + [x-x0 y-y0 z-z0].H.[x-x0; y-y0; z-z0],
        where H is the hessian matrix

        Args:
          X: Array with equilibrium lattice-parameters a,b,c at which free energy
          has to be computed
          x0: First unit-cell parameter length about which Taylor series is expanded
          y0: Second unit-cell parameter length about which Taylor series is expanded
          z0: Third unit-cell parameter length about which Taylor series is expanded
          hxx: Element of the hessian matrix
          hxy: Element of the hessian matrix (= hyx)
          hxz: Element of the hessian matrix (= hzx)
          hyy: Element of the hessian matrix
          hyz: Element of the hessian matrix (= hzy)
          hzz: Element of the hessian matrix
          
        Returns:
          An array containing the free energy of the unit-cell 
          for the given values of a,b,c.  
        """
        '''
        freeEnergy = np.zeros(len(X[0]))
        for i in range(0,len(X[0])):
                freeEnergy[i] = (e0 + 
                                np.dot(np.dot(np.subtract([X[0][i],X[1][i],X[2][i]],[x0,y0,z0]),
                                              [[hxx,hxy,hxz],[hxy,hyy,hyz],[hxz,hyz,hzz]]),
                                       np.subtract([[X[0][i]],[X[1][i]],[X[2][i]]],[[x0],[y0],[z0]])))
        return freeEnergy
        '''
        alpha = angles[0]
        beta = angles[1]
        gamma = angles[2]
        #a,b,c = X
        
        J = np.array([[1, np.cos(gamma), np.cos(beta)],
                      [0, np.sin(gamma), (np.cos(alpha) -
                                          np.cos(beta)*np.cos(gamma))/np.sin(gamma)],
                      [0, 0, 1/np.sin(gamma)*np.sqrt(-(np.cos(alpha))**2 -
                             (np.cos(beta))**2 +
                             2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) +
                             (np.sin(gamma))**2)]])
        H = np.array([[hxx,hxy,hxz],[hxy,hyy,hyz],[hxz,hyz,hzz]])
        freeEnergy = np.zeros(len(X[0]))
        for i in range(0,len(X[0])):
            freeEnergy[i] = (e0 +
                             np.dot(np.dot(np.subtract([X[0][i],X[1][i],X[2][i]],[x0,y0,z0]),
                                      np.dot(np.dot(np.transpose(J),H),J)),
                                    np.subtract([[X[0][i]],[X[1][i]],[X[2][i]]],[[x0],[y0],[z0]])))
        return freeEnergy
        '''
        freeEnergy = e0 + ((a - x0)*(hxx*(a - x0) + (b - y0)*(hxx*np.cos(gamma) + hxy*np.sin(gamma)) + 
        (c - z0)*(hxx*np.cos(beta) + hxy*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))*1/np.sin(gamma) + 
        hxz*1/np.sin(gamma)*np.sqrt(-np.cos(alpha)**2 - np.cos(beta)**2 + 2*np.cos(alpha)*np.cos(beta)*
        np.cos(gamma) + np.sin(gamma)**2))) + 
        (b - y0)*((a - x0)*(hxx*np.cos(gamma) + hxy*np.sin(gamma)) + 
        (b - y0)*(np.cos(gamma)*(hxx*np.cos(gamma) + hxy*np.sin(gamma)) + 
        np.sin(gamma)*(hxy*np.cos(gamma) + hyy*np.sin(gamma))) + 
        (c - z0)*(np.cos(gamma)*(hxx*np.cos(beta) + hxy*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))*
        1/np.sin(gamma) + hxz*1/np.sin(gamma)*np.sqrt(-np.cos(alpha)**2 - np.cos(beta)**2 + 
        2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) + np.sin(gamma)**2)) + 
        np.sin(gamma)*(hxy*np.cos(beta) + hyy*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))*
        1/np.sin(gamma) + hyz*1/np.sin(gamma)*np.sqrt(-np.cos(alpha)**2 - np.cos(beta)**2 + 
        2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) + np.sin(gamma)**2)))) + 
        (c - z0)*((a - x0)*(hxx*np.cos(beta) + hxy*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))*
        1/np.sin(gamma) + hxz*1/np.sin(gamma)*np.sqrt(-np.cos(alpha)**2 - np.cos(beta)**2 + 
        2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) + np.sin(gamma)**2)) + 
        (b - y0)*(np.cos(beta)*(hxx*np.cos(gamma) + hxy*np.sin(gamma)) + 
        (np.cos(alpha) - np.cos(beta)*np.cos(gamma))*1/np.sin(gamma)*(hxy*np.cos(gamma) + 
        hyy*np.sin(gamma)) + 1/np.sin(gamma)*(hxz*np.cos(gamma) + hyz*np.sin(gamma))*
        np.sqrt(-np.cos(alpha)**2 - np.cos(beta)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) + 
        np.sin(gamma)**2)) + (c - z0)*
        (np.cos(beta)*(hxx*np.cos(beta) + hxy*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))*1/np.sin(gamma) + 
        hxz*1/np.sin(gamma)*np.sqrt(-np.cos(alpha)**2 - np.cos(beta)**2 + 2*np.cos(alpha)*np.cos(beta)*
        np.cos(gamma) + np.sin(gamma)**2)) + (np.cos(alpha) - np.cos(beta)*np.cos(gamma))*
        1/np.sin(gamma)*(hxy*np.cos(beta) + hyy*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))*
        1/np.sin(gamma) + hyz*1/np.sin(gamma)*np.sqrt(-np.cos(alpha)**2 - np.cos(beta)**2 + 
        2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) + np.sin(gamma)**2)) + 
        1/np.sin(gamma)*np.sqrt(-np.cos(alpha)**2 - np.cos(beta)**2 + 2*np.cos(alpha)*np.cos(beta)*
        np.cos(gamma) + np.sin(gamma)**2)*(hxz*np.cos(beta) + 
        hyz*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))*1/np.sin(gamma) + 
        hzz*1/np.sin(gamma)*np.sqrt(-np.cos(alpha)**2 - np.cos(beta)**2 + 2*np.cos(alpha)*np.cos(beta)*
        np.cos(gamma) + np.sin(gamma)**2))))) 
        
        #print(freeEnergy)
        return(freeEnergy)
        '''
        #x, y, z = X
        #return e0 +  (x - x0)* (hxx * (x - x0) + hxy *  (y - y0) + hxz *  (z - z0)) + (y - y0) * (hxy * (x - x0) + hyy * (y - y0) + hyz * (z - z0)) + (hxz * (x - x0) +  hyz * (y - y0) + hzz * (z - z0)) * (z - z0) 

    return func

def fitData(data, angles, aGuess, bGuess, cGuess, temperature):
    """
    Gives the best-fit unit-cell parameters at
    a given temperature, based on input data of unit-cell 
    parameter and free-enrgy.

    Args:
      data: Numpy array containg unit-cell parameters a,b,c and the free energy from
            getData function.
      angles: Array with angles alpha, beta, gamma
      aGuess:Intitial guess for first unit-cell parameter 
      bGuess:Intitial guess for second unit-cell parameter 
      cGuess:Intitial guess for third unit-cell parameter 
      temperature: Temperature 
    
    Returns:
      Array with best fit parameters e0, x0, y0, z0, hxx, hxy, hxz, hyy, hyz, hzz
    """
    a = data[:,0]
    b = data[:,1]
    c = data[:,2]
    gibbsFreeEnergy = data[:,3] - min(data[:,3]) 
    np.set_printoptions(suppress=True)
    #print(gibbsFreeEnergy)
    # min(data[:,3] is the free energy offset to avoid numerical errors
    # that may arise because free energy values can be big numbers 

    # p0 is the initial guess for e0, x0, y0, z0, hxx, hxy, hxz, hyy, hyz, hzz
    p0 = 0., aGuess, bGuess, cGuess, 5., 5., 5., 5., 5., 5.

    fit,cov = curve_fit(makeFunc(angles), (a,b,c), gibbsFreeEnergy, p0)
    #print(cov)
    return fit
    

def getData(folders, temperature):
    """
    Collates all the data from phonopy output folder and forms an array which
    can then be used as the input for the fitting of the data.
    
    Args:
      folder: Name of the phonopy output folder 
      temperature: Temperature at which the data is extracted
    
    Returns:
      Numpy array with unit-cell parameters and the Gibbs free energy
    """
    
    data = []
    #os.chdir('/home/athavale/Desktop/molecularCrystals/paracetamol/f1')
    for folder in folders: 
        for folderName, subfolders, filenames in os.walk(folder):
            for filename in filenames:
                #print('File inside'+folderName+': '+filename)
                if(filename=='phonopy-FHI-aims-free_energy.dat'):
                    #print('File inside'+folderName+': '+filename)
                   
                    # Get electronic energy from FHI-Aims output file
                    #fhiOutputFile = open(folderName+'/'+'relax.out','r')
                    fhiOutputFile = open(folderName+'/'+'out_relax','r') # paracetamol II
                    result = fhiOutputFile.read()
                    fhiOutputFile.close()
                    electronicEnergy = float(re.findall('Total energy uncorrected      :         (.*) eV'
                                                  ,result)[-1])
                    
                    # Get vibrational free energy from phonopy output
                    phonopyFile = folderName+'/'+filename
                    freeEnergyData = np.loadtxt(phonopyFile)
                    freeEnergy = freeEnergyData[temperature/10][1]

                    # Get unit-cell parameters a,b,c 
                    geometryFile = folderName+'/'+'geometry.in.next_step'
                    geometry = np.genfromtxt(geometryFile, dtype=float, skip_header=5,
                                             usecols = (1,2,3), max_rows=3)
                    a = np.linalg.norm(geometry[0])
                    b = np.linalg.norm(geometry[1])
                    c = np.linalg.norm(geometry[2])
                    alpha = np.arccos(np.dot(geometry[1],geometry[2])/
                                      (np.linalg.norm(geometry[1])*np.linalg.norm(geometry[2])))
                    beta = np.arccos(np.dot(geometry[0],geometry[2])/
                                      (np.linalg.norm(geometry[0])*np.linalg.norm(geometry[2])))
                    gamma = np.arccos(np.dot(geometry[1],geometry[0])/
                                      (np.linalg.norm(geometry[1])*np.linalg.norm(geometry[0])))
                    angles =np.array([alpha, beta, gamma])
                    
                    # Get pressure/volume contribution to Gibbs free energy
                    # 1 bar = 6.24151E-7 ev/Angstrom^3
                    PV = ((6.24151E-7)*(a*b*c*
                          np.sqrt(1-(np.cos(alpha))**2 - (np.cos(beta))**2
                                  -(np.cos(gamma))**2 + 
                                  2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))))
                    writeData = [a, b, c, electronicEnergy + freeEnergy + PV]
                    data.append(writeData)
                    
    #print(np.asanyarray(data))
    return np.asanyarray(data), angles



def printOutput(fit,angles,outputFolder,temperature,i):
    """
    Prints the output of best-fit parameters and writes them into a file
    
    Args:
      fit: Array contating the fitted parameters
      angles: Array with lattice angles
      outputFolder: Folder where the output is to be written
      temperature: Temperature at which fitting was done
      i: Gives the index of the finalResult array to write the necessary
        data--temperature, x0, y0, z0  which can be printed 
        for easy glance at the resultsi
    """

    output = open(outputFolder+'/fittingAt'+str(temperature)+'K.dat', 'w')
    output.write('e =' + str(fit[0]) + '\nx0 =' + str(fit[1]) + 
                 '\ny0 =' + str(fit[2]) + '\nz0 =' + str(fit[3]) + 
                 '\nhxx =' +str(fit[4]) + '\nhxy =' +str(fit[5]) +
                 '\nhxz =' +str(fit[6]) + '\nhyy =' +str(fit[7]) +
                 '\nhyz =' + str(fit[8]) + '\nhzz =' + str(fit[9]))
    output.close()
            
    # Print output
    print ('e =', fit[0], '\nx0 =', fit[1] , 
           '\ny0 =', fit[2], '\nz0 =', fit[3], 
           '\nhxx =',fit[4], '\nhxy =' ,fit[5],
           '\nhxz =',fit[6], '\nhyy =' ,fit[7],
           '\nhyz =', fit[8], '\nhzz =', fit[9])
 
    H = np.array([[fit[4], fit[5], fit[6]], [fit[5], fit[7], fit[8]],[fit[6],fit[8],fit[9]]])
    print('Hessian matrix:')
    print(H)
    print('EigenValues and eigenvectors:')
    print(np.linalg.eigh(H))
    x0 = fit[1]
    y0 = fit[2]
    z0 = fit[3]
    alpha = angles[0]
    beta = angles[1]
    gamma = angles[2]
    volume = (x0*y0*z0*np.sqrt(1-(np.cos(alpha))**2 - (np.cos(beta))**2
                                  -(np.cos(gamma))**2 + 
                                  2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)))
#    print (volume)
    finalResult[i] = [temperature, float("{0:.4f}".format(x0)),
                      float("{0:.4f}".format(y0)),
                      float("{0:.4f}".format(z0)),float("{0:.4f}".format(volume))]

def plot_data(finalResult):
    """
    Plots lattice parameters and volume against temperature

    Args:
      finalResult: Array containing temperature and corresponding lattice
                   parameters at that temperature. 
    """
    T = [row[0] for row in finalResult]
    a0 = [row[1] for row in finalResult]
    b0 = [row[2] for row in finalResult]
    c0 = [row[3] for row in finalResult]
    V = [row[4] for row in finalResult]

    plt.figure(1)
    plt.plot(T,a0,'C0.-')
    plt.ylabel('Lattice parameter $a_0$ ($\AA$)')
    plt.xlabel('Temperature (K)')
    plt.savefig('aVsTemp.pdf')
    plt.show()
    plt.close()
    
    plt.figure(2)
    plt.plot(T,b0,'C1.-')
    plt.ylabel('Lattice parameter $b_0$ ($\AA$)')
    plt.xlabel('Temperature (K)')
    plt.savefig('bVsTemp.pdf')
    plt.show()
    plt.close()
    
    plt.figure(3)
    plt.plot(T,c0,'C2.-')
    plt.ylabel('Lattice parameter $c_0$ ($\AA$)')
    plt.xlabel('Temperature (K)')
    plt.savefig('cVsTemp.pdf')
    plt.show()
    plt.close()

    plt.figure(4)
    plt.plot(T,V,'C3.-')
    plt.ylabel('Unit-cell volume ($\AA^3$)')
    plt.xlabel('Temperature (K)')
    plt.savefig('volumeVsTemp.pdf')
    plt.show()
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '''Gives the best-fit unit-cell parameters at
a given temperature, based on input data of unit-cell 
parameter from geometry relaxation, and corresponding 
free-energy of the unit-cell from phonopy.
   
Input file:  The first three columns of the input file 
should be the lattice-parameters a, b, c, 
and fourth column should be the free energy for that configuration
computed using phonopy              ''',formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-f", "--folder", type =str, default =
                        ['../paracetamol/f1/'], 
                        nargs = '*', help =
     "List of the phonopy data folders (default is '../paracetamol/f1')")
    parser.add_argument("-o","--output", type = str, default ='.', help = 
                        "Name of the output folder (default is current directory)")
    parser.add_argument("-a",type = float, default = 7.18, help = 
                        'Initial guess for first unit-cell parameter "a"')   
    parser.add_argument("-b",type = float, default = 9.30, help =
                        'Initial guess for second unit-cell parameter "b"')   
    parser.add_argument("-c",type = float, default = 11.00, help =
                        'Initial guess for third unit-cell parameter "c"')   
    parser.add_argument("-t","--temperature", type = int, nargs  = '*',
                        default = [300], help =
                        '''Specify the list of temperatures anywhere in  0-300 K in
multiples of 10 K. Default is 300 K
(use -at for all temperatures)''')
    parser.add_argument("-at", "--allt", action = 'store_true', help =
                       '''Specify whether calculation is to be done at all
temperatures (don't use -t if you're using -at)''')
    parser.add_argument("-plt", "--plot", action = 'store_true', help =
                        "Specify whether plots of a,b,c and V have to be made")
    args = parser.parse_args()
    
    i=0
    if(args.allt==False):
        finalResult = [[]]*len(args.temperature)
        for tempr in args.temperature:
            print('Calculations at T = %r K' %tempr)
            data, angles = getData(args.folder, tempr)
            fit = fitData(data, angles, args.a, args.b, args.c, tempr)
            printOutput(fit, angles, args.output, tempr, i)
            i = i + 1
    else:
        finalResult = [[]]*31 # 300 K
        for tempr in range(0,310,10): #300 K
            data, angles = getData(args.folder, tempr)
            fit = fitData(data, angles, args.a, args.b, args.c, tempr)
            print('T = %r K' %(tempr))
            printOutput(fit,angles, args.output, tempr, i)
            i = i+1
    #print("\nTemperature \t a0 \t b0 \t c0")
    #summary = open('resultSummary.dat','w')
    #summary.write("\nTemperature a0 \t b0 \t c0 \t Volume\n")
    #print(finalResult)
    #for data in finalResult:
    #    print(data)
        #summary.write("\n".join(str(elem)) for elem in data)   
        #summary.writelines(["%r\t" % item  for item in data])
        #summary.write("\n")

    #summary.close()
    #print("\nThe angles are:")
    #print(angles)

    #plot_data(finalResult)
