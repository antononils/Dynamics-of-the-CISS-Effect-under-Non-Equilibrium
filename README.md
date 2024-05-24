# Dynamics of the CISS Effect under Non-Equilibrium

This repository contains the code used in the Bachelor thesis "Dynamics of the Chirality Induced Spin
Selectivity Effect under Non-Equilibrium", conducted at Uppsala university as a part of the program of Engineering Physics.


## Analytical Approximation
This folder contains the code used for the analytical approximation:
  * [CISS_analytical](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Analytical%20Approximation/CISS_analytical.m)

    _In this file the geometrical and energy parameters are chosen. It also contains the code for plotting the results._
  * [Perturbation_analytic](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Analytical%20Approximation/Perturbation_analytic.m)

    _This file contains the code that sets up the Hamiltonian, both the unperturbed Hamiltonian H<sub>0</sub> and the perturbation term V(t)._
  * [Probability_density](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Analytical%20Approximation/Probability_density.m)

    _An analytic solution to the first perturbation term in the Dyson series using Symbolic Matlab is implemented in this file._


## Numerical Approximation
This folder contains the code used for the numerical approximation:
  * [CISS](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Numerical%20Approximation/CISS.m)
    
    _This file contains the code that runs the numerical simulation of the CISS effect._
  * [ConvergenceTest](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Numerical%20Approximation/ConvergenceTest.m)
    
    _This file contains the code that tests the convergence of the order of perturbation._
  * [Distrubutions](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Numerical%20Approximation/Distributions.m)
    
    _This file contains the code that calculates the probability densities and the spin polarization._
  * [Hamiltonian](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Numerical%20Approximation/Hamiltonian.m)
    
    _This file contains the code that creates the unperturbed Hamiltonian given specified parameters._
  * [Perturbation](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Numerical%20Approximation/Perturbation.m)
    
    _This file contains the code that creates the potential matrix for different perturbations._
  * [Wavefunction](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Numerical%20Approximation/Wavefunction.m)
    
    _This file contains the code that calculates the Dyson Series expansion up to order n for a given initial condition._
  * [WavefunctionWithErrors](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Numerical%20Approximation/WavefunctionWithErrors.m)
    
    _This file contains the code that calculates the total error from integration in "Wavefunction"._

## Plotting 
This folder contains the code used for plotting the results:
  * [ColorPlot](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Plotting/ColorPlot.m)
    
    _This file contains the code that creates the color plots for the probability densities and the spin polarization._
  * [PolarizedPlotting](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Plotting/PolarizedPlotting.m)
  
    _This file contains the code that creates the color plots for the polarized probability density._
  * [CISSAnimation](https://github.com/antononils/Dynamics-of-the-CISS-Effect-under-Non-Equilibrium/blob/main/Plotting/CISSAnimation.m)
  
    _This file contains the code that creates animations of probability densities for each spin over time._


## Authors 
This project was done by Anton O'Nils, Christoffer Damsgaard Falck, Gustav Teglund, and Hannes Tjulin.
