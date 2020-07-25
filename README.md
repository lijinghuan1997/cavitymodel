  # Cavity Model
  Get the self-consistent profiles of magnetic cavities with circular cross-section.
  # Description
  This magnetic cavity model is an equilibrium solution to the Vlasov-Maxwell equations, in the cylindrical coordinates, to reconstruct the magnetic cavities with circular cross-section. In the model, the proton/electron distribution is taken as functions of the particle motion invariants, particle energy, canonical angular momentum, and μ (magnetic moment) are included. They are combined to construct the current-carrying and background populations, for both electrons and protons, respectively. Then, the distribution functions are integrated to derive particle and current densities, which are substituted into Maxwell equations, to get the profiles of electromagnetic fields. Finally, you can get the self-consistent profiles for the magnetic cavity, including the electromagnetic fields and the particle distributions. 

  We notice that, the μ invariant in the distribution function, covers the drift velocity of the guiding center. Therefore, a uniform magnetic field and zero electric field are offered to construct this invariant, (with Vd=0 at first), before running the iteration part. Then taking advantage of the iteration procedure, we should update the profiles until they converge.

  Although we are not good at coding, the source code is very easy to use. Just run the code_submit.m in Matlab, you can get the profiles of electromagnetic field data saved. Then the handles for calculating the momentums, number density, bulk velocity, and pressure tensor included, are saved in momentum.m and model.m. Just run them in order, you can get all the related information about the magnetic cavity saved in momentum_end.mat and model_end.mat, respectively. These saved files have one-to-one correspondence of variable name and quantity, distributed with the distance from the cavity center.
  # Requirements
  This code should be compatible with Windows, Mac, and Linux operating systems, with matlab installed.

  The version for my computer is matlab (R2016a), and my operaing system is Win10. 

  The main arithmetic is ODE45 to solve the differential equations, so the requirement for the versions is very low. 
  # Instructions
  Run the code_submit.m, at first. 

  Run the momentum_submit.m, secondly.

  Run the model_submit.m, at last.
  # Running time
  The expected result lasts about an hour (code_submit.m), with the recommended specs (16 GB RAM, 8 cores@3.6 GHz)

  # Results
  code_submit.m----data_end.mat  
  % source code, generating the profiles
  5*4 cell (data), including: 
  5 lines: first run + 4 iterations;
  4 columns: the raw data drawn from the ode functions; x (distance from the cavity center); E (electric field); B (magnetic field)

  momentum_submit.m----momentum_end.mat  
  calculating the bulk velocity and pressure;
  12*2 cell (momentum);
  variable name; quantity

  model_submit.m----model_end.mat 
  details, including the number density, velocity, pressure, Ohm's law terms and so on. 
  26*2 cell (line_data);
  variable name; quantity

  The expected results are saved in the data folder. 
  # Parameter settings
  The first section in the code_submit.m is about this part. 

  The parameters are mainly about number density, temperature for different population, the magnetic strength at the center, as well as the angular velocity. 

  The normalized variables and the relations among the parameters, can refer to [Shustov et al., 2016].
  Shustov, P. I. , Artemyev, A. V. , Vasko, I. Y. , & Yushkov, E. V. . (2016). Kinetic models of sub-ion cylindrical magnetic hole. Physics of Plasmas, 23(12), 122903.

  You can change them to get single scale or cross-scale magnetic cavities.
  # License
  This code is covered under the MIT License.
