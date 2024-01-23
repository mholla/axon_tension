# axon_tension
The readme document for codes to reproduce the results of the axon tension paper.

To run the simulation, place the Fortran file ‘UMAT_axon_tension.f’ and the Python file ‘Python_subroutine_axon_tension.py’ in the same folder. 

Next, place one of the Python scripts for a specific section in the same folder.

For example, to run one of the Python scripts, open the terminal and run the following command:

**abaqus cae noGUI=Parabola_geometry.py**

For other scripts:

**abaqus cae noGUI=Thickness_perturbation.py**

**abaqus cae noGUI=Three_axon_tracts.py**

**abaqus cae noGUI=Wiring_length_change_eta.py**

**abaqus cae noGUI=Wiring_length_keep_eta.py**


To plot post-processed heatmaps, place the data files (.cvs or .npy) with the plotting script in the same folder. The data for the corresponding plotting script is shown below:
Data file  | Plotting script 
------------- | -------------
data_critical-strain.csv  | plot_critical_strain.py 
data_geometrical.csv  | plot_geometrical_parametric.py
data_single-curve.csv | plot_geometrical_parametric.py
data_three_axon_tracts.csv  | plot_three_curves_parametric.py
Job-K1-100-K2-25-total.npy <br> Job-K1-100-K2-100-total.npy<br> Job-K1-100-K2-400-total.npy | plot_wiring_length_change_eta.py
Job-K1-100-K2-100-total.npy <br> Job-K1-200-K2-200-total.npy<br> Job-K1-300-K2-300-total.npy | plot_wiring_length_keep_eta.py


Then open the Python terminal and run the following command:

**python plot_critical_strain.py** 

**python plot_geometrical_parametric.py**

**python plot_thickness_perturbation_parametric.py**

**python plot_three_curves_parametric.py**

**python plot_wiring_length_change_eta.py**

**python plot_wiring_length_keep_eta.py**
