import Python_subroutine_axon_tension
import numpy as np

# preload all subroutines to the local python environment
execfile("Python_subroutine_axon_tension.py")

if __name__ == '__main__':

    # ======================================================
    # dimensions of the model
    Length = 80.0 #mm
    Height = 40.0 #mm
    Cortex_thickness = 2.0 #mm
    Dimensions = [Length,Height,Cortex_thickness]

    # ======================================================
    # material properties
    mu_cortex = 100.0 #kPa
    lame_cortex= 9.3*mu_cortex #kPa
    stiffness_ratio = 3 
    mu_subcortex = stiffness_ratio*mu_cortex #kPa
    lame_subcortex = 9.3*mu_subcortex #kPa
    growth_rate = 0.05 
    Materials = [mu_cortex,lame_cortex,mu_subcortex,lame_subcortex,growth_rate]

    # ======================================================
    # step timing parameters
    Totaltime = 4.6 
    Maxincnum = 1000 
    Defaultstabilization = 0.0002 
    Defaultdampingratio = 0.05 
    Mininc = 0.001 
    Incrementsize = 0.025
    Steppara = [Totaltime,Maxincnum,Defaultstabilization,Defaultdampingratio,Mininc,Incrementsize]

    # ======================================================
    # axon tract parabola shape parameters
    num = 10 
    tangent_min = 0.5 
    tangent_max = 5 
    span_min = 16 #mm
    span_max = 25 #mm
    [a_coeffs_array,b_coeffs_array] = Compute_parabola_coeffs(num,tangent_min,tangent_max,span_min,span_max)
    m_coeff = 0 #mm

    # only primary axon tract
    curve_num = 1 

    # axon tract stiffness  
    Axon_tract_stiffness = 150 #N/m
    InfluenceRadius = 1.0 #mm

    # parameters that control the segments
    Geometric_length = 0.130 #mm
    Stretch_ratio = 2 
    Axon_tract_property = [Geometric_length,Stretch_ratio]

    # ======================================================
    # naming of model parts
    PartName = 'Part-1'
    Step = 'Step-1'
    InstanceName = 'Part-1-1'
    UMAT = "./UMAT_axon_tension.f"

    for i in range(0,num):
        for j in range(0,num):

            ModelName = 'Model-span%d-tangent%d' %(i,j) 
            JobName = 'Job-span%d-tangent%d' %(i,j)

            a_coeff = a_coeffs_array[i,j]
            b_coeff = b_coeffs_array[i,j]

            # This is a special case that the axon tract is close to the bottom of the model, where the biased mesh is too sparse to apply the coupling constrain. 
            # We increase the influence radius to apply the coupling and it doesn't affect the folding morphology
            if ((i == 9) and (j == 9)):
                InfluenceRadius = 1.5

            Create_Bilayered_Rectangle(ModelName, PartName, Dimensions)
            Create_Material(ModelName,Materials)
            Create_Section(ModelName, PartName, Dimensions)
            Create_Assembly(ModelName,InstanceName)
            Create_Sets(ModelName,InstanceName,Dimensions)
            Create_Step(ModelName, Step, Steppara)
            Create_Contact(ModelName, Step)
            Create_Boundary_Conditions(ModelName, Step)
            Create_Mesh(ModelName, PartName, Dimensions)
            Create_Axon_Connection(a_coeff,b_coeff,m_coeff,curve_num,ModelName,InstanceName,Axon_tract_property,Axon_tract_stiffness,InfluenceRadius,Dimensions)
            Create_mesh_node_sets(ModelName,InstanceName)
            Modify_input(ModelName)
            Modify_input_for_initialize_growth_variable(ModelName)
            Create_Job(ModelName, JobName, UMAT)

    # ==============================================================
    # specify parameters for post-processing 
    # The benchmark is selected as the minimum axon tract span and the minimum axon tract tangent at one root, which is the case of Job-span0-tangent0.odb
    ODB_Name_Bench = 'Job-span%d-tangent%d.odb' %(0,0)
    NodesetName = 'TOPSURF_NODES'
    Variable = 'COORD'
    [x_coords_bench, y_coords_bench] = Post_processing_odbs(ODB_Name_Bench,NodesetName,Step,Variable)
    ODBtype = 0
    parameter1_range = a_coeffs_array[:,0]
    parameter2_range = b_coeffs_array[:,0]
    psi_array = Calculate_psi(ODBtype,NodesetName,Step,Variable,parameter1_range,parameter2_range,x_coords_bench,y_coords_bench)
    psi_array = np.flip(psi_array,axis=0)
    
    # ==============================================================
    # write data to the .csv file for plotting
    np.savetxt("psi_array_geometry.csv", psi_array, delimiter=",")