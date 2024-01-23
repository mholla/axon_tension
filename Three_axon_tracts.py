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
    Totaltime = 4.0 
    Maxincnum = 1000 
    Defaultstabilization = 0.0002 
    Defaultdampingratio = 0.05 
    Mininc = 0.001 
    Incrementsize = 0.025
    Steppara = [Totaltime,Maxincnum,Defaultstabilization,Defaultdampingratio,Mininc,Incrementsize]

    # ======================================================
    # axon tract parabola shape parameters
    a_coeff = 1./30 #/mm
    b_coeff = -7 #mm
    m_coeff = 0 #mm

    # axon tract numbers
    curve_num1 = 1
    curve_num2 = 2
    curve_num3 = 3

    # parameters that control the segments
    Geometric_length = 0.130 #mm
    Stretch_ratio = 2
    InfluenceRadius = 1.0 #mm
    Axon_tract_property = [Geometric_length,Stretch_ratio]

    # ======================================================
    # vary axon tract stiffness within a range  
    num = 10
    min_stiffness = 10 #N/m
    max_stiffness = 1000 #N/m
    Stiffness_range = np.linspace(min_stiffness,max_stiffness,num)

    # ======================================================
    # naming of model parts
    PartName = 'Part-1'
    Step = 'Step-1'
    InstanceName = 'Part-1-1'
    UMAT = "./UMAT_axon_tension.f"

    for i in range(0,len(Stiffness_range)):
        for j in range(0,len(Stiffness_range)):

            ModelName =  'Model-primary%d-secondary%d' %(i,j) 
            JobName = 'Job-primary%d-secondary%d' %(i,j)

            Stiffness_primary = Stiffness_range[i]
            Stiffness_secondary = Stiffness_range[j]

            Create_Bilayered_Rectangle(ModelName, PartName, Dimensions)
            Create_Material(ModelName,Materials)
            Create_Section(ModelName, PartName,Dimensions)
            Create_Assembly(ModelName,InstanceName)
            Create_Sets(ModelName,InstanceName,Dimensions)
            Create_Step(ModelName, Step, Steppara)
            Create_Contact(ModelName, Step)
            Create_Boundary_Conditions(ModelName, Step)
            Create_Mesh(ModelName, PartName, Dimensions)
            Create_Axon_Connection(a_coeff,b_coeff,m_coeff,curve_num1,ModelName,InstanceName,Axon_tract_property,Stiffness_primary,InfluenceRadius,Dimensions)
            Create_Axon_Connection(a_coeff,b_coeff,m_coeff,curve_num2,ModelName,InstanceName,Axon_tract_property,Stiffness_secondary,InfluenceRadius,Dimensions)
            Create_Axon_Connection(a_coeff,b_coeff,m_coeff,curve_num3,ModelName,InstanceName,Axon_tract_property,Stiffness_secondary,InfluenceRadius,Dimensions)
            Create_mesh_node_sets(ModelName,InstanceName)
            Modify_input(ModelName)
            Modify_input_for_initialize_growth_variable(ModelName)
            Create_Job(ModelName, JobName, UMAT)

    # ==============================================================
    # specify parameters for post-processing 
    # The benchmark is selected as the maximum primary axon stiffness and the minimum secondary axon stiffness, which is the case of Job-primary9-secondary0.odb
    ODB_Name_Bench = 'Job-primary%d-secondary%d.odb' %(9,0)
    NodesetName = 'TOPSURF_NODES'
    Variable = 'COORD'
    [x_coords_bench, y_coords_bench] = Post_processing_odbs(ODB_Name_Bench,NodesetName,Step,Variable)
    ODBtype = 2
    parameter1_range = Stiffness_range
    parameter2_range = Stiffness_range
    psi_array = Calculate_psi(ODBtype,NodesetName,Step,Variable,parameter1_range,parameter2_range,x_coords_bench,y_coords_bench)
    psi_array = np.flip(psi_array,axis=0)
    
    # ==============================================================
    # write data to the .csv file for plotting
    np.savetxt("psi_array_three_curves.csv", psi_array, delimiter=",")