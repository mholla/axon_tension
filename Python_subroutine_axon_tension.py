from abaqus import *
from abaqusConstants import *
import __main__

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import numpy as np
import os
import sys
import time
import string
import random
import math
import bisect
import subprocess
import scipy
from scipy import interpolate
import sympy as sp
from odbAccess import *
import odbAccess


def Create_Bilayered_Rectangle(ModelName, PartName, Dimensions):

    """ Create bilayered rectangle 

    Parameters
    ----------
    ModelName : string
        name of the model to be modified
    PartName : string
        name of the part to be modified
    Dimensions : list of floats
        geometric parameters of the model: Length, Height, Cortex_thickness
    
    """

    cliCommand("""session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)""")

    Length = Dimensions[0]
    Height = Dimensions[1]
    Cortex_thickness = Dimensions[2]

    # Define the model and modeltype
    mdb.Model(name=ModelName, modelType=STANDARD_EXPLICIT)   
    s = mdb.models[ModelName].ConstrainedSketch(name='__profile__', sheetSize=100.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)

    # Sketch the rectangle as the simulation domain
    s.rectangle(point1=(-Length/2, 0.0), point2=(Length/2, -Height))
    p = mdb.models[ModelName].Part(name=PartName, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    p = mdb.models[ModelName].parts[PartName]
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()

    # Partition rectangle into two layers: the first layer as the cortex layer and the second layer as the subcortex layer
    f, e, d = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(Length/8, -Height/4, 0.0), normal=(0.0, 0.0, 1.0)), sketchPlaneSide=SIDE1, origin=(0.0, -Height/2, 0.0))
    s = mdb.models[ModelName].ConstrainedSketch(name='__profile__', sheetSize=100.0, gridSpacing=4, transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

    # Sketch the line to differentiate two layers
    s.Line(point1=(-Length/2, Height/2 - Cortex_thickness), point2=(Length/2, Height/2 - Cortex_thickness))
    s.HorizontalConstraint(entity=g.findAt((0.0, Height/2 - Cortex_thickness)), addUndoState=False)
    p = mdb.models[ModelName].parts[PartName]
    e1, d2 = p.edges, p.datums
    p.PartitionFaceBySketch(faces=p.faces.findAt(((Length/8, -Height/4, 0.0), )), sketch=s)
    s.unsetPrimaryObject()


def Create_Material(ModelName,Materials):

    """ Create materials

    Parameters
    ----------
    ModelName : string
        name of the model to be modified
    Materials : list of floats
        material parameters of the model: mu_cortex, lame_cortex), mu_subcortex, lame_subcortex, growth_rate

    """ 

    mu_cortex = Materials[0]
    lame_cortex = Materials[1]
    mu_subcortex = Materials[2]
    lame_subcortex = Materials[3]
    growth_rate = Materials[4]

    # Cortex (gray matter) shearmodulus_cortex, lameconstant_cortex, shearmodulus_subcortex, lameconstant_subcortex, Gctx
    mdb.models[ModelName].Material(name='CORTEX',description='*****************************************************************\n  Specification Of Material Properties\n*****************************************************************\n\n\nCOMMENTS FROM *USER MATERIAL\n============================\n\n*      shearmodulus of cortex = props(1)\n*      lame constant of cortex = props(2)\n*      shearmodulus of subcortex  = props(3)\n*      lame constant of subcortex  = props(4)\n*      growth rate  = props(5)')
    mdb.models[ModelName].materials['CORTEX'].UserMaterial(mechanicalConstants=(mu_cortex, lame_cortex, mu_subcortex, lame_subcortex, growth_rate))
    mdb.models[ModelName].materials['CORTEX'].Depvar(n=1)

    # Subcortex (white matter) shearmodulus_cortex, lameconstant_cortex, shearmodulus_subcortex, lameconstant_subcortex, Gctx
    mdb.models[ModelName].Material(name='SUBCORTEX',description='*****************************************************************\n  Specification Of Material Properties\n*****************************************************************\n\n\nCOMMENTS FROM *USER MATERIAL\n============================\n\n*      shearmodulus of cortex = props(1)\n*      lame constant of cortex = props(2)\n*      shearmodulus of subcortex  = props(3)\n*      lame constant of subcortex  = props(4)\n*      growth rate  = props(5)')
    mdb.models[ModelName].materials['SUBCORTEX'].UserMaterial(mechanicalConstants=(mu_cortex, lame_cortex, mu_subcortex, lame_subcortex, growth_rate))
    mdb.models[ModelName].materials['SUBCORTEX'].Depvar(n=1)


def Create_Section(ModelName, PartName, Dimensions):

    """ Create solid sections

    Parameters
    ----------
    ModelName : string
        name of the model to be modified
    PartName : string
        name of the part to be modified
    Dimensions : list of floats
        geometric parameters of the model: Length, Height, Cortex_thickness

    """ 

    Length = Dimensions[0]
    Height = Dimensions[1]

    # Create Sets and different Sections
    p = mdb.models[ModelName].parts[PartName]
    f = p.faces

    cortex_face = f.findAt(((Length/8, -Height/80, 0.0), ))
    p.Set(faces=cortex_face, name='Set-cortex')

    subcortex_face = f.findAt(((Length/8, -Height/4, 0.0), ))
    p.Set(faces=subcortex_face, name='Set-subcortex')

    # Define Sections
    mdb.models[ModelName].HomogeneousSolidSection(name='Section-1-cortex', material='CORTEX', thickness=None)
    mdb.models[ModelName].HomogeneousSolidSection(name='Section-2-subcortex', material='SUBCORTEX', thickness=None)

    # Assign Sections
    region = p.sets['Set-cortex']
    p.SectionAssignment(region=region, sectionName='Section-1-cortex', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    region = p.sets['Set-subcortex']
    p.SectionAssignment(region=region, sectionName='Section-2-subcortex', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)


def Create_Assembly(ModelName,InstanceName):

    """ Create assembly

    Parameters
    ----------
    ModelName : string
        name of the model to be modified
    InstanceName : string
        name of the instance to be modified
    Dimensions : list of floats
        geometric parameters of the model: Length, Height, Cortex_thickness

    """ 

    # Create datum system
    a = mdb.models[ModelName].rootAssembly
    p = mdb.models[ModelName].parts[PartName]

    a.DatumCsysByDefault(CARTESIAN)
    a.Instance(name=InstanceName, part=p, dependent=ON)


def Create_Sets(ModelName,InstanceName,Dimensions):

    """ Create sets

    Parameters
    ----------
    ModelName : string
        name of the model to be modified
    InstanceName : string
        name of the instance to be modified
    Dimensions : list of floats
        geometric parameters of the model: Length, Height, Cortex_thickness

    """ 

    Length = Dimensions[0]
    Height = Dimensions[1]
    Cortex_thickness = Dimensions[2]

    # Define edge Sets
    a = mdb.models[ModelName].rootAssembly
    e1 = a.instances[InstanceName].edges

    bottom_edge = e1.findAt(((-Length/4, -Height, 0.0), ))
    a.Set(edges=bottom_edge, name='bottom')

    left_edge = e1.findAt(((-Length/2, -Height/4, 0.0), ), ((-Length/2, -Height/32, 0.0), ), ((-Length/2, -Height/160, 0.0), ))
    a.Set(edges=left_edge, name='left_bc')

    right_edge = e1.findAt(((Length/2, -Height/1.5, 0.0), ), ((Length/2, -Height/32, 0.0), ), ((Length/2, -Height/160, 0.0), ))
    a.Set(edges=right_edge, name='right_bc')

    interface_edge = e1.findAt(((-Length/4, -Cortex_thickness, 0.0), ))
    a.Set(edges=interface_edge, name='interface')

    topsurface_edge = e1.findAt(((Length/4, 0.0, 0.0), ))
    a.Surface(side1Edges=topsurface_edge, name='Surf-top')


def Create_Step(ModelName, Step, Steppara):

    """ Create step

    Parameters
    ----------
    ModelName : string
        name of the model to be modified
    Step : string
        name of the step to be created
    Steppara : list of floats
        time parameters of the step: Totaltime, Maxincnum (max increment number), Defaultstabilization (default 
        stabilization magnitude),  Defaultdampingratio (default adaptive damping ratio), Mininc (the minimum increment 
        size), Incrementsize

    """

    Totaltime = Steppara[0]
    Maxincnum = Steppara[1]
    Defaultstabilization = Steppara[2]
    Defaultdampingratio = Steppara[3]
    Mininc = Steppara[4]
    Incrementsize = Steppara[5]

    mdb.models[ModelName].StaticStep(name='Step-1', previous='Initial', 
        timePeriod=Totaltime, maxNumInc=Maxincnum, stabilizationMagnitude=Defaultstabilization, 
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION, 
        continueDampingFactors=False, adaptiveDampingRatio=Defaultdampingratio, 
        initialInc=Incrementsize, minInc=Mininc, maxInc=Incrementsize, nlgeom=ON)


def Create_Contact(ModelName, Step):

    """ Create contacts

    Parameters
    ----------
    ModelName : string
        name of the model to be modified
    Step : string
        name of step to create contact in
    """ 

    # Create the interaction property 
    mdb.models[ModelName].ContactProperty('IntProp-1')
    mdb.models[ModelName].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((10.0, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    mdb.models[ModelName].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)
    mdb.models[ModelName].interactionProperties['IntProp-1'].GeometricProperties(
        contactArea=1.0, padThickness=None)

    # Apply the self-contact to the top surface of the cortex layer
    a = mdb.models[ModelName].rootAssembly
    region1=a.surfaces['Surf-top']
    region2=a.surfaces['Surf-top']
    mdb.models[ModelName].SurfaceToSurfaceContactStd(name='Int-1', 
        createStepName=Step, main=region1, secondary=region2, sliding=FINITE, 
        thickness=ON, interactionProperty='IntProp-1', adjustMethod=NONE, 
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None)


def Create_Boundary_Conditions(ModelName, Step):

    """ Create boundary conditions

    Parameters
    ----------
    ModelName : string 
        name of the model to be modified
    Step : string
        name of step to create BC in
    
    """ 

    a = mdb.models[ModelName].rootAssembly

    region = a.sets['bottom']
    mdb.models[ModelName].EncastreBC(name='BC_bottom', createStepName=Step, region=region, localCsys=None)

    region = a.sets['right_bc']
    mdb.models[ModelName].XsymmBC(name='BC_right', createStepName=Step, region=region, localCsys=None)

    region = a.sets['left_bc']
    mdb.models[ModelName].XsymmBC(name='BC_left', createStepName=Step, region=region, localCsys=None)


def Create_Mesh(ModelName, PartName, Dimensions):

    """ Create bilayered rectangle 

    Parameters
    ----------
    ModelName : string
        name of the model to be modified
    PartName : string
        name of the part to be modified
    Dimensions: list of floats
        geometric parameters of the model: Length, Height, Cortex_thickness

    Notes
    -----
    Changes in the mesh parameters or the dimensions will change hard-coded node numbers in node sets

    """

    Length = Dimensions[0]
    Height = Dimensions[1]

    p = mdb.models[ModelName].parts[PartName]
    domain = p.faces.findAt(((Length/6, -Height/3, 0.0), ), ((Length/6, -Height/30, 0.0), ), ((-Length/6, -Height/60, 0.0), ))
    p.setMeshControls(regions=domain, technique=STRUCTURED)

    elemType1 = mesh.ElemType(elemCode=CPE4H, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD)

    p.setElementType(regions=(domain, ), elemTypes=(elemType1, elemType2))

    top_edge = p.edges.findAt(((-Length/4, -Height/40, 0.0), ), ((Length/2, -Height/40, 0.0), ), ((Length/4, 0.0, 0.0), ), ((-Length/2, -Height/160, 0.0), ))
    p.seedEdgeBySize(edges=top_edge, size=0.5, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)

    interface_edge = p.edges.findAt(((-Length/4, -Height/40, 0.0), ), ((-Length/2, -Height/30, 0.0), ), ((-Length/4, -Height/20, 0.0), ), ((Length/2, -Height/25, 0.0), ))
    p.seedEdgeBySize(edges=interface_edge, size=0.5, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)

    left_edge = p.edges.findAt(((-Length/2, -Height/4, 0.0), ))
    right_edge = p.edges.findAt(((Length/2, -Height/4, 0.0), ))
    p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=left_edge, end2Edges=right_edge, minSize=0.5, maxSize=6.0, constraint=FINER)

    bottom_edge = p.edges.findAt(((-Length/4, -Height, 0.0), ))
    p.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=bottom_edge, minSize=0.4, maxSize=12.0, constraint=FINER)
    p.generateMesh()


def Create_Axon_Connection(a_coeff,b_coeff,m_coeff,curve_num,ModelName,InstanceName,Axon_tract_property,Axon_tract_stiffness,InfluenceRadius,Dimensions):

    """ Create Axon Connections

    Parameters
    ----------
    a_coeff, b_coeff, m_coeff : float
        parameters that describe the shape and the location of the parabola
    curve_num: number integer
        number of the wiring curve and wiring set
    ModelName : string
        name of the model to be modified
    InstanceName : string
        name of the instance to be modified
    Axon_tract_property : list of floats
        parameters of the axon tract: Geometric_length (undeformed length of line segments in the axon bundle) and 
        Stretch_ratio (pre-stretch of the axon bundle)
    Axon_tract_stiffness : float  
        stiffness of the primary/secondary axon bundle
    InfluenceRadius : float 
        influence radius of the coupling area
    Dimensions : list of floats
        geometric parameters of the model: Length, Height, Cortex_thickness
    
    """

    Geometric_length = Axon_tract_property[0]
    Stretch_ratio = Axon_tract_property[1]
    Length = Dimensions[0]

    # Generate original coordinates with the given parabola coefficients
    num = 200
    xmax = np.abs(np.sqrt((-2-b_coeff)/a_coeff))

    # curve_num = 1 means the primary axon tract on the middle
    # curve_num = 2 means the secondary axon tract on the right
    # curve_num = 3 means the secondary axon tract on the left
    if curve_num == 1:
        x_coords_parabola = np.linspace(-xmax,xmax,num)
        y_coords_parabola = (x_coords_parabola)**2*a_coeff+b_coeff
    elif curve_num == 2:
        x_coords_parabola = np.linspace(-xmax,xmax,num) + xmax
        y_coords_parabola = (x_coords_parabola - xmax)**2*a_coeff+b_coeff
    elif curve_num == 3:
        x_coords_parabola = np.linspace(-xmax,xmax,num) - xmax
        y_coords_parabola = (x_coords_parabola + xmax)**2*a_coeff+b_coeff

    # Phase shift parabola m_coeff distance to the right
    x_coords_parabola = x_coords_parabola + m_coeff

    # pxy is the 2 by num array for the original x-y coordinates
    pxy = np.asarray([x_coords_parabola,y_coords_parabola]).transpose()

    # Calculate the cumulative distance
    diffs = pxy[1:, :] - pxy[:-1, :]
    dist = np.linalg.norm(diffs, axis=1)
    u = np.cumsum(dist)
    
    # Calculate the segmentation number based on Geometric_length
    new_num = round(u.max()/Geometric_length)
    Referenece_length = Geometric_length/Stretch_ratio

    # Interpolate the new x-y coordinates based on equal segmentation length
    u = np.hstack([[0], u])
    t = np.linspace(0, u[-1], new_num)
    resampled = scipy.interpolate.interpn((u,), pxy, t)

    # Reshape the data [x_coords,y_coords,z_coords] of all points as a list 
    coords = np.zeros((len(resampled),3))
    coords[:,0:2] = resampled
    print(np.shape(coords))
    print(coords.size)
    coords = coords[1:,:]
    coords_inline = coords.reshape(coords.size)
    print(np.shape(coords_inline))
    coords_inline_list = [float(i) for i in coords_inline]

    # Pre-creating ConnectorSection with axon tract Stiffness and Referenece_length
    mdb.models[ModelName].ConnectorSection(name='ConnSect-'+str(curve_num), translationalType=AXIAL)
    axon_behavior = connectorBehavior.ConnectorElasticity(components=(1, ), table=((Axon_tract_stiffness, ), ))
    mdb.models[ModelName].sections['ConnSect-'+str(curve_num)].setValues(behaviorOptions =(axon_behavior, ), u1ReferenceLength=Referenece_length )

    # Select the faces and the region
    a = mdb.models[ModelName].rootAssembly
    f1 = a.instances[InstanceName].faces
    faces1 = f1.findAt(((Length/8, -Length/6, 0.0), ))
    region2=a.Set(faces=faces1, name='Set-subcortex-for-coupling')

    # Read and preallocate coords data 
    linex = coords_inline_list
    s_x = []

    # The coordinates for each point are organized as (x_coords, y_coords, z_coords)
    for i in range(len(linex)//3):
        s_x.append(linex[0+3*i:3*i+3])
    s_x = np.array(s_x)
    s = tuple(map(tuple,s_x))

    # Preallocate tuple sss for creating wiring set
    s1 = np.zeros((len(s_x[:,0])-1,len(s_x[0,:])))
    for i in range(len(s1[:,0])):
        s1[i,:] = (s_x[i+1,:]-s_x[i,:])*3/4 + s_x[i,:]
    sss = list(map(tuple,s1))

    # Create attachment points Set
    a.AttachmentPoints(name='Attachment Points-1-Set-'+str(curve_num), points=(s),setName='Attachment Points-1-Set-'+str(curve_num))

    print(len(s_x)) 
    v = []
    for k in range(len(s_x[:,0])-1):
        v.append(a.vertices.findAt(coordinates=(s_x[k,:])))
        k = k + 1
    p = []
    for k in range(len(v)-1):
        p.append((v[k],v[k+1]))
        k = k + 1
    print(p)
    a.WirePolyLine(points=(p), mergeType=IMPRINT, meshable=OFF) 

    edges1 = a.edges.findAt(coordinates = sss)
    a.Set(edges=edges1, name='Wire-1-Set-'+str(curve_num))

    region=a.sets['Wire-1-Set-'+str(curve_num)]
    csa = a.SectionAssignment(sectionName='ConnSect-'+str(curve_num), region=region)
    a.ConnectorOrientation(orient2sameAs1=OFF, region=csa.getSet())

    region1=a.sets['Attachment Points-1-Set-'+str(curve_num)]

    mdb.models[ModelName].Coupling(name='Constraint-'+str(curve_num), controlPoint=region1, 
        surface=region2, influenceRadius=InfluenceRadius, couplingType=DISTRIBUTING, 
        weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, ur3=ON)


def Create_mesh_node_sets(ModelName,InstanceName):
    
    """ Create node sets

    Parameters
    ----------
    ModelName : string
        name of the model to be modified
    InstanceName : string
        name of the instance to be modified

    """ 

    a = mdb.models[ModelName].rootAssembly
    a.regenerate()

    # Create the node set for applying thickness perturabtion
    n1 = a.instances[InstanceName].nodes
    nodes1 = n1[84:87]
    a.Set(nodes=nodes1, name='Perturbation_nodes_set')

    # Create the node set for calculating the psi value
    n1 = a.instances[InstanceName].nodes
    nodes1 = n1[4:6]+n1[220:379]
    a.Set(nodes=nodes1, name='Topsurf_nodes')


def Compute_parabola_coeffs(num,tangent_min,tangent_max,span_min,span_max): 

    """ Compute parabola coefficients given the tangent and the span length

    Parameters
    ----------
    num : int
        number of parameter values to generate within the given ranges
    tangent_min : float
        minimal value of the range of parabola tangents at the cortical interface
    tangent_max : float
        maximum value of the range of parabola tangents at the cortical interface
    span_min : float
        minimal value of the range of parabola spans at the cortical interface
    span_max : float
        maximum value of the range of parabola spans at the cortical interface
    
    Returns
    -------
    [a_array,b_array] : list of arrays
        a_array((num,num)): coefficients array for the second-order parameter a 
        b_array((num,num)): coefficients array for the first-order parameter b
    
    """

    # Preallocate the coefficient array
    a_array = np.zeros((num,num))
    b_array = np.zeros((num,num))

    # Create the range of tangent and span
    tangent = np.linspace(tangent_min, tangent_max, num)
    span_half = np.linspace(span_min, span_max, num)/2
    x, y = sp.symbols('x y')

    # Calculate the coefficient: a and b
    for i in range(num):
        for j in range(num):
            eq1 = 2 * sp.sqrt(x) * sp.sqrt(-2 - y) - tangent[i]
            eq2 = sp.sqrt((-2 - y) / x) - span_half[j]
            sol = sp.solve((eq1, eq2), (x, y))

            a_array[j,i] =  np.array(sol)[0,0]
            b_array[j,i] =  np.array(sol)[0,1]

    return [a_array,b_array]


def Create_thickness_perturbation(ModelName,Length,Waves,Perturbation_percentage,Cortex_thickness):

    """ Create thickness perturbation

    Parameters
    ---------- 
    ModelName : string
        name of the model to be modified
    Length : float
        length of the domain
    Waves : float
        wavenumber of the perturbation, chosen so that the perturbed area has one wave
    Perturbation_percentage : float
        magnitude of the perturbation as a percentage of the cortex thickness
    Cortex_thickness : float
        cortex thickness

    """

    # Specify the node set name
    a = mdb.models[ModelName].rootAssembly
    nodes = a.sets['Perturbation_nodes_set'].nodes
    
    # Preallocate the coords vector
    coors = []
    wavelength = Length/Waves
    perturbation = Perturbation_percentage * Cortex_thickness

    for n in nodes:
        x = n.coordinates[0]
        y = n.coordinates[1]
        z = n.coordinates[2]

        p_n = cos(x/wavelength*2.*math.pi)

        y_p = y - perturbation * p_n
        
        coors.append((x, y_p, z))

    a.editNode(nodes=nodes, coordinates=coors) 


def GetKeywordPosition(blockPrefix,ModelName):
    
    """ Get the keyword position from the inp file to request variables history output

    Parameters
    ---------- 
    blockPrefix : string
        string to search for; text will be added after this line
    ModelName : string
        name of the model to be modified

    """

    m = mdb.models[ModelName]
    m.keywordBlock.synchVersions() 
    pos = 0
    positions = []
    for block in m.keywordBlock.sieBlocks:
        if string.lower(block[0:len(blockPrefix)])==string.lower(blockPrefix):
            positions.append(pos)
        pos=pos+1
    return positions


def Modify_input(ModelName):

    """ Modify the command line to request top surface nodal coordinates for calculating psi values

    Parameters
    ---------- 
    ModelName : string
        name of the model to be modified

    """

    blocks = GetKeywordPosition('*Output, field, variable=PRESELECT',ModelName)
    for b in reversed(blocks):
        mdb.models[ModelName].keywordBlock.insert(position=b, text= '*Output, field, frequency=999\n*Node Output, nset=Topsurf_nodes\nCOORD,  ')


def Modify_input_for_wiring_nodeSet(ModelName):

    """ Modify the command line to request connector nodal coordinates for calculating wiring length

    Parameters
    ---------- 
    ModelName : string
        name of the model to be modified

    """

    blocks = GetKeywordPosition('*Output, field, variable=PRESELECT',ModelName)
    for b in reversed(blocks):
        mdb.models[ModelName].keywordBlock.insert(position=b, text= '*Output, field\n*Node Output, nset="Attachment Points-1-Set-1"\nCOORD,\n*Node Output,nset="Attachment Points-1-Set-2"\nCOORD, ')


def Modify_input_for_initialize_growth_variable(ModelName):
   
    """ Modify the command line to initialize growth variable

    Parameters
    ---------- 
    ModelName : string
        name of the model to be modified

    """

    blocks = GetKeywordPosition('*User Material, constants=5\n100., 930., 300.,2790., 0.05',ModelName)
    for b in reversed(blocks):
        mdb.models[ModelName].keywordBlock.insert(position=b, text= '*Initial Conditions, Type=Solution, User')

def Create_Job(ModelName, JobName, USERSUB):

    """ Create Job and execute the simulation

    Parameters
    ----------
    ModelName : string
        name of the model to be modified
    JobName : string
        name of the job file to be created
    USERSUB : string
        filename of the user subroutine to be used for the job

    """

    mdb.Job(name=JobName, model=ModelName, description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, 
        userSubroutine=USERSUB, scratch='', 
        resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=8, 
        numDomains=8, numGPUs=0)
    mdb.jobs[JobName].writeInput(consistencyChecking=OFF)
    try: 
        mdb.jobs[JobName].submit(consistencyChecking=OFF)
        mdb.jobs[JobName].waitForCompletion()
    except: 
        print "Job %s crashed" % JobName


def Post_processing_odbs(ODB_Name,NodesetName,Step,Variable):

    """ Post-processing odb files: extract x and y coordinates for the given node set

    Parameters
    ----------
    ODB_Name : string
        name of .odb file to be read
    NodesetName : string
        name of node set to be read
    Step : string
        name of the step to read
    Variable : string
        name of variable to be read from the output file: 'COORD'

    Returns
    -------
    [x_coords_odb, y_coords_odb] : arrays
        nx2 array containing the x and y coordinates of the specified node set
        
    """

    # Specify the odb file namne    
    odb = openOdb(ODB_Name)
    mySet = odb.rootAssembly.nodeSets[NodesetName]

    # Specify the last increment
    lastFrame = odb.steps[Step].frames[-1]

    # Specify the variables 
    coords=lastFrame.fieldOutputs[Variable]
    mySetCoord = coords.getSubset(region=mySet)
    mySetCoordValues=mySetCoord.values
    
    # Preallocate the coords vector
    coords_empty = []

    for curNode in mySetCoordValues:
        coords_empty.append(curNode.data)

    # Reshape the coords array      
    coords_empty = np.array(coords_empty)
    num_rows = np.shape(coords_empty)[0]
    num_cols = np.shape(coords_empty)[1]
    new_array = coords_empty.reshape(num_rows,num_cols)

    # Adjust the position of the second row and the last row (Abaqus has different node numbering for the nodes at the boundary)
    index = 1
    v = new_array[index,:].reshape(1,2)
    new_array = np.delete(new_array, (index), axis=0)
    new_array = np.concatenate((new_array, v), axis=0)

    # Obtain the x and y coordinates
    x_coords_odb = new_array[:,0]
    y_coords_odb = new_array[:,1]

    return [x_coords_odb, y_coords_odb]


def Calculate_psi(ODBtype,NodesetName,Step,Variable,parameter1_range,parameter2_range,x_coords_bench,y_coords_bench):

    """ Calculate mean squared displacements psi for different odb files
 
    Parameters
    ----------
    ODBtype : integer (0, 1, 2)
        which type of ODB file is to be read
    NodesetName : string
        name of node set to be read
    Step : string
        name of the step to be read
        Variable : string
        name of variable to be read from the output file: 'COORD'
    parameter1_range : array
        range of the first parameter for the resulting heatmap
    parameter2_range : array
        range of the second parameter for the resulting heatmap
    x_coords_bench : list of floats
        deformed x coordinates from the benchmark case
    y_coords_bench : list of floats
        deformed y coordinates from the benchmark case

    Returns
    -------
    psi : number
        psi value 

    """

    # Specify the odb file type
    if ODBtype == 0:
        ODB_Name_template = 'Job-span%d-tangent%d.odb'
    elif ODBtype == 1:
        ODB_Name_template = 'Job-perturbation%d-axonstiffness%d.odb'
    elif ODBtype == 2:
        ODB_Name_template = 'Job-primary%d-secondary%d.odb'
    
    # Preallocate the psi value array
    psi = np.zeros((len(parameter1_range),len(parameter2_range)))

    # Loop through different odb files
    for i in range(0,len(parameter1_range)):
        for j in range(0,len(parameter2_range)):

            ODB_Name = ODB_Name_template %(i,j)
            [x_coords_odb,y_coords_odb] = Post_processing_odbs(ODB_Name,NodesetName,Step,Variable)

            # Compute the mean squared displacements based on benchmark coordinates
            psi[i,j]  = np.sum(np.sqrt((x_coords_odb-x_coords_bench)**2 + (y_coords_odb-y_coords_bench)**2)) / len(x_coords_odb)

    return psi


def Calculate_wiring_length(ODB_Name,NodeSetName,FrameNumber):

    """ Calculate wiring length of the axon tract at the specific increment
 
    Parameters
    ----------
    ODB_Name : string
        name of .odb file to be read
    NodesetName : string
        name of node set to be read
    FrameNumber : number
        increment (time step) to be read

    Returns
    -------
    length : number
        wiring length based on coordinates of every attachment point

    """

    # Specify the odb file namne
    odb = openOdb(ODB_Name)

    # Specify the Node set name and the frame number
    mySet = odb.rootAssembly.nodeSets[NodeSetName]
    Frame = odb.steps['Step-1'].frames[FrameNumber]

    # Specify the variables 
    coords = Frame.fieldOutputs['COORD']
    mySetCoord = coords.getSubset(region=mySet)
    mySetCoordValues = mySetCoord.values

    # Preallocate the coords vector
    coords_empty = []

    for curNode in mySetCoordValues:
        coords_empty.append(curNode.data)
    
    # Reshape the coords array  
    coords_empty = np.array(coords_empty)
    num_rows = np.shape(coords_empty)[0]
    num_cols = np.shape(coords_empty)[1]
    coords_array = coords_empty.reshape(num_rows,num_cols)
    coords_array = coords_array[1:,:]

    # Obtain the x and y coordinates
    x_coords = coords_array[:,0]
    y_coords = coords_array[:,1]

    # Calculate the difference
    x_coords_diff = np.diff(x_coords)
    y_coords_diff = np.diff(y_coords)

    # Calculate the total length
    length = np.sum(np.sqrt(np.square(x_coords_diff)+np.square(y_coords_diff)))

    return length