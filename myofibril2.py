"""
The output data from the first script is used in this code which is run
in Abaqus to create the myofibrils.
Both scripts need to be located in the same folder. Use this folder as
the work directory in Abaqus.
"""

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


#------------------------------------------------------------------------------------------------------------
#--- HIERARCHICAL LEVEL 2 - MYOFIBRIL -----------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


#------ OUTPUT OF FIRST SCRIPT ------------------------------------------------------------------------------

"""
From the first script, read the output data of the text file where the names of
all text files for the input variables are saved.
"""

all_txt_input="all_txt_input.txt"
txt_input = open(all_txt_input,'r')
input_entries = txt_input.read().splitlines()
txt_input.close


#------ GENERATION OF THE STRUCTURES ------------------------------------------------------------------------

# loop to create all structures
for i_structure in range(len(input_entries)):

    # read the input values from each structure of the first script
    input_entry = open(input_entries[i_structure],'r')
    read_input_entry = input_entry.readlines()
    l_sarcomere = float(read_input_entry[0])
    number_sarcomeres = int(read_input_entry[1])
    r_unit_cell = float(read_input_entry[2])
    number_unit_cells_per_sarcomere = int(read_input_entry[3])
    input_entry.close()


    #------ INITIAL VALUES ----------------------------------------------------------------------------------

    arrangement_angle=0
    unit_cell_circle=0
    calc_variable=0
    number_unit_cells_in_circle=6
    distance_x_unit_cell=0
    distance_y_unit_cell=0


    #------ GEOMETRY NAME -----------------------------------------------------------------------------------

    # prepare a newly named model and delete the default model
    model = mdb.Model(name='myofibril')
    assembly = model.rootAssembly
    if mdb.Model(name='Model-1'):
        del mdb.models['Model-1']


    #--------- UNIT CELLS -----------------------------------------------------------------------------------

    '''
    The hexagonal unit cells per sarcomere are generated and placed
    circumferentially in the sarcomere. Inner unit cells are completely
    surrounded circumferentially by further unit cells. If the number
    of the remaining unit cells for another circle are not sufficient
    to completely surround the inner unit cells, the loop will be
    terminated after the last complete unit cell circle.
    '''

    sketch_unit_cell = model.ConstrainedSketch(name='sketch_unit_cell',sheetSize=20.0)
    sketch_unit_cell.Line(point1=(0, r_unit_cell), point2=(sin(pi/3)*r_unit_cell, cos(pi/3)*r_unit_cell))
    sketch_unit_cell.Line(point1=(sin(pi/3)*r_unit_cell, cos(pi/3)*r_unit_cell), point2=(sin(pi/3)*r_unit_cell, -cos(pi/3)*r_unit_cell))
    sketch_unit_cell.Line(point1=(sin(pi/3)*r_unit_cell, -cos(pi/3)*r_unit_cell), point2=(0, -r_unit_cell))
    sketch_unit_cell.Line(point1=(0, -r_unit_cell), point2=(-sin(pi/3)*r_unit_cell, -cos(pi/3)*r_unit_cell))
    sketch_unit_cell.Line(point1=(-sin(pi/3)*r_unit_cell, -cos(pi/3)*r_unit_cell), point2=(-sin(pi/3)*r_unit_cell, cos(pi/3)*r_unit_cell))
    sketch_unit_cell.Line(point1=(-sin(pi/3)*r_unit_cell, cos(pi/3)*r_unit_cell), point2=(0, r_unit_cell))

    unit_cell = model.Part(name='unit_cell', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    unit_cell.BaseSolidExtrude(sketch=sketch_unit_cell, depth=l_sarcomere)

    # iterate over all unit cells per sarcomere
    for i_cell in range(1,number_unit_cells_per_sarcomere+1):

        # place the unit cells in the sarcomere
        if (round(pi/6,3)==round(arrangement_angle,3)) or (round(pi/6+2*pi/6,3)==round(arrangement_angle,3)) or (round(pi/6+2*pi/6*2,3)==round(arrangement_angle,3)) \
        or (round(pi/6+2*pi/6*3,3)==round(arrangement_angle,3)) or (round(pi/6+2*pi/6*4,3)==round(arrangement_angle,3)) or (round(pi/6+2*pi/6*5,3)==round(arrangement_angle,3)):
            distance_x_unit_cell=sin(arrangement_angle)*(cos(pi/6)*r_unit_cell)*2*unit_cell_circle
            distance_y_unit_cell=cos(arrangement_angle)*(cos(pi/6)*r_unit_cell)*2*unit_cell_circle
        else:
            if pi/6<=arrangement_angle<=pi/2:
                distance_x_unit_cell=distance_x_unit_cell+cos(pi/6)*r_unit_cell
                distance_y_unit_cell=distance_y_unit_cell-(r_unit_cell+sin(pi/6)*r_unit_cell)
            if pi/2<=arrangement_angle<=5./6*pi:
                distance_x_unit_cell=distance_x_unit_cell-cos(pi/6)*r_unit_cell
                distance_y_unit_cell=distance_y_unit_cell-(r_unit_cell+sin(pi/6)*r_unit_cell)
            if 5./6*pi<=arrangement_angle<=7./6*pi:
                distance_x_unit_cell=distance_x_unit_cell-cos(pi/6)*r_unit_cell*2
                distance_y_unit_cell=distance_y_unit_cell
            if 7./6*pi<=arrangement_angle<=3./2*pi:
                distance_x_unit_cell=distance_x_unit_cell-cos(pi/6)*r_unit_cell
                distance_y_unit_cell=distance_y_unit_cell+(r_unit_cell+sin(pi/6)*r_unit_cell)
            if 3./2*pi<=arrangement_angle<=11./6*pi:
                distance_x_unit_cell=distance_x_unit_cell+cos(pi/6)*r_unit_cell
                distance_y_unit_cell=distance_y_unit_cell+(r_unit_cell+sin(pi/6)*r_unit_cell)
            if 11./6*pi<=arrangement_angle<=13./6*pi:
                distance_x_unit_cell=distance_x_unit_cell+cos(pi/6)*r_unit_cell*2
                distance_y_unit_cell=distance_y_unit_cell

        assembly.Instance(name='unit_cell'+str(i_cell), part=unit_cell, dependent=ON)
        assembly.translate(instanceList=('unit_cell'+str(i_cell), ), vector=(distance_x_unit_cell, distance_y_unit_cell, 0))

        arrangement_angle=arrangement_angle+2*pi/number_unit_cells_in_circle

        # if inner unit cells are completely surrounded, the number of
        # unit cells for the next circle are calculated
        if i_cell==1+calc_variable*6:
            calc_variable=calc_variable+(1+unit_cell_circle)
            number_unit_cells_in_circle=(unit_cell_circle+1)*6
            unit_cell_circle=unit_cell_circle+1
            arrangement_angle=pi/6
            # terminate the loop if the number of the remaining unit cells
            # are not sufficient for a complete further circle
            if number_unit_cells_per_sarcomere-i_cell<number_unit_cells_in_circle:
                break

    # print the number of the generated unit cells per sarcomere
    print i_cell,"unit cell/s per sarcomere is/are generated."


    #--------- SARCOMERES -----------------------------------------------------------------------------------

    """
    The required number of sarcomeres are created. Each sarcomere consists of
    the set of previously generated unit cells.
    """

    instances_unit_cells_per_sarcomere = assembly.instances.keys()
    assembly.InstanceFromBooleanMerge(name='sarcomere', instances=([assembly.instances[instances_unit_cells_per_sarcomere[i]]
        for i in range(i_cell)] ), keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)
    assembly.features.changeKey(fromName='sarcomere-1', toName='sarcomere1')
    sarcomere=model.parts['sarcomere']

    if number_sarcomeres>1:
        # iterate over the number of sarcomeres
        for i_sarcomere in range(2,number_sarcomeres+1):
            assembly.Instance(name='sarcomere'+str(i_sarcomere), part=sarcomere, dependent=ON)
            assembly.translate(instanceList=('sarcomere'+str(i_sarcomere), ), vector=(0.0, 0.0, (i_sarcomere-1)*l_sarcomere))


    #--------- MERGE TO ONE INSTANCE ------------------------------------------------------------------------

    # merge the different instances of the myofibril to one instance
    all_instances = assembly.instances.keys()
    assembly.InstanceFromBooleanMerge(name='myofibril', instances=([assembly.instances[all_instances[i]]
        for i in range(len(assembly.instances))] ), keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
    assembly.features.changeKey(fromName='myofibril-1', toName='myofibril')


    #--------- SAVE MODEL AS CAE-FILE -----------------------------------------------------------------------

    file_ending=input_entries[i_structure].replace('input_values_myofibril','')
    file_ending=file_ending.replace(str(number_unit_cells_per_sarcomere),str(i_cell))
    file_ending=file_ending.replace('.txt','')

    file='myofibril'+str(file_ending)
    mdb.saveAs(str(file)+'.cae')
