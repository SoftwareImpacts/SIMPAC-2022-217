"""
The output data from the first script is used in this code which is run
in Abaqus to create the sarcomeres.
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
#--- HIERARCHICAL LEVEL 1 - SARCOMERE -----------------------------------------------------------------------
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
    cross_section = str(read_input_entry[0]).rstrip('\n')
    r_myosin = float(read_input_entry[1])
    l_myosin = float(read_input_entry[2])
    r_actin = float(read_input_entry[3])
    l_actin = float(read_input_entry[4])
    thickness_z_disc = float(read_input_entry[5])
    l_sarcomere = float(read_input_entry[6])
    l_crossbridge = float(read_input_entry[7])
    number_crossbridges = int(read_input_entry[8])
    input_entry.close()


    #------ GEOMETRY NAME -----------------------------------------------------------------------------------

    # prepare a newly named model and delete the default model
    model = mdb.Model(name='sarcomere')
    assembly = model.rootAssembly
    if mdb.Model(name='Model-1'):
        del mdb.models['Model-1']


    #--------- Z-DISC ---------------------------------------------------------------------------------------

    """
    Two thin z-discs are created which define both ends of the sarcomere. If only level
    1 is modeled, a cylindrical z-disc should be created. For multiscale modeling purposes,
    a hexagonal z-disc cross-section is modeled with only half of the z-disc's thickness.
    """

    r_unit_cell = r_myosin+l_crossbridge+2.0*r_actin
    if cross_section=='cir':
        sketch_z_disc = model.ConstrainedSketch(name='sketch_z-disc', sheetSize=200.0)
        sketch_z_disc.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(r_unit_cell, 0.0))
        z_disc = model.Part(dimensionality=THREE_D, name='z-disc', type=DEFORMABLE_BODY)
        z_disc.BaseSolidExtrude(depth=thickness_z_disc, sketch=sketch_z_disc)
        assembly.Instance(name='z-disc1', part=z_disc, dependent=ON)
        assembly.Instance(name='z-disc2', part=z_disc, dependent=ON)
        assembly.translate(instanceList=('z-disc2', ), vector=(0.0, 0.0, l_sarcomere-thickness_z_disc))
    elif cross_section=='hex':
        sketch_z_disc = model.ConstrainedSketch(name='sketch_z-disc', sheetSize=200.0)
        sketch_z_disc.Line(point1=(0, r_unit_cell-r_actin), point2=(sin(pi/3)*(r_unit_cell-r_actin), cos(pi/3)*(r_unit_cell-r_actin)))
        sketch_z_disc.Line(point1=(sin(pi/3)*(r_unit_cell-r_actin), cos(pi/3)*(r_unit_cell-r_actin)), point2=(sin(pi/3)*(r_unit_cell-r_actin), -cos(pi/3)*(r_unit_cell-r_actin)))
        sketch_z_disc.Line(point1=(sin(pi/3)*(r_unit_cell-r_actin), -cos(pi/3)*(r_unit_cell-r_actin)), point2=(0, -(r_unit_cell-r_actin)))
        sketch_z_disc.Line(point1=(0, -(r_unit_cell-r_actin)), point2=(-sin(pi/3)*(r_unit_cell-r_actin), -cos(pi/3)*(r_unit_cell-r_actin)))
        sketch_z_disc.Line(point1=(-sin(pi/3)*(r_unit_cell-r_actin), -cos(pi/3)*(r_unit_cell-r_actin)), point2=(-sin(pi/3)*(r_unit_cell-r_actin), cos(pi/3)*(r_unit_cell-r_actin)))
        sketch_z_disc.Line(point1=(-sin(pi/3)*(r_unit_cell-r_actin), cos(pi/3)*(r_unit_cell-r_actin)), point2=(0, r_unit_cell-r_actin))
        z_disc = model.Part(dimensionality=THREE_D, name='z-disc', type=DEFORMABLE_BODY)
        z_disc.BaseSolidExtrude(depth=thickness_z_disc/2.0, sketch=sketch_z_disc)
        assembly.Instance(name='z-disc1', part=z_disc, dependent=ON)
        assembly.Instance(name='z-disc2', part=z_disc, dependent=ON)
        assembly.translate(instanceList=('z-disc2', ), vector=(0.0, 0.0, l_sarcomere-thickness_z_disc/2.0))


    #--------- MYOSIN ---------------------------------------------------------------------------------------

    """
    One cylindrical myosin filament is generated.
    """

    sketch_myosin = model.ConstrainedSketch(name='sketch_myosin', sheetSize=20.0)
    sketch_myosin.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(r_myosin, 0.0))
    myosin = model.Part(name='myosin', dimensionality=THREE_D,type=DEFORMABLE_BODY)
    myosin.BaseSolidExtrude(sketch=sketch_myosin, depth=l_myosin)

    assembly.Instance(name='myosin', part=myosin, dependent=ON)
    assembly.translate(instanceList=('myosin', ), vector=(0.0, 0.0, (l_sarcomere-l_myosin)/2.0))


    #--------- ACTIN ----------------------------------------------------------------------------------------

    """
    Twelve actin filaments are generated and moved to the right place. Whole
    cylindrical actin filaments are only created for a unit cell with
    cylindrical cross-section. For multiscale modeling purposes, hexagonal
    unit cell cross-sections are used where the actin filaments are cutted.
    """

    sketch_actin = model.ConstrainedSketch(name='sketch_actin', sheetSize=20.0)
    sketch_actin.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(r_actin, 0.0))
    actin = model.Part(name='actin', dimensionality=THREE_D,type=DEFORMABLE_BODY)
    actin.BaseSolidExtrude(sketch=sketch_actin, depth=l_actin)

    if cross_section=='hex':
        # create one xy-plane (d1[2]) and one axis (d1[5]) through two points (d1[3] and d1[4])
        d1 = actin.datums
        actin.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0)
        actin.DatumPointByCoordinate(coords=(0, 0, 0))
        actin.DatumPointByCoordinate(coords=(0, r_actin, 0))
        actin.DatumAxisByTwoPoint(point1=d1[3], point2=d1[4])

        # cut out actin
        t = actin.MakeSketchTransform(sketchPlane=d1[2], sketchUpEdge=d1[5],
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0))
        sketch_actin_cut = model.ConstrainedSketch(name='sketch_actin_cut',
            sheetSize=r_actin*2, gridSpacing=r_actin/50, transform=t)
        sketch_actin_cut.Line(point1=(0.0, 0.0), point2=(-sin(pi/3.0)*r_actin, -cos(pi/3.0)*r_actin))
        sketch_actin_cut.Line(point1=(0.0, 0.0), point2=(sin(pi/3.0)*r_actin, -cos(pi/3.0)*r_actin))
        sketch_actin_cut.ArcByCenterEnds(center=(0.0, 0.0), point1=(-sin(pi/3.0)*r_actin, -cos(pi/3.0)*r_actin), point2=(
            sin(pi/3.0)*r_actin, -cos(pi/3.0)*r_actin), direction=CLOCKWISE)
        actin.CutExtrude(sketchPlane=d1[2], sketchUpEdge=d1[5], sketchPlaneSide=SIDE1,
            sketchOrientation=RIGHT, sketch=sketch_actin_cut, flipExtrudeDirection=ON)

    assembly.Instance(name='actin1', part=actin, dependent=ON)
    assembly.Instance(name='actin2', part=actin, dependent=ON)
    assembly.Instance(name='actin3', part=actin, dependent=ON)
    assembly.Instance(name='actin4', part=actin, dependent=ON)
    assembly.Instance(name='actin5', part=actin, dependent=ON)
    assembly.Instance(name='actin6', part=actin, dependent=ON)
    assembly.Instance(name='actin7', part=actin, dependent=ON)
    assembly.Instance(name='actin8', part=actin, dependent=ON)
    assembly.Instance(name='actin9', part=actin, dependent=ON)
    assembly.Instance(name='actin10', part=actin, dependent=ON)
    assembly.Instance(name='actin11', part=actin, dependent=ON)
    assembly.Instance(name='actin12', part=actin, dependent=ON)

    if cross_section=='hex':
        assembly.rotate(instanceList=('actin2', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=-60.0)
        assembly.rotate(instanceList=('actin3', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=-120.0)
        assembly.rotate(instanceList=('actin4', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=-180.0)
        assembly.rotate(instanceList=('actin5', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=-240.0)
        assembly.rotate(instanceList=('actin6', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=-300.0)
        assembly.rotate(instanceList=('actin8', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=-60.0)
        assembly.rotate(instanceList=('actin9', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=-120.0)
        assembly.rotate(instanceList=('actin10', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=-180.0)
        assembly.rotate(instanceList=('actin11', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=-240.0)
        assembly.rotate(instanceList=('actin12', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=-300.0)

        assembly.translate(instanceList=('actin1', ), vector=(0.0, r_unit_cell-r_actin, thickness_z_disc/2.0))
        assembly.translate(instanceList=('actin2', ), vector=(sin(pi/3.0)*(r_unit_cell-r_actin), cos(pi/3.0)*(r_unit_cell-r_actin), thickness_z_disc/2.0))
        assembly.translate(instanceList=('actin3', ), vector=(sin(pi/3.0)*(r_unit_cell-r_actin), -cos(pi/3.0)*(r_unit_cell-r_actin), thickness_z_disc/2.0))
        assembly.translate(instanceList=('actin4', ), vector=(0.0, -(r_unit_cell-r_actin), thickness_z_disc/2.0))
        assembly.translate(instanceList=('actin5', ), vector=(-sin(pi/3.0)*(r_unit_cell-r_actin), -cos(pi/3.0)*(r_unit_cell-r_actin), thickness_z_disc/2.0))
        assembly.translate(instanceList=('actin6', ), vector=(-sin(pi/3.0)*(r_unit_cell-r_actin), cos(pi/3.0)*(r_unit_cell-r_actin), thickness_z_disc/2.0))
        assembly.translate(instanceList=('actin7', ), vector=(0.0, r_unit_cell-r_actin, l_sarcomere-thickness_z_disc/2.0-l_actin))
        assembly.translate(instanceList=('actin8', ), vector=(sin(pi/3.0)*(r_unit_cell-r_actin), cos(pi/3.0)*(r_unit_cell-r_actin), l_sarcomere-thickness_z_disc/2.0-l_actin))
        assembly.translate(instanceList=('actin9', ), vector=(sin(pi/3.0)*(r_unit_cell-r_actin), -cos(pi/3.0)*(r_unit_cell-r_actin), l_sarcomere-thickness_z_disc/2.0-l_actin))
        assembly.translate(instanceList=('actin10', ), vector=(0.0, -(r_unit_cell-r_actin), l_sarcomere-thickness_z_disc/2.0-l_actin))
        assembly.translate(instanceList=('actin11', ), vector=(-sin(pi/3.0)*(r_unit_cell-r_actin), -cos(pi/3.0)*(r_unit_cell-r_actin), l_sarcomere-thickness_z_disc/2.0-l_actin))
        assembly.translate(instanceList=('actin12', ), vector=(-sin(pi/3.0)*(r_unit_cell-r_actin), cos(pi/3.0)*(r_unit_cell-r_actin), l_sarcomere-thickness_z_disc/2.0-l_actin))

    else:
        assembly.translate(instanceList=('actin1', ), vector=(0.0, r_unit_cell-r_actin, thickness_z_disc))
        assembly.translate(instanceList=('actin2', ), vector=(sin(pi/3.0)*(r_unit_cell-r_actin), cos(pi/3.0)*(r_unit_cell-r_actin), thickness_z_disc))
        assembly.translate(instanceList=('actin3', ), vector=(sin(pi/3.0)*(r_unit_cell-r_actin), -cos(pi/3.0)*(r_unit_cell-r_actin), thickness_z_disc))
        assembly.translate(instanceList=('actin4', ), vector=(0.0, -(r_unit_cell-r_actin), thickness_z_disc))
        assembly.translate(instanceList=('actin5', ), vector=(-sin(pi/3.0)*(r_unit_cell-r_actin), -cos(pi/3.0)*(r_unit_cell-r_actin), thickness_z_disc))
        assembly.translate(instanceList=('actin6', ), vector=(-sin(pi/3.0)*(r_unit_cell-r_actin), cos(pi/3.0)*(r_unit_cell-r_actin), thickness_z_disc))
        assembly.translate(instanceList=('actin7', ), vector=(0.0, r_unit_cell-r_actin, l_sarcomere-thickness_z_disc-l_actin))
        assembly.translate(instanceList=('actin8', ), vector=(sin(pi/3.0)*(r_unit_cell-r_actin), cos(pi/3.0)*(r_unit_cell-r_actin), l_sarcomere-thickness_z_disc-l_actin))
        assembly.translate(instanceList=('actin9', ), vector=(sin(pi/3.0)*(r_unit_cell-r_actin), -cos(pi/3.0)*(r_unit_cell-r_actin), l_sarcomere-thickness_z_disc-l_actin))
        assembly.translate(instanceList=('actin10', ), vector=(0.0, -(r_unit_cell-r_actin), l_sarcomere-thickness_z_disc-l_actin))
        assembly.translate(instanceList=('actin11', ), vector=(-sin(pi/3.0)*(r_unit_cell-r_actin), -cos(pi/3.0)*(r_unit_cell-r_actin), l_sarcomere-thickness_z_disc-l_actin))
        assembly.translate(instanceList=('actin12', ), vector=(-sin(pi/3.0)*(r_unit_cell-r_actin), cos(pi/3.0)*(r_unit_cell-r_actin), l_sarcomere-thickness_z_disc-l_actin))



    #--------- CROSSBRIDGES ---------------------------------------------------------------------------------

    """
    The number of crossbridges are equally distributed between the left and right
    actin filaments in the sarcomere. At each crossbridge connection along the
    myosin filament, two opposing crossbridges are linked to two actin filaments.
    Here, the two crossbridges at neighboring connections along myosin are rotated
    by an angle of 60 degrees relative to the previous crossbridge pair. Depending
    on the required number of crossbridges, fewer crossbridge connections may be
    generated.
    """

    sketch_crossbridge = model.ConstrainedSketch(name='sketch_crossbridge', sheetSize=200.0)
    sketch_crossbridge.Line(point1=(0.0, r_myosin), point2=(0.0, r_myosin+l_crossbridge))
    crossbridge = model.Part(name='crossbridge', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    crossbridge.BaseWire(sketch=sketch_crossbridge)

    rows_crossbridges=int(number_crossbridges//4.0*2.0)
    if cross_section=='hex':
        length_for_crossbridges=l_myosin-(l_sarcomere-thickness_z_disc-2.0*l_actin)
    else:
        length_for_crossbridges=l_myosin-(l_sarcomere-2.0*thickness_z_disc-2.0*l_actin)
    if rows_crossbridges>2:
        distance_between_crossbridges=length_for_crossbridges/2.0/(rows_crossbridges/2.0-1.0)
    else:
        distance_between_crossbridges=0

    angle_crossbridge_right=0.0
    angle_crossbridge_left=0.0

    # iterate over all connection sites along the myosin filament
    for i_row in range(rows_crossbridges):

        # connection sites on the right side of the myosin filament
        if i_row<rows_crossbridges/2.0:
            assembly.Instance(name='crossbridge'+str(i_row*2+1), part=crossbridge, dependent=ON)
            assembly.Instance(name='crossbridge'+str(i_row*2+2), part=crossbridge, dependent=ON)

            assembly.translate(instanceList=('crossbridge'+str(i_row*2+1), ), vector=(0.0, 0.0, (l_sarcomere-l_myosin)/2.0+i_row*distance_between_crossbridges))
            assembly.translate(instanceList=('crossbridge'+str(i_row*2+2), ), vector=(0.0, 0.0, (l_sarcomere-l_myosin)/2.0+i_row*distance_between_crossbridges))
            if angle_crossbridge_right == 180.0:
                angle_crossbridge_right = 0.0
            assembly.rotate(instanceList=('crossbridge'+str(i_row*2+1), ), axisPoint=(0.0, 0.0, 0), axisDirection=(0.0, 0.0, 1), angle=0.0+angle_crossbridge_right)
            assembly.rotate(instanceList=('crossbridge'+str(i_row*2+2), ), axisPoint=(0.0, 0.0, 0), axisDirection=(0.0, 0.0, 1), angle=180.0+angle_crossbridge_right)
            angle_crossbridge_right=angle_crossbridge_right+60.0

        # connection sites on the left side of the myosin filament
        else:
            assembly.Instance(name='crossbridge'+str(i_row*2+1), part=crossbridge, dependent=ON)
            assembly.Instance(name='crossbridge'+str(i_row*2+2), part=crossbridge, dependent=ON)

            assembly.translate(instanceList=('crossbridge'+str(i_row*2+1), ), vector=(0.0, 0.0, (l_sarcomere+l_myosin)/2.0-(i_row-rows_crossbridges/2.0)*distance_between_crossbridges))
            assembly.translate(instanceList=('crossbridge'+str(i_row*2+2), ), vector=(0.0, 0.0, (l_sarcomere+l_myosin)/2.0-(i_row-rows_crossbridges/2.0)*distance_between_crossbridges))
            if angle_crossbridge_left == 180.0:
                angle_crossbridge_left = 0.0
            assembly.rotate(instanceList=('crossbridge'+str(i_row*2+1), ), axisPoint=(0.0, 0.0, 0), axisDirection=(0.0, 0.0, 1), angle=0.0+angle_crossbridge_left)
            assembly.rotate(instanceList=('crossbridge'+str(i_row*2+2), ), axisPoint=(0.0, 0.0, 0), axisDirection=(0.0, 0.0, 1), angle=180.0+angle_crossbridge_left)
            angle_crossbridge_left=angle_crossbridge_left+60.0

    # print the number of the generated crossbridges
    print rows_crossbridges*2,"crossbridges are generated."


    #--------- TITIN ----------------------------------------------------------------------------------------

    '''
    Six titin filaments are created on each side of the sarcomere which connect
    myosin to the z-disc and to each actin filament.
    '''

    if cross_section=='hex':
        l_titin = sqrt((l_sarcomere-thickness_z_disc-l_myosin)**2/4.0+(r_myosin+l_crossbridge)**2)
        angle_titin_rad = atan((r_myosin+l_crossbridge)/((l_sarcomere-thickness_z_disc-l_myosin)/2.0))
        angle_titin_deg = atan((r_myosin+l_crossbridge)/((l_sarcomere-thickness_z_disc-l_myosin)/2.0))*(180/pi)
    else:
        l_titin = sqrt((l_sarcomere-thickness_z_disc*2.0-l_myosin)**2/4.0+(r_myosin+l_crossbridge)**2)
        angle_titin_rad = atan((r_myosin+l_crossbridge)/((l_sarcomere-thickness_z_disc*2.0-l_myosin)/2.0))
        angle_titin_deg = atan((r_myosin+l_crossbridge)/((l_sarcomere-thickness_z_disc*2.0-l_myosin)/2.0))*(180/pi)

    sketch_titin = model.ConstrainedSketch(name='sketch_titin', sheetSize=200.0)
    sketch_titin.Line(point1=(0.0, 0.0), point2=(l_titin, 0.0))
    titin = model.Part(name='titin', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    titin.BaseWire(sketch=sketch_titin)

    # iterate over all titin filaments
    for i_titin in range(12):

        # titin filaments on the left side of the sarcomere
        if i_titin<6:
            assembly.Instance(name='titin'+str(i_titin+1), part=titin, dependent=ON)
            assembly.rotate(instanceList=('titin'+str(i_titin+1), ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 1.0, 0.0), angle=90.0)
            assembly.rotate(instanceList=('titin'+str(i_titin+1), ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(-cos(pi/3.0*i_titin)*r_myosin, sin(pi/3.0*i_titin)*r_myosin, 0.0), angle=angle_titin_deg)
            assembly.translate(instanceList=('titin'+str(i_titin+1), ), vector=(0.0, 0.0, (l_sarcomere-l_myosin)/2.0))

        # titin filaments on the right side of the sarcomere
        else:
            assembly.Instance(name='titin'+str(i_titin+1), part=titin, dependent=ON)
            assembly.rotate(instanceList=('titin'+str(i_titin+1), ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 1.0, 0.0), angle=90.0)
            assembly.rotate(instanceList=('titin'+str(i_titin+1), ), axisPoint=(0.0, 0.0, -l_titin), axisDirection=(-cos(pi/3.0*i_titin)*r_myosin, sin(pi/3.0*i_titin)*r_myosin, 0.0), angle=angle_titin_deg)
            if cross_section=='hex':
                assembly.translate(instanceList=('titin'+str(i_titin+1), ), vector=(0.0, 0.0, l_sarcomere-thickness_z_disc/2.0+l_titin*(1.0-cos(angle_titin_rad))))
            else:
                assembly.translate(instanceList=('titin'+str(i_titin+1), ), vector=(0.0, 0.0, l_sarcomere-thickness_z_disc+l_titin*(1.0-cos(angle_titin_rad))))


    #--------- MERGE TO ONE INSTANCE ------------------------------------------------------------------------

    # merge the different instances of the sarcomere to one instance
    all_instances = assembly.instances.keys()
    assembly.InstanceFromBooleanMerge(name='sarcomere', instances=([assembly.instances[all_instances[i]]
        for i in range(len(assembly.instances))] ), keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
    assembly.features.changeKey(fromName='sarcomere-1', toName='sarcomere')


    #--------- SAVE MODEL AS CAE-FILE -----------------------------------------------------------------------

    file_ending=input_entries[i_structure].replace('input_values_sarcomere','')
    file_ending=file_ending.replace(str(number_crossbridges),str(rows_crossbridges*2))
    file_ending=file_ending.replace('.txt','')

    file='sarcomere'+str(file_ending)
    mdb.saveAs(str(file)+'.cae')
