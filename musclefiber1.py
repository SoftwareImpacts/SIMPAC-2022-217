"""
In the first of two scripts, voronoi diagrams generate polygonal
myofibrils which are used in the second script to create 3D muscle
fiber structures in Abaqus. The second script is automatically started
at the end of the first script.
Two separate scripts are required, since Abaqus Python has no tool to
install Python libraries which are necessary to create voronoi structures.
Both files need to be located in the same folder to ensure access to the output
data of the first script.
---
Install the following extending python-libraries to make sure, the code
is going to work properly. Required are: matplotlib, numpy, scipy and shapely
(those can be freely downloaded for the preferred python compiler).
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import subprocess
import math
from scipy.spatial import Voronoi
from shapely.geometry import Polygon
from numpy import pi, cos, sin, sqrt


#------------------------------------------------------------------------------------------------------------
#--- HIERARCHICAL LEVEL 3 - MUSCLE FIBER --------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


#------ FUNCTIONS -------------------------------------------------------------------------------------------

def voronoi_finite_polygons_2d(vor, radius=None):

    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions (The function is freely available online).

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.
    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()*2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge
            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


#------------------------------------------------------------------------------------------------------------
#------ START USER INPUT ------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


#------ GENERAL INPUT ---------------------------------------------------------------------------------------

# (1) number of myofibrils
number_myofibrils = 100

# (2) myofibril volume fraction [%] (0<vf_myofibril<100)
vf_myofibril = 80

# (3) sarcoplasmic reticulum volume fraction [%] (0<vf_sarcoplasmic_reticulum<100)
vf_sarcoplasmic_reticulum = 15

# (4) sarcolemma volume fraction [%] (0<vf_sarcolemma<100)
vf_sarcolemma = 5
# note that all volume fractions must sum to 100%!

# (5) muscle fiber radius
r_muscle_fiber = 30

# (6) muscle fiber length
l_muscle_fiber = 60


#------ ADDITIONAL INPUT FOR MULTIPLE STRUCTURES ------------------------------------------------------------

"""
Multiple structures are generated by changing one input value. First, choose
the changing input value (enter 0 to change nothing, enter 1 to change
the number of myofibrils, enter 2 to change the myofibril volume fraction,
enter 3 to change the sarcoplasmic reticulum volume fraction, enter 4 to
change the sarcolemma volume fraction, enter 5 to change the muscle fiber
radius or enter 6 to change the muscle fiber length). After choosing the number
of structures, a minimum and maximum value for the changing input value needs
to be defined. These values replace the upper input value. By the desired number
of structures, the values of the varying parameter are automatically calculated
for the structures and are evenly distributed in the desired range.
"""

# choose the changing input value (enter 0 if none, otherwise enter 1, 2, 3, 4, 5 or 6)
vary = 0

# number of structures
number_structures = 1

# enter the minimum and maximum value for the changing input value
min_value = 10
max_value = 30


#------ OPTIONAL INPUT --------------------------------------------------------------------------------------

"""
This factor is additionally used for the sunflower seed arrangement to optimize
the seed distribution in edge regions for voronoi tessellation (higher values
lead to more circular seed arrangement in edge regions but can also result in a
more edge-oriented seed arrangement in terms of all seeds; alpha=1 is
recommended!).
"""

alpha = 1


#------------------------------------------------------------------------------------------------------------
#------ END USER INPUT --------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


#------ TEXT FILES ------------------------------------------------------------------------------------------

"""
Create two text files where the names of all text files for (1) the input variables
and (2) the vertice coordinates are saved. Both text files are used for the second
script to create all structures.
"""

all_txt_input = "all_txt_input.txt"
all_txt_coordinates = "all_txt_coordinates.txt"

open(all_txt_input, 'w').close()
open(all_txt_coordinates, 'w').close()


#------ GENERATION OF THE STRUCTURES ------------------------------------------------------------------------

# sort input values and their labeling in lists
input_list = [0 , number_myofibrils, vf_myofibril, vf_sarcoplasmic_reticulum,
              vf_sarcolemma, r_muscle_fiber, l_muscle_fiber]
input_list_string = ['', '_number_myofibrils_', '_vf_myofibril_', '_vf_sarc_ret_',
                     '_vf_sarcolemma_', '_radius_', '_length_']

# check the type of the input values
if (type(vary) == str or type(number_structures) == str or type(min_value) == str or type(max_value) == str or
type(number_myofibrils) == str or type(vf_myofibril) == str or type(vf_sarcoplasmic_reticulum) == str or
type(vf_sarcolemma) == str or type(r_muscle_fiber) == str or type(l_muscle_fiber) == str or type(alpha) == str):
    sys.exit('+ + + INPUT ERROR - NO STRINGS + + +')

# check the input for multiple structures and stop the programm if the input values are incorrect
if not 0<=vary<=6:
    sys.exit('+ + + INPUT ERROR - VARY + + +')
if number_structures <=0:
    sys.exit('+ + + INPUT ERROR - NUMBER OF STRUCTURES + + +')
if vary == 0:
    number_structures = 1
if number_structures == 1:
    vary = 0
if number_structures>1:
    input_list[vary] = min_value

# loop to create all structures
for var in np.linspace(input_list[vary], max_value, num=number_structures):

    if vary == 1:
        input_list[vary] = round(var)
    else:
        input_list[vary] = round(var,4)

    # assign the input values for each structure.
    # the ratio between the non-varying volume fractions remains the same.
    if vary == 0:
        vf_myofibril = input_list[2]/100
        vf_sarcoplasmic_reticulum = input_list[3]/100
        vf_sarcolemma = input_list[4]/100
    elif vary == 1:
        number_myofibrils = input_list[1]
        vf_myofibril = input_list[2]/100
        vf_sarcoplasmic_reticulum = input_list[3]/100
        vf_sarcolemma = input_list[4]/100
    elif vary == 2:
        vf_myofibril = input_list[2]/100
        vf_sarcoplasmic_reticulum = input_list[3]/(input_list[3]+input_list[4])*(100-input_list[2])/100
        vf_sarcolemma = 1 - (vf_myofibril + vf_sarcoplasmic_reticulum)
    elif vary == 3:
        vf_sarcoplasmic_reticulum = input_list[3]/100
        vf_myofibril = input_list[2]/(input_list[2]+input_list[4])*(100-input_list[3])/100
        vf_sarcolemma = 1 - (vf_sarcoplasmic_reticulum + vf_myofibril)
    elif vary == 4:
        vf_sarcolemma = input_list[4]/100
        vf_myofibril = input_list[2]/(input_list[2]+input_list[3])*(100-input_list[4])/100
        vf_sarcoplasmic_reticulum = 1 - (vf_sarcolemma + vf_myofibril)
    elif vary == 5:
        r_muscle_fiber = input_list[5]
        vf_myofibril = input_list[2]/100
        vf_sarcoplasmic_reticulum = input_list[3]/100
        vf_sarcolemma = input_list[4]/100
    elif vary ==6:
        l_muscle_fiber = input_list[6]
        vf_myofibril = input_list[2]/100
        vf_sarcoplasmic_reticulum = input_list[3]/100
        vf_sarcolemma = input_list[4]/100

    # check the input values and stop the programm if the input values are incorrect
    if number_myofibrils <= 2:
        sys.exit('+ + + INPUT ERROR - NUMBER OF MYOFIBRILS + + +')
    if not (0<vf_myofibril<1):
        sys.exit('+ + + INPUT ERROR - VOLUME FRACTION OF MYOFIBRIL + + +')
    if not (0<vf_sarcoplasmic_reticulum<1):
        sys.exit('+ + + INPUT ERROR - VOLUME FRACTION OF SARCOPLASMIC RETICULUM + + +')
    if not (0<vf_sarcolemma<1):
        sys.exit('+ + + INPUT ERROR - VOLUME FRACTION OF SARCOLEMMA + + +')
    if not math.isclose(vf_myofibril+vf_sarcoplasmic_reticulum+vf_sarcolemma, 1):
        sys.exit("+ + + INPUT ERROR - SUM OF VOLUME FRACTIONS MUST BE 100%")
    if r_muscle_fiber <= 0:
        sys.exit('+ + + INPUT ERROR - RADIUS OF MUSCLE FIBER + + +')
    if l_muscle_fiber <= 0:
        sys.exit('+ + + INPUT ERROR - LENGTH OF MUSCLE FIBER + + +')

    # define name of txt-files for the vertice coordinates and input values for each structure
    if vary == 0:
        filename = "coordinates_musclefiber.txt"
        filename2 = "input_values_musclefiber.txt"
    else:
        filename = "coordinates_musclefiber"+input_list_string[vary]+str(input_list[vary])+".txt"
        filename2 = "input_values_musclefiber"+input_list_string[vary]+str(input_list[vary])+".txt"


    #------ ADJUST SIZES TO INPUT VOLUME FRACTIONS ----------------------------------------------------------

    """
    First, the structure is considered without the sarcolemma. For this, the radius of
    the muscle fiber is redetermined and the myofibril volume fraction has to be
    adjusted as well. The sarcolemma is inserted afterwards in script 2, which
    generates the structure with the required input values.
    """

    # determine sarcolemma thickness
    thickness_sarcolemma = r_muscle_fiber*(1-sqrt(1-vf_sarcolemma))

    # redefined muscle fiber radius without sarcolemma
    r_muscle_fiber = r_muscle_fiber-thickness_sarcolemma

    # redefined myofibril volume fraction adjusted to the new muscle fiber radius
    vf_myofibril = vf_myofibril/(1-vf_sarcolemma)


    #------ SEED DISTRIBUTION FOR VORONOI TESSELLATION ------------------------------------------------------

    """
    The sunflower seed arrangement is used to define the seeds for voronoi
    tessellation. This procedure allows to generate homogeneous circular voronoi
    diagrams due to the uniform seed distribution in circular structures. The
    number of myofibrils per muscle fiber should be >10 for improved seed
    distribution, whereas, in reality, the number is strongly exceeded. Additionally,
    the factor alpha is applied to optimize the distribution of the seeds at the
    edge.
    """

    n = number_myofibrils + 1
    s = (1,n-1)
    x = np.zeros(s)
    y = np.zeros(s)

    b = round(alpha*sqrt(n))
    phi = (sqrt(5)+1)/2
    for k in range(0,n):
        if k > 1+(n-b):
            r = 1
        else:
            r = sqrt(k/2)/sqrt((n-(b+1))/2)
        if r != 0:
            theta = 2*pi*k/phi**2
            x[0,k-1] = r*cos(theta)
            y[0,k-1] = r*sin(theta)

    points_x = np.multiply(x,r_muscle_fiber)
    points_y = np.multiply(y,r_muscle_fiber)
    points = (np.vstack((points_x,points_y))).T


    #------ VORONOI TESSELLATION ----------------------------------------------------------------------------

    """
    The voronoi diagram is generated without considering the sarcoplasmic
    reticulum and all polygons with the associated vertices are saved in a list
    whereas the coordinates of the vertices are written sorted into a separate
    array. This ensures individual access to each polygon.
    ---
    The function 'voronoi_finite_polygons_2d' is used to reconstruct infinite
    Voronoi structures to a finite 2D diagram. Since the diagram is finite now,
    but still undefined, a square box is created to define the rectangular boundary
    of the diagram and a circle is generated to represent the muscle fiber boundary.
    The box and the circle are used in "FINAL VERTICE CALCULATION".
    """

    vor = Voronoi(points)

    # create finite polygons with undefined diagram size from infinite voronoi tessellation
    regions, vertices = voronoi_finite_polygons_2d(vor)

    # generate a box to define the boundary of the voronoi structure
    box_size = r_muscle_fiber*1.3
    box = Polygon([[-box_size, -box_size], [-box_size, box_size], [box_size, box_size], [box_size, -box_size]])

    # generate a circle consisting of a polygon with many vertices to define the outer muscle fiber boundary
    step_range = np.arange(0, 2*pi, 0.1*pi/number_myofibrils)
    outer_poly = [[]]*len(step_range)
    for i in range(len(step_range)):
        outer_poly[i] = (r_muscle_fiber*cos(step_range[i]), r_muscle_fiber*sin(step_range[i]))
    outer_boundary = Polygon(outer_poly)


    #------ FINAL VERTICE CALCULATION -----------------------------------------------------------------------

    """
    The following part scales down all voronoi polygons in order to insert
    sarcoplasmic reticulum between the myofibrils. The new coordinates of the
    vertices are written in a txt-file.
    ---
    Reducing the size by scaling down each polygon to create sarcoplasmic
    reticulum between myofibrils results in an outer edge with relatively small
    sarcoplasmic reticulum width. A second scaling of the entire structure is
    done with an experimentally determined edge-width factor which improves the
    thickness ratio of the sarcoplasmic reticulum between the inner part and the
    edge of the muscle fiber. This slightly changes the myofibril volume fraction
    of the structure. Therefore, the additionally defined rescale factor adjusts
    the polygons to the required myofibril volume fraction afterwards. The
    coordinates of their vertices are saved in a text-file.
    ---
    Cutting the round muscle fiber from the rectangular voronoi diagram is done in
    the second script.
    """

    area_myofibrils = 0
    # define edge-width factor
    edge_fac = (1 + 0.3*(1-vf_myofibril)*(1-4e-3*number_myofibrils))

    with open(filename, "w") as c:
        for x in range(2):
            counter = 0
            # process polygon-wise
            for region in regions:
                # step 1: the previously defined box is used to cut off unused outer regions
                polygon = vertices[region]
                poly = Polygon(polygon)
                poly = poly.intersection(box)
                polygon = np.array([p for p in poly.exterior.coords])

                # step 2: scale polygons according to required myofibril volume fraction
                if x == 0:
                    for i in range(len(polygon)):
                        polygon[i] = polygon[i]*sqrt(vf_myofibril*edge_fac) + (1-sqrt(vf_myofibril*edge_fac))*points[counter]

                    # calculate area of newly scaled polygons
                    poly_scaled = Polygon(polygon)
                    poly_scaled = poly_scaled.intersection(outer_boundary)
                    area_myofibrils += poly_scaled.area

                # step 3: 2. scaling by using the edge-width factor; the rescale factor adjusts the polygons to the required volume fraction
                else:
                    for i in range(len(polygon)):
                        polygon[i] = polygon[i]*sqrt(vf_myofibril*edge_fac) + (1-sqrt(vf_myofibril*edge_fac))*points[counter]

                        polygon[i] = polygon[i]*sqrt(rescale_factor) + (1-sqrt(rescale_factor))*np.array((0, 0))

                        # write coordinates to txt-file
                        c.write(str(polygon[i][0]) + " ")
                        c.write(str(polygon[i][1]) + "\n")

                    plt.fill(*zip(*polygon), alpha=0.5)

                    c.write("\n")

                counter += 1
            # define rescale factor
            rescale_factor = (r_muscle_fiber**2*pi*vf_myofibril)/area_myofibrils
        # determine the radius of the muscle fiber without the surrounding sarcoplasmic reticulum at the muscle fiber edge
        r_muscle_fiber_inner = r_muscle_fiber*sqrt(rescale_factor)
    c.close()


    #------ PLOT RESULTS ------------------------------------------------------------------------------------

    """
    Preview of the voronoi structure in a colored 2D plot.
    """

#    plt.plot()
#    plt.axis('equal')
#    plt.xlim(-r_muscle_fiber*1.3, r_muscle_fiber*1.3)
#    plt.ylim(-r_muscle_fiber*1.3, r_muscle_fiber*1.3)
#    plt.show()


    #------ WRITE TEXT FILES --------------------------------------------------------------------------------

    """
    Add the input values and text file names of the vertice coordinates and
    input values to the corresponding text files. These files are used in the
    second script.
    """

    with open(filename2, "w") as i:
        i.write(str(number_myofibrils) + "\n")
        i.write(str(r_muscle_fiber) + "\n")
        i.write(str(r_muscle_fiber_inner) + "\n")
        i.write(str(l_muscle_fiber) + "\n")
        i.write(str(thickness_sarcolemma) + "\n")
    i.close()

    with open(all_txt_input, "a") as i:
        i.write(str(filename2) + "\n")
    i.close

    with open(all_txt_coordinates, "a") as i:
        i.write(str(filename) + "\n")
    i.close


#------ START THE SECOND SCRIPT -----------------------------------------------------------------------------

"""
A batch-file is created which is used to start the second python script.
"""

bat = open('musclefiber.bat','w+')
bat.write('''@echo off
          %start abaqus cae script=musclefiber2.py
          ''')
bat.close()

subprocess.call(['musclefiber.bat'])
