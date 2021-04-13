from schrodinger.application.desmond.packages import traj, topo
from schrodinger.structutils.build import delete_hydrogens
from schrodinger.structutils.analyze import center_of_mass
from schrodinger.application.desmond.packages import msys
from mpl_toolkits.mplot3d.axes3d import Axes3D
from schrodinger import structure
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from math import pi
from numpy import *
import numpy as np # ugh
import pickle

#-------------------------------------------------------------------------------
## ANALYSIS TOOLS DUMP FOR DNA AMPHIPHILE CHEM PAPER ##
## many of these functions were not utilized in analysis for the final paper ##
## functions used to fit a circle to backbone atoms were adapted from: ## 
## https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/ ##
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Generate points on circle
# P(t) = r*cos(t)*u + r*sin(t)*(n x u) + C
def generate_circle_by_vectors(t, C, r, n, u):
    n = n/linalg.norm(n)
    u = u/linalg.norm(u)
    P_circle = r*cos(t)[:,newaxis]*u + r*sin(t)[:,newaxis]*cross(n,u) + C
    return P_circle

#-------------------------------------------------------------------------------
# Generate circle by angles
def generate_circle_by_angles(t, C, r, theta, phi):
    # Orthonormal vectors n, u, <n,u>=0
    n = array([cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)])
    u = array([-sin(phi), cos(phi), 0])

    # P(t) = r*cos(t)*u + r*sin(t)*(n x u) + C
    P_circle = r*cos(t)[:,newaxis]*u + r*sin(t)[:,newaxis]*cross(n,u) + C
    return P_circle

#-------------------------------------------------------------------------------
# FIT CIRCLE 2D
# - Find center [xc, yc] and radius r of circle fitting to set of 2D points
# - Optionally specify weights for points
#
# - Implicit circle function:
#   (x-xc)^2 + (y-yc)^2 = r^2
#   (2*xc)*x + (2*yc)*y + (r^2-xc^2-yc^2) = x^2+y^2
#   c[0]*x + c[1]*y + c[2] = x^2+y^2
#
# - Solution by method of least squares:
#   A*c = b, c' = argmin(||A*c - b||^2)
#   A = [x y 1], b = [x^2+y^2]
def fit_circle_2d(x, y, w=[]):

    A = array([x, y, ones(len(x))]).T
    b = x**2 + y**2

    # Modify A,b for weighted least squares
    if len(w) == len(x):
        W = diag(w)
        A = dot(W,A)
        b = dot(W,b)

    # Solve by method of least squares
    c, res = linalg.lstsq(A,b,rcond=None)[0:2]
    # print(c)
    # Get circle parameters from solution c
    xc = c[0]/2
    yc = c[1]/2
    r = sqrt(c[2] + xc**2 + yc**2)
    return xc, yc, r, res

#-------------------------------------------------------------------------------
# RODRIGUES ROTATION
# - Rotate given points based on a starting and ending vector
# - Axis k and angle of rotation theta given by vectors n0,n1
#   P_rot = P*cos(theta) + (k x P)*sin(theta) + k*<k,P>*(1-cos(theta))
def rodrigues_rot(P, n0, n1):

    # If P is only 1d array (coords of single point), fix it to be matrix
    if P.ndim == 1:
        P = P[newaxis,:]

    # Get vector of rotation k and angle theta
    n0 = n0/linalg.norm(n0)
    n1 = n1/linalg.norm(n1)
    k = cross(n0,n1)
    k = k/linalg.norm(k)
    theta = arccos(dot(n0,n1))

    # Compute rotated points
    P_rot = zeros((len(P),3))
    for i in range(len(P)):
        P_rot[i] = P[i]*cos(theta) + cross(k,P[i])*sin(theta) + k*dot(k,P[i])*(1-cos(theta))

    return P_rot

#-------------------------------------------------------------------------------
# ANGLE BETWEEN
# - Get angle between vectors u,v with sign based on plane with unit normal n
def angle_between(u, v, n=None):
    if n is None:
        return arctan2(linalg.norm(cross(u,v)), dot(u,v))
    else:
        return arctan2(dot(n,cross(u,v)), dot(u,v))

#-------------------------------------------------------------------------------
# - Make axes of 3D plot to have equal scales
# - This is a workaround to Matplotlib's set_aspect('equal') and axis('equal')
#   which were not working for 3D
def set_axes_equal_3d(ax):
    limits = array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    spans = abs(limits[:,0] - limits[:,1])
    centers = mean(limits, axis=1)
    radius = 0.5 * max(spans)
    ax.set_xlim3d([centers[0]-radius, centers[0]+radius])
    ax.set_ylim3d([centers[1]-radius, centers[1]+radius])
    ax.set_zlim3d([centers[2]-radius, centers[2]+radius])

#-------------------------------------------------------------------------------
# function for getting distance of a point from the plane
def shortest_distance(x1, y1, z1, a, b, c, d):
    d = abs((a * x1 + b * y1 + c * z1 + d))
    e = (math.sqrt(a * a + b * b + c * c))
    return d/e

#-------------------------------------------------------------------------------
# main function that calls others in order to generate deviation from perfect circle
def get_rmsd(p_atoms, c_atoms=None):
    #-------------------------------------------------------------------------------
    # (1) Get 3D points and create P Nx3 matrix
    #-------------------------------------------------------------------------------
    # center in space
    P = p_atoms
    P = array(P)
    P_mean = P.mean(axis=0)
    P_centered = P - P_mean
    U,s,V = linalg.svd(P_centered)

    # Normal vector of fitting plane is given by 3rd column in V
    # Note linalg.svd returns V^T, so we need to select 3rd row from V^T
    normal = V[2,:]
    d = -dot(P_mean, normal)  # d = -<p,n>

    #-------------------------------------------------------------------------------
    # (2) Project points to coords X-Y in 2D plane
    #-------------------------------------------------------------------------------
    P_xy = rodrigues_rot(P_centered, normal, [0,0,1])
    # C_xy = rodrigues_rot(C_centered, normal, [0,0,1])

    #-------------------------------------------------------------------------------
    # (3) Fit circle in new 2D coords
    #-------------------------------------------------------------------------------
    xc, yc, r, res = fit_circle_2d(P_xy[:,0], P_xy[:,1])

    # get plane for atoms
    # atoms = P_centered
    # com = average(P_centered[:,:], axis=0, weights=P_centered[:,])
    point = array([0, 0, 0])


    # a plane is a*x+b*y+c*z+d=0
    # [a,b,c] is the normal. Thus, we have to calculate
    # d and we're set
    d = -point.dot(normal)

    # create x,y
    xx, yy = meshgrid(range(10), range(10))

    # calculate corresponding z
    z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]

    # iterate through atoms and get the sum of ditances between plane and point
    plane_res = 0
    for atom in P_centered:
        plane_res += shortest_distance(atom[0], atom[1], atom[2], normal[0], normal[1], normal[2], d)**2

    # combine residual from fit circle and from plane of circle
    # if your monomeric units do not have backbones in a plane, you can weight plane_res more heavily by multiplying it
    to_return = sqrt(res + plane_res)
    return to_return

#-------------------------------------------------------------------------------
# moving average of N
def moving_avg(list):
    N = 5
    cumsum, moving_aves = [0], []
    for val in list[:6]:
        moving_aves.append(val)

    for i, x in enumerate(list[5:], 1):
        # print(i, x)
        cumsum.append(cumsum[i-1] + x)
        if i>=N:
            moving_ave = (cumsum[i] - cumsum[i-N])/N
            # print(moving_ave)
            #can do stuff with moving_ave here
            moving_aves.append(moving_ave)

    return moving_aves

#-------------------------------------------------------------------------------
# function to get atom ids
def get_atom_ids(st, ele):
    ids = []
    for atom in st.atom:
        if atom.element == ele:
            ids.append(atom.index)
    return ids

#-------------------------------------------------------------------------------
# function to get terminal and connecting carbon ids for alkyl chains
def get_alkyl_ids(st):
    term_ids = []
    connect_ids = []
    for atom in st.atom:
        if atom.element == "C" and len(atom.bond) == 1:
            term_ids.append(atom.index-1)
            tmp_atom, i = get_alkyl_ids_helper(atom)
            connect_ids.append(tmp_atom.index-1)

    return term_ids, connect_ids, i

#-------------------------------------------------------------------------------
# helper function for finding connected carbon
def get_alkyl_ids_helper(atom):
    # initialize visited atom list and list of atom bonds
    visited = [atom]
    bonds = atom.bond
    # backbone carbon will have three bonds, all others one or two
    i = 0
    while len(bonds) != 3:
        i += 1
        # iterate through bonds and get bonded atom
        for bond in bonds:
            bonded_atom = bond.atom2
            # if we haven't visited bonded atom (this might be redundant) then move to it
            if bonded_atom not in visited:
                visited.append(atom)
                atom = bonded_atom
                bonds = atom.bond
    return atom, (i+2)

#-------------------------------------------------------------------------------
# function for removing protruding chains from data
def remove_protruding(list, dist=10):
    tmp = []
    for i in list:
        if i < dist:
            tmp.append(i)
    return tmp

#-------------------------------------------------------------------------------
# function to get normalized vector
def calc_vector(coords1, coords2):
    vec = np.array(coords1-coords2)
    vec = vec / np.sqrt(np.sum(vec**2))
    return vec

#-------------------------------------------------------------------------------
# FUNCTIONS FOR SLICE-ANGLE CALCULATIONS
#-------------------------------------------------------------------------------
# function used to find terminal c atom
def find_c_index(st, p_index):
    # find p atom
    for atom in st.atom:
        if atom.index == p_index:
            # get bonded atom
            bonded_atoms = atom.bond
            # find correct oxygen
            for bonded_atom in bonded_atoms:
                bonded_atom = bonded_atom.atom2
                # pass correct oxygen to recursive helper function
                if bonded_atom.index < p_index and bonded_atom.element == "O":
                    bonded_atoms = bonded_atom.bond
                    for bonded_atom in bonded_atoms:
                        bonded_atom = bonded_atom.atom2
                        if bonded_atom.element == "C":
                            c_index = find_c_index_helper(st, bonded_atom.index)
    return c_index

#-------------------------------------------------------------------------------
# recursive helper function
def find_c_index_helper(st, atom_index):
    for atom in st.atom:
        if atom.index == atom_index:
            bonded_atoms = atom.bond
            if len(bonded_atoms) == 1:
                for bond in bonded_atoms:
                    to_return = bond.atom1.index
            else:
                for bonded_atom in bonded_atoms:
                    bonded_atom = bonded_atom.atom2
                    if bonded_atom.element == "C" and bonded_atom.index > atom_index:
                        to_return = find_c_index_helper(st, bonded_atom.index)
    return to_return

#-------------------------------------------------------------------------------
# function to get angle between two P backbone atoms and COM of alkyl chains terminal carbons
def get_angles(st):
    # iterate through tails and get P and terminal carbons
    atom_pairs = []
    prev_p_indices = []
    angles = []
    c_coords = []
    p_coords = []
    for atom in st.atom:
        if atom.element == "P" and atom.chain_name == "A" and atom.index not in prev_p_indices:
            pair = []
            # add to pair and prev found list
            pair.append(atom.index)
            prev_p_indices.append(atom.index)

            # find paired c terminal
            c_index = find_c_index(st, atom.index)
            pair.append(c_index)
            atom_pairs.append(pair)

    # iterate through neighboring pairs and get pizza-slice angle
    for i in range(len(atom_pairs)-1):
        coord_set = []
        # define pairs
        pair1 = atom_pairs[i]
        pair2 = atom_pairs[i + 1]

        # get coordinates
        p1_coords = np.asarray(st.atom[pair1[0]].xyz)
        p2_coords = np.asarray(st.atom[pair2[0]].xyz)
        c1_coords = st.atom[pair1[1]].xyz
        c2_coords = st.atom[pair2[1]].xyz
        p_coords.append(p1_coords)
        p_coords.append(p2_coords)
        # coords.append(c1_coords)
        # coords.append(c2_coords)

        c_com_coords = np.asarray([(c1_coords[0] + c2_coords[0])/2, (c1_coords[1] + c2_coords[1])/2, (c1_coords[2] + c2_coords[2])/2])
        c_coords.append(c_com_coords)
        # get vectors
        vec1 = np.subtract(p1_coords, c_com_coords)
        vec2 = np.subtract(p2_coords, c_com_coords)

        # get angle in degrees
        vec_cos = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
        angle = np.arccos(vec_cos)*180/pi
        angles.append(angle)

    coords = [p_coords, c_coords]
    return angles, coords
