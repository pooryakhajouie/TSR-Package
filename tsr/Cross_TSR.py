# Program to calculate the cross keys between two chains of one proteins
# You should have drug_atom_lexical_txt.csv in you directory to be able to run the code
# The code can operate on single PDB, multiple PDB as a list, or a csv file which has the protein and two corresponding chains

import os
import csv
import pandas as pd
import Bio.PDB
import numpy as np
import math
from joblib import Parallel, delayed
import warnings
import multiprocessing

warnings.filterwarnings("ignore")

dTheta = 29
dLen = 17
numOfLabels = 112
num_cores = multiprocessing.cpu_count()

# Read atom sequences from a CSV file
def load_atom_seq(filepath):
    atomSeq = {}
    with open(filepath, 'r') as atomSeqFile:
        reader = csv.reader(atomSeqFile)
        next(reader)
        for row in reader:
            atomSeq[row[2]] = row[1]
    return atomSeq

# Theta Bin for 3D
def thetaClass_(Theta):
    theta_bins = [12.11, 17.32, 21.53, 25.21, 28.54, 31.64, 34.55, 37.34, 40.03, 42.64,
                  45.17, 47.64, 50.05, 52.43, 54.77, 57.08, 59.38, 61.64, 63.87, 66.09,
                  68.30, 70.5, 72.69, 79.2, 81.36, 83.51, 85.67, 87.80, 90.00]
    return next((i+1 for i, val in enumerate(theta_bins) if Theta < val), len(theta_bins))

# MaxDist bin for 3D
def dist12Class_(dist12):
    dist_bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1000]
    return next((i + 1 for i, val in enumerate(dist_bins) if dist12 < val), len(dist_bins))


# Distance calculation function
def calDist(x1, y1, z1, x2, y2, z2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)

def remove_hydrogen_atoms(data_dir, pdb_file):
    input_path = f"{data_dir}/{pdb_file}.pdb"
    output_path = f"{data_dir}/no_hydrogen/{pdb_file}.pdb"

    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Read and filter lines without hydrogen atoms
    with open(input_path, 'r') as input_file, open(output_path, 'w') as output_file:
        for line in input_file:
            atom_type = line[0:6].strip()
            if (atom_type != 'ATOM' and atom_type != 'HETATM'):
                output_file.write(line)
            elif (atom_type == 'ATOM' and line[75:78].strip() != 'H'):
                output_file.write(line)
            elif (atom_type == 'HETATM' and line[75:78].strip() != 'H'):
                output_file.write(line)
    return output_path

# Process individual PDB file
def process_pdb_file(data_dir, pdb_file, chain1, chain2, output_subdir, output_option='both'):
    no_hydrogen_file = remove_hydrogen_atoms(data_dir, pdb_file)
    pdb_parser = Bio.PDB.PDBParser()
    Structure = pdb_parser.get_structure("PrimaryStructureChain", no_hydrogen_file)
    model = Structure[0]
    atomSeq = load_atom_seq("drug_atom_lexical_txt.csv")

    Chain1_At_Coord_Dict = {}
    Chain2_At_Coord_Dict = {}
    Chain1_Res_Dict = {}
    Chain2_Res_Dict = {}
    Chain1_Atom_Dict = {}
    Chain2_Atom_Dict = {}
    Chain1_Res_Id_Dict = {}
    Chain2_Res_Id_Dict = {}

    keyDict3D = {}

    # Create output files based on output_options
    if output_option in ['triplets', 'both']:
        outputFile1 = open(f'{data_dir}/{output_subdir}/{os.path.basename(pdb_file)}_{chain1}_{chain2}.cross_triplets_29_17', 'w')
        outputFile1.writelines('Residue1   Residue2   Residue3   Edge1  Edge2  Edge3\t   Coor_R1\t           Coor_R2\t         CoorR3\t         Theta\tmax_dist\td_3\tkey3D label1 label2 label3 BinLength BinTheta d12 d13\n')

    if output_option in ['keys', 'both']:
        outputFile2 = open(f'{data_dir}/{output_subdir}/{os.path.basename(pdb_file)}_{chain1}_{chain2}.cross_keys_29_17', 'w')
        outputFile2.writelines('key\t\tfreq\n')

    Chain1_counter = Chain2_counter = 0

    for chain in model:
        if chain.id == chain1:
            for res in chain:
                if str(res)[17:18] != 'H' and str(res)[17:18] != 'W':
                    Res_Name = str(res)[9:12]
                    numeric_filter = filter(str.isdigit, str(res.id)[6:10])
                    Res_Id = "".join(numeric_filter)
                    for atom in res:
                        Chain1_At_Coord_Dict[Chain1_counter] = atom.get_vector()
                        Chain1_Res_Dict[Chain1_counter] = Res_Name
                        Chain1_Res_Id_Dict[Chain1_counter] = Res_Id
                        Chain1_Atom_Dict[Chain1_counter] = atom.get_name()
                        Chain1_counter += 1

        elif chain.id == chain2:
            for res in chain:
                if str(res)[17:18] != 'H' and str(res)[17:18] != 'W':
                    Res_Name = str(res)[9:12]
                    numeric_filter = filter(str.isdigit, str(res.id)[6:10])
                    Res_Id = "".join(numeric_filter)
                    for atom in res:
                        Chain2_At_Coord_Dict[Chain2_counter] = atom.get_vector()
                        Chain2_Res_Dict[Chain2_counter] = Res_Name
                        Chain2_Res_Id_Dict[Chain2_counter] = Res_Id
                        Chain2_Atom_Dict[Chain2_counter] = atom.get_name()
                        Chain2_counter += 1

    InterRes_Chain1 = []  # Interacting residue list- Chain1
    InterRes_Chain2 = []  # Interacting residue list -Chain2

    chain1_Inter_Atom_Pos = []  # Interacting atom Pos -Chain1
    chain2_Inter_Atom_Pos = []  # Interacting atom pos -chain2

    for i in range(len(Chain1_At_Coord_Dict)):
        for j in range(len(Chain2_At_Coord_Dict)):
            X1 = Chain1_At_Coord_Dict[i][0]
            Y1 = Chain1_At_Coord_Dict[i][1]
            Z1 = Chain1_At_Coord_Dict[i][2]

            X2 = Chain2_At_Coord_Dict[j][0]
            Y2 = Chain2_At_Coord_Dict[j][1]
            Z2 = Chain2_At_Coord_Dict[j][2]

            Distance = calDist(X1, Y1, Z1, X2, Y2, Z2)
            if Distance <= 5:
                data_Chain1 = f'{Chain1_Res_Dict[i]}_{Chain1_Res_Id_Dict[i]}'
                if data_Chain1 not in InterRes_Chain1:
                    InterRes_Chain1.append(data_Chain1)
                if i not in chain1_Inter_Atom_Pos:
                    chain1_Inter_Atom_Pos.append(i)  # gets the interacting atom position no

                data_Chain2 = f'{Chain2_Res_Dict[j]}_{Chain2_Res_Id_Dict[j]}'
                if data_Chain2 not in InterRes_Chain2:
                    InterRes_Chain2.append(data_Chain2)
                if j not in chain2_Inter_Atom_Pos:
                    chain2_Inter_Atom_Pos.append(j)  # gets the Non-Interacting atom position no

    Antib_Num = len(chain1_Inter_Atom_Pos)
    Spike_Num = len(chain2_Inter_Atom_Pos)

    # Calculates keys between two spike Protein  atom  and one antibody atom
    c = 0
    for i in range(Antib_Num):
        c = Spike_Num - 1
        l = 1
        for j in range(Spike_Num - 1):

            for k in range(c):
                L1 = calDist(Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[i]][0],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[i]][1],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[i]][2],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[j]][0],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[j]][1],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[j]][2])
                L2 = calDist(Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[j]][0],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[j]][1],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[j]][2],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[k + l]][0]
                             , Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[k + l]][1],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[k + l]][2])
                L3 = calDist(Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[i]][0],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[i]][1],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[i]][2],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[k + l]][0],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[k + l]][1],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[k + l]][2])

                l1 = atomSeq[Chain2_Atom_Dict[chain2_Inter_Atom_Pos[j]]]

                l2 = atomSeq[Chain2_Atom_Dict[chain2_Inter_Atom_Pos[k+l]]]
                l3 = atomSeq[Chain1_Atom_Dict[chain1_Inter_Atom_Pos[i]]]


                Med1 = (1 / 2) * math.sqrt(2 * (L1 ** 2) + 2 * (L2 ** 2) - L3 ** 2)
                Med2 = (1 / 2) * math.sqrt(2 * (L2 ** 2) + 2 * (L3 ** 2) - L1 ** 2)
                Med3 = (1 / 2) * math.sqrt(2 * (L3 ** 2) + 2 * (L1 ** 2) - L2 ** 2)
                Median = [Med1, Med2, Med3]
                Label = [l1, l2, l3]
                index1 = [L3, L1, L2]
                # 1st condition
                if l1 != l2 != l3:
                    X = [l1, l2, l3]
                    b3 = Median[Label.index(min(l1, l2, l3))]
                    d12 = index1[Label.index(min(l1, l2, l3))]
                    if d12 == L3 and max(l1, l2, l3) == l2:
                        d13 = L2
                    elif d12 == L3 and max(l1, l2, l3) == l3:
                        d13 = L1

                    elif d12 == L2 and max(l1, l2, l3) == l1:
                        d13 = L1
                    elif d12 == L2 and max(l1, l2, l3) == l2:
                        d13 = L3
                    elif d12 == L1 and max(l1, l2, l3) == l1:
                        d13 = L2
                    elif d12 == L1 and max(l1, l2, l3) == l3:
                        d13 = L3
                    X.remove(max(X))
                    X.remove(min(X))
                    Label1 = max(l1, l2, l3)
                    Label2 = X[0]
                    Label3 = min(l1, l2, l3)

                # 2nd condition
                elif l1 > l2 == l3:
                    Label1 = l1
                    if L2 > L1:
                        b3 = Med3
                        d13 = L1
                        d12 = L2
                        Label2 = l2
                        Label3 = l3
                    else:
                        b3 = Med2
                        d13 = L2
                        d12 = L1
                        Label2 = l3
                        Label3 = l2

                elif l2 > l1 == l3:
                    Label1 = l2
                    if L3 > L2:
                        b3 = Med1
                        d13 = L2
                        d12 = L3
                        Label2 = l3
                        Label3 = l1
                    else:
                        b3 = Med3
                        d13 = L3
                        d12 = L2
                        Label2 = l1
                        Label3 = l3

                elif l3 > l1 == l2:
                    Label1 = l3
                    if L1 > L3:
                        b3 = Med2
                        d13 = L3
                        d12 = L1
                        Label2 = l1
                        Label3 = l2
                    else:
                        b3 = Med1
                        d13 = L1
                        d12 = L3
                        Label2 = l2
                        Label3 = l1
                # 3rd condition
                elif l1 == l2 > l3:
                    b3 = Med3
                    Label3 = l3
                    if L1 > L3:
                        d13 = L1
                        d12 = L2
                        Label1 = l1
                        Label2 = l2
                    else:
                        d13 = L3
                        d12 = L2
                        Label1 = l2
                        Label2 = l1

                elif l1 == l3 > l2:
                    Label3 = l2
                    b3 = Med2
                    if L2 > L3:
                        d13 = L2
                        d12 = L1
                        Label1 = l1
                        Label2 = l3
                    else:
                        d13 = L3
                        d12 = L1
                        Label1 = l3
                        Label2 = l1
                elif l2 == l3 > l1:
                    Label3 = l1
                    b3 = Med1
                    if L2 > L1:
                        d13 = L2
                        d12 = L3
                        Label1 = l2
                        Label2 = l3
                    else:
                        d13 = L1
                        d12 = L3
                        Label1 = l3
                        Label2 = l2

                # 4th condition
                if l1 == l2 == l3:
                    if L2 >= max(L1, L2, L3):
                        b3 = Med3
                        d13 = L1
                        d12 = L2
                        Label1 = l1
                        Label2 = l2
                        Label3 = l3
                    if L1 >= max(L1, L2, L3):
                        b3 = Med2
                        d13 = L2
                        d12 = L1
                        Label1 = l1
                        Label2 = l3
                        Label3 = l2

                    if L3 >= max(L1, L2, L3):
                        b3 = Med1
                        d13 = L1
                        d12 = L3
                        Label1 = l3
                        Label2 = l2
                        Label3 = l1

                a = (d13 ** 2 - (d12 / 2) ** 2 - b3 ** 2)
                b = (2 * (d12 / 2) * b3)
                Theta1 = (math.acos(a / b)) * (180 / math.pi)

                if Theta1 <= 90:
                    Theta = Theta1
                else:
                    Theta = abs(180 - Theta1)
                maxDist = max(L1, L2, L3)
                ClassT1 = thetaClass_(Theta)
                ClassL1 = dist12Class_(maxDist)

                # GENERATE 3D key
                key3D = dLen * dTheta * (numOfLabels ** 2) * (int(Label1) - 1) + \
                        dLen * dTheta * (numOfLabels) * (int(Label2) - 1) + \
                        dLen * dTheta * (int(Label3) - 1) + \
                        dTheta * (ClassL1 - 1) + \
                        (ClassT1 - 1)
                if key3D in keyDict3D:
                    keyDict3D[key3D] += 1
                else:
                    keyDict3D[key3D] = 1

                if output_option in ['triplets', 'both']:
                  outputFile1.write("{}_{}_{}_{}  ".format(Chain1_Res_Dict[chain1_Inter_Atom_Pos[i]], chain1,
                                                          Chain1_Res_Id_Dict[chain1_Inter_Atom_Pos[i]],
                                                          Chain1_Atom_Dict[chain1_Inter_Atom_Pos[i]]))
                  outputFile1.write("{}_{}_{}_{}  ".format(Chain2_Res_Dict[chain2_Inter_Atom_Pos[j]], chain2,
                                                          Chain2_Res_Id_Dict[chain2_Inter_Atom_Pos[j]],
                                                          Chain2_Atom_Dict[chain2_Inter_Atom_Pos[j]]))
                  outputFile1.write("{}_{}_{}_{}  ".format(Chain2_Res_Dict[chain2_Inter_Atom_Pos[k + l]], chain2,
                                                          Chain2_Res_Id_Dict[chain2_Inter_Atom_Pos[k + l]],
                                                          Chain2_Atom_Dict[chain2_Inter_Atom_Pos[k + l]]))
                  outputFile1.write(" {:.2f}  {:.2f}  {:.2f} ".format(L1, L2, L3))
                  outputFile1.write("  {:.2f},{:.2f},{:.2f}   {:.2f},{:.2f},{:.2f}  ".format(
                      Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[i]][0],
                      Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[i]][1],
                      Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[i]][2],
                      Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[j]][0],
                      Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[j]][1],
                      Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[j]][2]))
                  outputFile1.write(
                      " {:.2f},{:.2f},{:.2f}  ".format(Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[k + l]][0],
                                                      Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[k + l]][1],
                                                      Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[k + l]][2]))

                  outputFile1.write("{:.2f}   {:.2f}   {:.2f}   {:.0f} {} {} {} {} {} {} {}\n".format(Theta, maxDist, b3, key3D,
                                                                  int(Label1), int(Label2),
                                                                  int(Label3), ClassL1, ClassT1,d12,d13))
            l += 1
            c -= 1

    # Calcultes keys between one drug and two proteins
    for i in range(Spike_Num):
        c = Antib_Num - 1
        l = 1
        for j in range(Antib_Num - 1):

            for k in range(c):
                L1 = calDist(Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[i]][0],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[i]][1],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[i]][2],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[j]][0],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[j]][1],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[j]][2])

                L2 = calDist(Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[j]][0],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[j]][1],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[j]][2],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[k + l]][0],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[k + l]][1],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[k + l]][2])

                L3 = calDist(Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[i]][0],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[i]][1],
                             Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[i]][2],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[k + l]][0],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[k + l]][1],
                             Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[k + l]][2])

                l1 = atomSeq[Chain1_Atom_Dict[chain1_Inter_Atom_Pos[j]]]
                l2 = atomSeq[Chain1_Atom_Dict[chain1_Inter_Atom_Pos[k+l]]]
                l3 = atomSeq[Chain2_Atom_Dict[chain2_Inter_Atom_Pos[i]]]

                Med1 = (1 / 2) * math.sqrt(2 * (L1 ** 2) + 2 * (L2 ** 2) - L3 ** 2)
                Med2 = (1 / 2) * math.sqrt(2 * (L2 ** 2) + 2 * (L3 ** 2) - L1 ** 2)
                Med3 = (1 / 2) * math.sqrt(2 * (L3 ** 2) + 2 * (L1 ** 2) - L2 ** 2)
                Median = [Med1, Med2, Med3]
                Label = [l1, l2, l3]
                index1 = [L3, L1, L2]
                # 1st condition
                if l1 != l2 != l3:
                    X = [l1, l2, l3]
                    b3 = Median[Label.index(min(l1, l2, l3))]
                    d12 = index1[Label.index(min(l1, l2, l3))]
                    if d12 == L3 and max(l1, l2, l3) == l2:
                        d13 = L2
                    elif d12 == L3 and max(l1, l2, l3) == l3:
                        d13 = L1

                    elif d12 == L2 and max(l1, l2, l3) == l1:
                        d13 = L1
                    elif d12 == L2 and max(l1, l2, l3) == l2:
                        d13 = L3
                    elif d12 == L1 and max(l1, l2, l3) == l1:
                        d13 = L2
                    elif d12 == L1 and max(l1, l2, l3) == l3:
                        d13 = L3
                    X.remove(max(X))
                    X.remove(min(X))
                    Label1 = max(l1, l2, l3)
                    Label2 = X[0]
                    Label3 = min(l1, l2, l3)

                # 2nd condition
                elif l1 > l2 == l3:
                    Label1 = l1
                    if L2 > L1:
                        b3 = Med3
                        d13 = L1
                        d12 = L2
                        Label2 = l2
                        Label3 = l3
                    else:
                        b3 = Med2
                        d13 = L2
                        d12 = L1
                        Label2 = l3
                        Label3 = l2

                elif l2 > l1 == l3:
                    Label1 = l2
                    if L3 > L2:
                        b3 = Med1
                        d13 = L2
                        d12 = L3
                        Label2 = l3
                        Label3 = l1
                    else:
                        b3 = Med3
                        d13 = L3
                        d12 = L2
                        Label2 = l1
                        Label3 = l3

                elif l3 > l1 == l2:
                    Label1 = l3
                    if L1 > L3:
                        b3 = Med2
                        d13 = L3
                        d12 = L1
                        Label2 = l1
                        Label3 = l2
                    else:
                        b3 = Med1
                        d13 = L1
                        d12 = L3
                        Label2 = l2
                        Label3 = l1
                # 3rd condition
                elif l1 == l2 > l3:
                    b3 = Med3
                    Label3 = l3
                    if L1 > L3:
                        d13 = L1
                        d12 = L2
                        Label1 = l1
                        Label2 = l2
                    else:
                        d13 = L3
                        d12 = L2
                        Label1 = l2
                        Label2 = l1

                elif l1 == l3 > l2:
                    Label3 = l2
                    b3 = Med2
                    if L2 > L3:
                        d13 = L2
                        d12 = L1
                        Label1 = l1
                        Label2 = l3
                    else:
                        d13 = L3
                        d12 = L1
                        Label1 = l3
                        Label2 = l1
                elif l2 == l3 > l1:
                    Label3 = l1
                    b3 = Med1
                    if L2 > L1:
                        d13 = L2
                        d12 = L3
                        Label1 = l2
                        Label2 = l3
                    else:
                        d13 = L1
                        d12 = L3
                        Label1 = l3
                        Label2 = l2

                # 4th condition
                if l1 == l2 == l3:
                    if L2 >= max(L1, L2, L3):
                        b3 = Med3
                        d13 = L1
                        d12 = L2
                        Label1 = l1
                        Label2 = l2
                        Label3 = l3
                    if L1 >= max(L1, L2, L3):
                        b3 = Med2
                        d13 = L2
                        d12 = L1
                        Label1 = l1
                        Label2 = l3
                        Label3 = l2

                    if L3 >= max(L1, L2, L3):
                        b3 = Med1
                        d13 = L1
                        d12 = L3
                        Label1 = l3
                        Label2 = l2
                        Label3 = l1

                a = (d13 ** 2 - (d12 / 2) ** 2 - b3 ** 2)
                b = (2 * (d12 / 2) * b3)
                Theta1 = (math.acos(a / b)) * (180 / math.pi)

                if Theta1 <= 90:
                    Theta = Theta1
                else:
                    Theta = abs(180 - Theta1)
                maxDist = max(L1, L2, L3)
                ClassT1 = thetaClass_(Theta)
                ClassL1 = dist12Class_(maxDist)

                # Generates 3D key
                key3D = dLen * dTheta * (numOfLabels ** 2) * (int(Label1) - 1) + \
                        dLen * dTheta * (numOfLabels) * (int(Label2) - 1) + \
                        dLen * dTheta * (int(Label3) - 1) + \
                        dTheta * (ClassL1 - 1) + \
                        (ClassT1 - 1)
                if key3D in keyDict3D:
                    keyDict3D[key3D] += 1
                else:
                    keyDict3D[key3D] = 1

                if output_option in ['triplets', 'both']:
                  outputFile1.write("{}_{}_{}_{}  ".format(Chain2_Res_Dict[chain2_Inter_Atom_Pos[i]], chain2,
                                                          Chain2_Res_Id_Dict[chain2_Inter_Atom_Pos[i]],
                                                          Chain2_Atom_Dict[chain2_Inter_Atom_Pos[i]]))
                  outputFile1.write("{}_{}_{}_{}  ".format(Chain1_Res_Dict[chain1_Inter_Atom_Pos[j]], chain1,
                                                          Chain1_Res_Id_Dict[chain1_Inter_Atom_Pos[j]],
                                                          Chain1_Atom_Dict[chain1_Inter_Atom_Pos[j]]))
                  outputFile1.write("{}_{}_{}_{}  ".format(Chain1_Res_Dict[chain1_Inter_Atom_Pos[k + l]], chain1,
                                                          Chain1_Res_Id_Dict[chain1_Inter_Atom_Pos[k + l]],
                                                          Chain1_Atom_Dict[chain1_Inter_Atom_Pos[k + l]]))
                  outputFile1.write(" {:.2f}  {:.2f}  {:.2f} ".format(L1, L2, L3))
                  outputFile1.write(" {:.2f},{:.2f},{:.2f}   {:.2f},{:.2f},{:.2f}  ".format(
                      Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[i]][0],
                      Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[i]][1],
                      Chain2_At_Coord_Dict[chain2_Inter_Atom_Pos[i]][2],
                      Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[j]][0],
                      Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[j]][1],
                      Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[j]][2]))
                  outputFile1.write(
                      " {:.2f},{:.2f},{:.2f}  ".format(Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[k + l]][0],
                                                      Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[k + l]][1],
                                                      Chain1_At_Coord_Dict[chain1_Inter_Atom_Pos[k + l]][2]))

                  outputFile1.write( "{:.2f}   {:.2f}   {:.2f}   {:.0f} {} {} {} {} {} {} {} \n".format(Theta, maxDist, b3, key3D,
                                                                  int(Label1), int(Label2),
                                                                  int(Label3), ClassL1, ClassT1,d12,d13))
            l += 1
            c -= 1
    if output_option in ['keys', 'both']:
      for value_ in keyDict3D:
          outputFile2.writelines([str(value_), '\t', str(keyDict3D[value_]), '\n'])

    if output_option in ['triplets', 'both']:
        outputFile1.close()

    if output_option in ['keys', 'both']:
        outputFile2.close()

# Main function to handle input and output for CrossTSR
def CrossTSR(data_dir, input_files, chain_1=None, chain_2=None, output_option='both', output_subdir='lexicographic'):
    os.makedirs(os.path.join(data_dir, output_subdir), exist_ok=True)
    chain_dict_1, chain_dict_2 = {}, {}

    # Handle single file input
    if isinstance(input_files, str):
        if input_files.endswith('.csv'):
            # Read CSV
            df = pd.read_csv(input_files)
            chain_dict_1 = dict(zip(df['protein'].str.upper(), df['chain_s']))
            chain_dict_2 = dict(zip(df['protein'].str.upper(), df['chain_h']))
            input_files = df['protein'].str.upper().tolist()
        else:
            input_files = [input_files]

    if chain_1 and chain_2:
        if isinstance(chain_1, list) and isinstance(chain_2, list):
            chain_dict_1 = {f.upper(): c for f, c in zip(input_files, chain_1)}
            chain_dict_2 = {f.upper(): c for f, c in zip(input_files, chain_2)}
        elif isinstance(chain_1, str) and isinstance(chain_2, str):
            chain_dict_1 = {f.upper(): chain_1 for f in input_files}
            chain_dict_2 = {f.upper(): chain_2 for f in input_files}

    # Parallel processing
    Parallel(n_jobs=num_cores, verbose=50)(
        delayed(process_pdb_file)(data_dir, file_name.upper(), chain_dict_1.get(file_name.upper()), chain_dict_2.get(file_name.upper()), output_subdir, output_option)
        for file_name in input_files
    )

'''
Example usage:
data_dir = 'Dataset'
pdb_files = ['6W41', '6XC3']
chain_1=['C', 'C']
chain_2=['H', 'C']
PDB_DL(pdb_files)
CrossTSR(data_dir, pdb_files, chain_1=chain_1, chain_2=chain_2, output_option='both')

csv_file = 'sample_details.csv'
PDB_DL(csv_file)
CrossTSR(data_dir, csv_file, output_option='triplets')'''