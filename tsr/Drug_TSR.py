# Program to calculate the cross keys between a protein and a drug
# You should have drug_atom_lexical_txt.csv in you directory to be able to run the code
# The code can operate on single PDB, multiple PDB as a list, or a csv file which has the protein name, chain, drug name, and drug id

import os
import csv
import pandas as pd
from Bio.PDB import PDBParser
import numpy as np
import math
from collections import defaultdict
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

def process_pdb_file(data_dir, pdb_file, chain_name, drug_name, drug_id, output_subdir, output_option='both'):
    no_hydrogen_file = remove_hydrogen_atoms(data_dir, pdb_file)
    pdb_parser = PDBParser()
    Structure = pdb_parser.get_structure("PrimaryStructureChain", no_hydrogen_file)
    model = Structure[0]
    atomSeq = load_atom_seq("drug_atom_lexical_txt.csv")

    drugAtom, xCordDrug, yCordDrug, zCordDrug = {}, {}, {}, {}
    protAtom, prot_res_No, prot_res, xCordProt, yCordProt, zCordProt = {}, {}, {}, {}, {}, {}
    keyDict3D = {}

    if output_option in ['triplets', 'both']:
        outputFile1 = open(f'{data_dir}/{output_subdir}/{os.path.basename(pdb_file)}_{chain_name}.drug_triplets_29_17', 'w')
        outputFile1.writelines('Residue1   Residue2   Residue3   Edge1  Edge2  Edge3\t   Coor_R1\t           Coor_R2\t         CoorR3\tTheta\tmax_dist\td_3\tkey3D\n')
    if output_option in ['keys', 'both']:
        outputFile2 = open(f'{data_dir}/{output_subdir}/{os.path.basename(pdb_file)}_{chain_name}.drug_keys_29_17', 'w')
        outputFile2.writelines('key\t\tfreq\n')

    counter2 = 0
    counter1 = 0

    for chain in model:
        if chain.id == chain_name:
            for residue in chain:
                drug_Identifier1 = str(residue)[17:18].strip()
                drug_Identifier2 = str(residue)[16:17].strip()
                resName = str(residue)[9:12].strip()  # residue Name

                numeric_filter = filter(str.isdigit, str(residue.id))
                Res_Id = "".join(numeric_filter)  # Residue ID

                if drug_Identifier1 == 'H' and resName == drug_name and Res_Id == drug_id:

                    for atom1 in residue:
                        atomCoord = atom1.get_vector()
                        drugAtom[counter2] = atom1.get_name()
                        xCordDrug[counter2] = atomCoord[0]
                        yCordDrug[counter2] = atomCoord[1]
                        zCordDrug[counter2] = atomCoord[2]

                        counter2 += 1

                elif drug_Identifier1 != 'H' and drug_Identifier2 != 'H' and drug_Identifier1 != 'W':

                    for atom2 in residue:
                        atomCoord = atom2.get_vector()
                        protAtom[counter1] = atom2.get_name()
                        xCordProt[counter1] = atomCoord[0]
                        yCordProt[counter1] = atomCoord[1]
                        zCordProt[counter1] = atomCoord[2]

                        prot_res[counter1] = resName
                        prot_res_No[counter1] = Res_Id
                        counter1 += 1
    protRegion = {}
    xCordProt_New = {}
    yCordProt_New = {}
    zCordProt_New = {}

    protAtom_New = {}
    prot_res_No_New = {}
    prot_res_New = {}

    counter3 = 0
    for i in range(len(xCordDrug)):
        for j in range(len(xCordProt)):
            disProtDrug = calDist(xCordDrug[i], yCordDrug[i], zCordDrug[i], xCordProt[j], yCordProt[j],
                                  zCordProt[j])
            if disProtDrug <= 3.5:

                if j not in protRegion:
                    xCordProt_New[counter3] = xCordProt[j]
                    yCordProt_New[counter3] = yCordProt[j]
                    zCordProt_New[counter3] = zCordProt[j]

                    protAtom_New[counter3] = protAtom[j]
                    prot_res_No_New[counter3] = prot_res_No[j]
                    prot_res_New[counter3] = prot_res[j]

                    counter3 += 1
                protRegion[j] = xCordProt[j]

    Protlen = len(xCordProt_New)
    Druglen = len(xCordDrug)

    c = 0
    # Calculates keys between two drugs and one protein
    Total_key=0
    for i in range(Protlen):
        c = Druglen - 1
        l = 1
        for j in range(Druglen - 1):

            for k in range(c):
                L1 = calDist(xCordProt_New[i], yCordProt_New[i], zCordProt_New[i], xCordDrug[j], yCordDrug[j],
                             zCordDrug[j])
                L2 = calDist(xCordDrug[j], yCordDrug[j], zCordDrug[j], xCordDrug[k + l], yCordDrug[k + l],
                             zCordDrug[k + l])
                L3 = calDist(xCordProt_New[i], yCordProt_New[i], zCordProt_New[i], xCordDrug[k + l],
                             yCordDrug[k + l],
                             zCordDrug[k + l])

                l1 = atomSeq[drugAtom[j]]
                l2 = atomSeq[drugAtom[k + l]]
                l3 = atomSeq[protAtom_New[i]]

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
                  outputFile1.write(
                      "{}_{}_{}_{}  ".format(prot_res_New[i], chain_name, prot_res_No_New[i], protAtom_New[i]))
                  outputFile1.write("{}_{}_{}_{}  ".format(drug_name, chain_name, drug_id, drugAtom[j]))
                  outputFile1.write("{}_{}_{}_{}  ".format(drug_name, chain_name, drug_id, drugAtom[k + l]))
                  outputFile1.write(" {:.2f}  {:.2f}  {:.2f} ".format(L1, L2, L3))
                  outputFile1.write(
                      "  {:.2f},{:.2f},{:.2f}   {:.2f},{:.2f},{:.2f}  ".format(xCordProt_New[i], yCordProt_New[i],
                                                                              zCordProt_New[i], xCordDrug[j],
                                                                              yCordDrug[j], zCordDrug[j]))
                  outputFile1.write(
                      " {:.2f},{:.2f},{:.2f}  ".format(xCordDrug[k + l], yCordDrug[k + l], zCordDrug[k + l]))
                  outputFile1.write("{:.2f}   {:.2f}   {:.2f}   {:.0f}\n".format(Theta, maxDist, b3, key3D))
                  Total_key+=1

            l += 1
            c -= 1

    # Calcultes keys between one drug and two proteins
    for i in range(Druglen):
        c = Protlen - 1
        l = 1
        for j in range(Protlen - 1):

            for k in range(c):
                L1 = calDist(xCordDrug[i], yCordDrug[i], zCordDrug[i], xCordProt_New[j], yCordProt_New[j],
                             zCordProt_New[j])
                L2 = calDist(xCordProt_New[j], yCordProt_New[j], zCordProt_New[j], xCordProt_New[k + l],
                             yCordProt_New[k + l],
                             zCordProt_New[k + l])
                L3 = calDist(xCordDrug[i], yCordDrug[i], zCordDrug[i], xCordProt_New[k + l], yCordProt_New[k + l],
                             zCordProt_New[k + l])

                l1 = atomSeq[protAtom_New[j]]
                l2 = atomSeq[protAtom_New[k + l]]
                l3 = atomSeq[drugAtom[i]]

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
                  outputFile1.write("{}_{}_{}_{}  ".format(drug_name, chain_name, drug_id, drugAtom[i]))
                  outputFile1.write(
                      "{}_{}_{}_{}  ".format(prot_res_New[j], chain_name, prot_res_No_New[j], protAtom_New[j]))
                  outputFile1.write("{}_{}_{}_{}  ".format(prot_res_New[k + l], chain_name, prot_res_No_New[k + l],
                                                          protAtom_New[k + l]))
                  outputFile1.write(" {:.2f}  {:.2f}  {:.2f} ".format(L1, L2, L3))
                  outputFile1.write(" {:.2f},{:.2f},{:.2f}   {:.2f},{:.2f},{:.2f}  ".format(xCordDrug[i], yCordDrug[i],
                                                                                            zCordDrug[i],
                                                                                            xCordProt_New[j],
                                                                                            yCordProt_New[j],
                                                                                            zCordProt_New[j]))
                  outputFile1.write(" {:.2f},{:.2f},{:.2f}  ".format(xCordProt_New[k + l], yCordProt_New[k + l],
                                                                    zCordProt_New[k + l]))
                  outputFile1.write("{:.2f}    {:.2f}    {:.2f}   {:.0f}\n".format(Theta, maxDist, b3, key3D))
                  Total_key+=1


            l += 1
            c -= 1
    if output_option in ['keys', 'both']:
      for value_ in keyDict3D:
          outputFile2.writelines([str(value_), '\t', str(keyDict3D[value_]), '\n'])

    if output_option in ['triplets', 'both']:
      outputFile1.close()

    if output_option in ['keys', 'both']:
      outputFile2.close()

def DrugTSR(data_dir, input_files, chain=None, drug_name=None, drug_id=None, output_option='both', output_subdir='lexicographic'):
    os.makedirs(os.path.join(data_dir, output_subdir), exist_ok=True)

    # Use defaultdict to store lists for multiple chains per protein
    chain_dict = defaultdict(list)
    drug_name_dict = defaultdict(list)
    drug_id_dict = defaultdict(list)

    # Handle single file input
    if isinstance(input_files, str):
        if input_files.endswith('.csv'):
            # Read CSV and populate dictionaries for chain, drug name, and drug ID
            print("Reading CSV file...")
            df = pd.read_csv(input_files)
            for _, row in df.iterrows():
                protein = row['protein'].upper()
                chain_dict[protein].append(row['chain'])
                drug_name_dict[protein].append(row['drug_name'])
                drug_id_dict[protein].append(str(row['drug_id']))
            input_files = df['protein'].str.upper().unique().tolist()

        else:
            input_files = [input_files]

    # Initialize chain, drug_name, and drug_id dictionaries if input arguments are provided
    if chain and drug_name and drug_id:
        if isinstance(chain, list) and isinstance(drug_name, list) and isinstance(drug_id, list):
            for f, c, d, i in zip(input_files, chain, drug_name, drug_id):
                chain_dict[f.upper()].append(c)
                drug_name_dict[f.upper()].append(d)
                drug_id_dict[f.upper()].append(i)
        elif isinstance(chain, str) and isinstance(drug_name, str) and isinstance(drug_id, str):
            for f in input_files:
                chain_dict[f.upper()].append(chain)
                drug_name_dict[f.upper()].append(drug_name)
                drug_id_dict[f.upper()].append(drug_id)

    Parallel(n_jobs=num_cores, verbose=50)(
        delayed(process_pdb_file)(
            data_dir,
            protein.upper(),
            chain,
            drug_name.upper(),
            drug_id,
            output_subdir,
            output_option
        )
        for protein in input_files
        for chain, drug_name, drug_id in zip(
            chain_dict.get(protein.upper(), []),
            drug_name_dict.get(protein.upper(), []),
            drug_id_dict.get(protein.upper(), [])
        )
    )

'''
Example usage:
data_dir = 'Dataset'
pdb_files = ['4CI2', '4ci1']
chain=['B', 'B']
drug_name=['LVY', 'EF2']
drug_id=['1429', '21429']
PDB_DL(pdb_files)
DrugTSR(data_dir, input_files=pdb_files, chain=chain, drug_name=drug_name, drug_id=drug_id, output_option='both')

csv_file = 'sample_details_six_cancer_drugs.csv'
PDB_DL(csv_file)
DrugTSR(data_dir, input_files=csv_file, output_option='both')'''