import math
import os
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing

# Define constants
dTheta = 29
dLen = 35
numOfLabels = 20
num_cores = multiprocessing.cpu_count()

# Amino acid label mapping (hardcoded from the file)
aminoAcidLabel = {
    'PHE': 4, 'CYS': 5, 'GLN': 6, 'GLU': 7, 'LEU': 8, 'HIS': 9,
    'ILE': 10, 'GLY': 11, 'LYS': 12, 'MET': 13, 'ASP': 14, 'PRO': 15,
    'SER': 16, 'THR': 17, 'TRP': 18, 'TYR': 19, 'VAL': 20, 'ALA': 21,
    'ARG': 22, 'ASN': 23
}

aminoAcidGroup = {
    'ASP':	4, 'GLU':	4, 'ASN':	5, 'GLN':	5, 'GLY':	6, 'ALA':	7, 'VAL':	7,
    'LEU':	8, 'ILE':	8, 'PRO':	9, 'PHE':	10, 'TRP':	10, 'TYR':	11, 'SER':	12,
    'THR':	12, 'CYS':	13, 'MET':	14, 'HIS':	15, 'LYS':	16, 'ARG':	16
}

# Classify theta angles into predefined classes
def thetaClass_(Theta):
    theta_bins = [12.11, 17.32, 21.53, 25.21, 28.54, 31.64, 34.55, 37.34, 40.03, 42.64, 
                  45.17, 47.64, 50.05, 52.43, 54.77, 57.08, 59.38, 61.64, 63.87, 66.09, 
                  68.30, 70.5, 72.69, 79.2, 81.36, 83.51, 85.67, 87.80, 90.00]
    return next((i+1 for i, val in enumerate(theta_bins) if Theta < val), len(theta_bins))

# Classify distances into predefined classes
def dist12Class_(dist12):
    dist_bins = [3.83, 7.00, 9.00, 11.00, 14.00, 17.99, 21.25, 23.19, 24.8, 26.26, 
                 27.72, 28.9, 30.36, 31.62, 32.76, 33.84, 35.13, 36.26, 37.62, 38.73, 
                 40.12, 41.8, 43.41, 45.55, 47.46, 49.69, 52.65, 55.81, 60.2, 64.63, 
                 70.04, 76.15, 83.26, 132.45]
    return next((i+1 for i, val in enumerate(dist_bins) if dist12 < val), len(dist_bins))

# Calculate distance between two points
def calcDist(indexLabel1,indexLabel2):
    x1=xCord[indexLabel1]
    x2=xCord[indexLabel2]
    y1=yCord[indexLabel1]
    y2=yCord[indexLabel2]
    z1=zCord[indexLabel1]
    z2=zCord[indexLabel2]
    distance=(((x1-x2)**2+(y2-y1)**2+(z2-z1)**2)**0.5)
    return distance

def indexFind(index_of_2,i1,j1,k1):
    if index_of_2==i1:
        indexOf0=j1
        indexOf1=k1
    elif index_of_2==j1:
        indexOf0=i1
        indexOf1=k1
    elif index_of_2==k1:
        indexOf0=i1
        indexOf1=j1

    return indexOf0, indexOf1

# Parallel key and triplet file generation function
def generate_keys_and_triplets(data_dir, file_name, chain, output_subdir, output_option='both',
                               aa_grouping=False, mirror_image=False, size_filter=10000):

    if aa_grouping:
        aminoAcidDict = aminoAcidGroup
        print("Using amino acid group labels.")
    else:
        aminoAcidDict = aminoAcidLabel
        print("Using standard amino acid labels.")

    inFile = open(f'{data_dir}{file_name.upper()}.pdb', 'r')
    outFile2 = open(f'{data_dir}{output_subdir}/{file_name.upper()}.3Dkeys_29_35', "w") if output_option in ['both', 'keys'] else None
    fileTriplets = open(f'{data_dir}{output_subdir}/{file_name.upper()}.triplets_29_35', "w") if output_option in ['both', 'triplets'] else None

    global xCord, yCord, zCord
    aminoAcidName={}
    xCord={}
    yCord={}
    zCord={}
    seq_number={}
    counter=0

    # Read through the PDB file and extract necessary details
    for line in inFile:
        if line.startswith("ENDMDL") or (line.startswith("TER") and line[21].strip() == chain):
            break
        if line.startswith("MODEL") and int(line[10:14].strip()) > 1:
            break
        if line.startswith("ATOM") and line[13:15].strip() == "CA" and (line[16] == 'A' or line[16] == ' ') and line[21:22].strip() == chain and line[17:20] != "UNK":
            aminoAcidName[counter] = aminoAcidDict[line[17:20].strip()]
            xCord[counter] = float(line[30:38])
            yCord[counter] = float(line[38:46])
            zCord[counter] = float(line[46:54])
            seq_number[counter] = line[22:27].strip()
            counter += 1

    protLen = len(yCord)
    filesDict = {}

    # Process the triplets
    initialLabel=[]
    sortedLabel=[]
    sortedIndex=[]
    for m in range(0,3):
        initialLabel.append(0)
        sortedLabel.append(0)
        sortedIndex.append(0)
    for i in range(0,protLen-2):
        for j in range(i+1,protLen-1):
            for k in range(j+1, protLen):
                global i1,j1,k1
                i1=i
                j1=j
                k1=k
                keepLabelIndex={}
                keepLabelIndex[aminoAcidName[i]]=i
                keepLabelIndex[aminoAcidName[j]]=j
                keepLabelIndex[aminoAcidName[k]]=k
                initialLabel[0]=aminoAcidName[i]
                initialLabel[1]=aminoAcidName[j]
                initialLabel[2]=aminoAcidName[k]
                sortedLabel=list(initialLabel)
                sortedLabel.sort(reverse=True)
                if (sortedLabel[0]==sortedLabel[1])and(sortedLabel[1]==sortedLabel[2]):
                    dist1_2Temp=calcDist(i,j)
                    dist1_3Temp=calcDist(i,k)
                    dist2_3Temp=calcDist(j,k)
                    if dist1_2Temp>=(max(dist1_2Temp,dist1_3Temp,dist2_3Temp)):
                        indexOf0=i
                        indexOf1=j
                        indexOf2=k
                    elif dist1_3Temp>=(max(dist1_2Temp,dist1_3Temp,dist2_3Temp)):
                        indexOf0=i
                        indexOf1=k
                        indexOf2=j
                    else:
                        indexOf0=j
                        indexOf1=k
                        indexOf2=i
                elif(aminoAcidName[i]!=aminoAcidName[j])and(aminoAcidName[i]!=aminoAcidName[k])and(aminoAcidName[j]!=aminoAcidName[k]):
                    for index_ in range(0,3):
                        sortedIndex[index_]=keepLabelIndex[sortedLabel[index_]]
                    indexOf0=sortedIndex[0]
                    indexOf1=sortedIndex[1]
                    indexOf2=sortedIndex[2]
                elif(sortedLabel[0]==sortedLabel[1])and(sortedLabel[1]!=sortedLabel[2]):
                    indexOf2=keepLabelIndex[sortedLabel[2]]
                    indices=indexFind(indexOf2,i,j,k)
                    a=indexOf2
                    b=indices[0]
                    c=indices[1]
                    dist1_3Temp=calcDist(b,a)
                    dist2_3Temp=calcDist(c,a)
                    if dist1_3Temp>=dist2_3Temp:
                        indexOf0=indices[0]
                        indexOf1=indices[1]	
                    else:
                        indexOf0=indices[1]
                        indexOf1=indices[0]
                elif(sortedLabel[0]!=sortedLabel[1])and(sortedLabel[1]==sortedLabel[2]):
                    indexOf0=keepLabelIndex[sortedLabel[0]]
                    indices=indexFind(indexOf0,i,j,k)
                    if calcDist(indexOf0,indices[0])>= calcDist(indexOf0,indices[1]):
                        indexOf1=indices[0]
                        indexOf2=indices[1]	
                    else:
                        indexOf2=indices[0]
                        indexOf1=indices[1]
                dist01=calcDist(indexOf0,indexOf1)
                s2=dist01/2
                dist02=calcDist(indexOf0,indexOf2)
                s1=dist02
                dist12=dist01
                dist03=calcDist(indexOf1,indexOf2)
                maxDist=max(dist01,dist02,dist03)

                if maxDist<size_filter:
                    s3=(((xCord[indexOf0]+xCord[indexOf1])/2-xCord[indexOf2])**2+((yCord[indexOf0]+yCord[indexOf1])/2-yCord[indexOf2])**2+((zCord[indexOf0]+zCord[indexOf1])/2-zCord[indexOf2])**2)**0.5
                    Theta1=180*(math.acos((s1**2-s2**2-s3**2)/(2*s2*s3)))/3.14
                    if Theta1<=90:
                        Theta=Theta1
                    else:
                        Theta=abs(180-Theta1)
                    classT1=thetaClass_(Theta)
                    classL1=dist12Class_(maxDist)

                    position0 = str(list(seq_number.values())[indexOf0])
                    position1 = str(list(seq_number.values())[indexOf1])
                    position2 = str(list(seq_number.values())[indexOf2])

                    aacd0 = list(aminoAcidDict.keys())[list(aminoAcidDict.values()).index(aminoAcidName[indexOf0])]
                    aacd1 = list(aminoAcidDict.keys())[list(aminoAcidDict.values()).index(aminoAcidName[indexOf1])]
                    aacd2 = list(aminoAcidDict.keys())[list(aminoAcidDict.values()).index(aminoAcidName[indexOf2])]

                    x0 = str(xCord.get(indexOf0))
                    y0 = str(yCord.get(indexOf0))
                    z0 = str(zCord.get(indexOf0))

                    x1 = str(xCord.get(indexOf1))
                    y1 = str(yCord.get(indexOf1))
                    z1 = str(zCord.get(indexOf1))

                    x2 = str(xCord.get(indexOf2))
                    y2 = str(yCord.get(indexOf2))
                    z2 = str(zCord.get(indexOf2))

                    key_2=dLen*dTheta*(numOfLabels**2)*(aminoAcidName[indexOf0]-1)+dLen*dTheta*(numOfLabels)*(aminoAcidName[indexOf1]-1)+dLen*dTheta*(aminoAcidName[indexOf2]-1)+dTheta*(classL1-1)+(classT1-1)
                    
                    if mirror_image and Theta1 > 90:
                        key_2 = (-1) * key_2

                    if key_2 in filesDict:
                        filesDict[key_2]+=1
                    else:
                        filesDict[key_2]=1

                    # Write to triplets file if needed
                    if output_option in ['both', 'triplets'] and fileTriplets:
                        line = f"{key_2}\t{aacd0}\t{position0}\t{aacd1}\t{position1}\t{aacd2}\t{position2}\t{classT1}\t{Theta}\t{classL1}\t{maxDist}\t{xCord[indexOf0]}\t{yCord[indexOf0]}\t{zCord[indexOf0]}\t{xCord[indexOf1]}\t{yCord[indexOf1]}\t{zCord[indexOf1]}\t{xCord[indexOf2]}\t{yCord[indexOf2]}\t{zCord[indexOf2]}\n"
                        fileTriplets.writelines(line)

    # Write to keys file if needed
    if output_option in ['both', 'keys'] and outFile2:
        for value_ in filesDict:
            outFile2.writelines([str(value_), '\t', str(filesDict[value_]), '\n'])

    # Close the files
    if outFile2:
        outFile2.close()
    if fileTriplets:
        fileTriplets.close()
    inFile.close()

# Main function to handle input and output
def TSR(data_dir, input_files, chain=None, output_option='both', output_subdir='lexiographic', aa_grouping=False, mirror_image=False, size_filter=10000):
    os.makedirs(os.path.join(data_dir, output_subdir), exist_ok=True)
    chain_dict = {}
    # Handle single file input
    if isinstance(input_files, str):
        if input_files.endswith('.csv'):
            # Read CSV
            df = pd.read_csv(input_files)
            chain_dict = dict(zip(df['protein'].str.upper(), df['chain']))
            input_files = df['protein'].str.upper().tolist()
        else:
            input_files = [input_files]

    if chain:
        if isinstance(chain, list):
            chain_dict = {f.upper(): c for f, c in zip(input_files, chain)}
        elif isinstance(chain, str):
            chain_dict = {f.upper(): chain for f in input_files}

    # Parallel processing
    Parallel(n_jobs=num_cores, verbose=50)(
        delayed(generate_keys_and_triplets)(data_dir, file_name.upper(), chain_dict.get(file_name.upper(), chain), output_subdir, output_option, aa_grouping, mirror_image, size_filter)
        for file_name in input_files
    )