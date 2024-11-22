import os
import math
import multiprocessing
from joblib import Parallel, delayed
import pandas as pd
import pickle

dTheta = 29
dLen = 35
numOfLabels = 20
num_cores = multiprocessing.cpu_count()

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

# SEPERATE HELIX and SHEET PARTS FOR EACH PROTEIN, SAVE IN A PICKLE FILE
def seperate_helix_sheet_chain(file, data_dir):
    helix_dict={}
    sheet_dict={}

    # read file for helix and sheet
    with open(f'{data_dir}/{file.upper()}.pdb', 'r') as f:
        for line in f:
            if line[0:6].strip()=='HELIX':
                serNum=line[7:10].strip()
                chain=line[19:20].strip()
                start=int(line[21:25].strip())
                stop=int(line[33:37].strip())
                h_range=list(range(start,stop+1))
                for pos in h_range:
                    pos1 = str(pos)+chain
                    if pos1 not in helix_dict:
                        helix_dict[pos1]='HELIX_'+str(serNum)

            if line[0:6].strip()=='SHEET':
                serNum=line[7:10].strip()
                chain=line[21:22].strip()
                start=int(line[22:26].strip())
                stop=int(line[33:37].strip())
                s_range=list(range(start,stop+1))
                for pos in s_range:
                    pos1 = str(pos)+chain
                    if pos1 not in sheet_dict:
                        sheet_dict[pos1]='SHEET_'+str(serNum)

    fh=open(f'{data_dir}/{file.upper()}_helix_dict.pkl',"wb")
    fs=open(f'{data_dir}/{file.upper()}_sheet_dict.pkl',"wb")
    pickle.dump(helix_dict,fh)
    pickle.dump(sheet_dict,fs)
    fh.close()
    fs.close()
    f.close()
    print("Helix and Sheet seperation done.")

# FIND TYPE ASSIGNMENT FOR EACH TRIPLETS (HELIX, SHEET, NONE)
def determine_type(triplets):
    code=0
    t0 = triplets[0]
    t1 = triplets[1]
    t2 = triplets[2]

    # 3 SAME
    if (t0==t1==t2):
        if t0.split("_")[0]=='HELIX': # 3a1 - 3 vertices from same helix
            code = '1_3a1'
        if t0.split("_")[0]=='SHEET': # 3b1 - 3 vertices from same sheet
            code = '7_3b1'
        if t0.split("_")[0] =='NONE': # 3c - all vertices from none, no helix or sheet
            code = '17_3c'

    # 3 DIFF h/s/n
    if(t0.split("_")[0]!=t1.split("_")[0]!=t2.split("_")[0]):
        code = '18_1a1b1c' # 1a1b1c - all different, I helix, 1 sheet, 1 none

    # 3 SAME TYPE, DIFF NUM
    if (t0.split("_")[0]==t1.split("_")[0]==t2.split("_")[0]):
        if(t0.split("_")[0]=='HELIX' and t0.split("_")[1]!=t1.split("_")[1]!=t2.split("_")[1]):
            code = '3_3a3' # 3a3 - 3 vertices from three different helices
        if(t0.split("_")[0]=='SHEET' and t0.split("_")[1]!=t1.split("_")[1]!=t2.split("_")[1]):
            code = '9_3b3' # 3a3 - 3 vertices from three different sheets

    # 2 SAME, 1 DIFF.
    if(t0==t1 and t1!=t2) or (t1==t2 and t2!=t0) or (t0==t2 and t1!=t2):
        if(t0==t1 and t1!=t2):
            one=t0
            two=t1
            three=t2
        if(t1==t2 and t2!=t0):
            one=t1
            two=t2
            three=t0
        if(t0==t2 and t1!=t2):
            one=t0
            two=t2
            three=t1
        if one.split("_")[0]=='HELIX' and three.split("_")[0]=='HELIX': #3a2 - 2 vertices from same helix, 1 vertex from another helix
            code = '2_3a2'
        if one.split("_")[0]=='HELIX' and three.split("_")[0]=='SHEET': #2a11b 	2 vertices from same helix, 1 vertex from sheet
            code = '14_2a11b'
        if one.split("_")[0]=='HELIX' and three.split("_")[0]=='NONE': #2a11c 	2 vertices from same helix, 1 vertex from none
            code = '4_2a11c'

        if one.split("_")[0]=='SHEET' and three.split("_")[0]=='SHEET': #3a2 - 2 vertices from same sheet, 1 vertex from another sheet
            code = '8_3b2'
        if one.split("_")[0]=='SHEET' and three.split("_")[0]=='HELIX': #2b11a 	2 vertices from same sheet, 1 vertex from helix
            code = '16_2b11a'
        if one.split("_")[0]=='SHEET' and three.split("_")[0]=='NONE': # 2b11c 	2 vertices from same sheet, 1 vertex from none
            code = '10_2b11c'

        if one.split("_")[0]=='NONE' and three.split("_")[0]=='HELIX': #1a2c 	1 vertex from helix, two vertices from none
            code = '6_1a2c'
        if one.split("_")[0]=='NONE' and three.split("_")[0]=='SHEET': #1b2c 	1 vertex from sheet, two vertices from none
            code = '12_1b2c'

    # 2 SAMEDIFF, 1 DIFF.
    if((t0.split("_")[0]==t1.split("_")[0] and t1.split("_")[0]!=t2.split("_")[0]) or \
       (t1.split("_")[0]==t2.split("_")[0] and t2.split("_")[0]!=t0.split("_")[0]) or \
       (t2.split("_")[0]==t0.split("_")[0] and t0.split("_")[0]!=t1.split("_")[0])):

        if (t0.split("_")[0]==t1.split("_")[0] and t1.split("_")[0]!=t2.split("_")[0]):
            one_ = t0
            two_ = t1
            three_ = t2
        if (t1.split("_")[0]==t2.split("_")[0] and t2.split("_")[0]!=t0.split("_")[0]):
            one_ = t1
            two_ = t2
            three_ = t0
        if (t2.split("_")[0]==t0.split("_")[0] and t0.split("_")[0]!=t1.split("_")[0]):
            one_ = t2
            two_ = t0
            three_ = t1

        if(one_.split("_")[0]=='HELIX' and one_.split("_")[1]!=two_.split("_")[1] and three_.split("_")[0]=='SHEET'):
            code = '13_2a21b'  # 2a21b 	2 vertices from different helices, 1 vertex from sheet
        if(one_.split("_")[0]=='HELIX' and one_.split("_")[1]!=two_.split("_")[1] and three_.split("_")[0]=='NONE'):
            code = '5_2a21c'  # 2a21c 	2 vertices from two different helices, 1 vertex from none

        if(one_.split("_")[0]=='SHEET' and one_.split("_")[1]!=two_.split("_")[1] and three_.split("_")[0]=='HELIX'):
            code = '15_2b21a'  # 2b21a 	2 vertices from different sheets, 1 vertex from helix
        if(one_.split("_")[0]=='SHEET' and one_.split("_")[1]!=two_.split("_")[1] and three_.split("_")[0]=='NONE'):
            code = '11_2b21c' # 2b21c 	2 verticed from two different sheets, 1 vertex from none

    return code

#for fileName in files:
def generate_keys_and_triplets_ss_info(data_dir, file_name, chain, output_subdir, output_option='both',
                               aa_grouping=False, mirror_image=False, size_filter=10000):

    seperate_helix_sheet_chain(file_name, data_dir)
    if aa_grouping:
        aminoAcidDict = aminoAcidGroup
        print("Using amino acid group labels.")
    else:
        aminoAcidDict = aminoAcidLabel
        print("Using standard amino acid labels.")
    filesDict3D={}
    filesDict1D={}
    types = [0] * 18  # array of len 18 to save 18 type combination of helix/sheet for each key, each element of array is freq of key i of type j
    if not os.path.exists(f'{data_dir}/{file_name.upper()}.pdb'):
        pass
    inFile=open(f'{data_dir}/{file_name.upper()}.pdb','r')

    fileKey3D = open(f'{data_dir}/{output_subdir}/{file_name.upper()}.3Dkeys_29_35_SSE', "w") if output_option in ['both', 'keys'] else None
    fileTriplets = open(f'{data_dir}/{output_subdir}/{file_name.upper()}.triplets_29_35_SSE', "w") if output_option in ['both', 'triplets'] else None

    # Write Header !
    fileKey3D.writelines("key3D\t1_3a1\t2_3a2\t3_3a3\t4_2a11c\t5_2a21c\t6_1a2c \t7_3b1\t8_3b2\t9_3b3\t10_2b11c\t11_2b21c\t12_1b2c\t13_2a21b\t14_2a11b\t15_2b21a\t16_2b11a\t17_3c\t18_1a1b1c\n")
    fileTriplets.writelines("key3D\taa0\tpos0\taa1\tpos1\taa2\tpos2\tclassT1\tTheta\tclassL1\tmaxDist\tx0\ty0\tz0\tx1\ty1\tz1\tx2\ty2\tz2\tTheta1\theight\ttype\n")

    global xCord, yCord, zCord
    aminoAcidName={}
    xCord={}
    yCord={}
    zCord={}
    seq_number={}
    counter=0
    helix_dict=pickle.load(open(f'{data_dir}/{file_name.upper()}_helix_dict.pkl',"rb"))
    sheet_dict=pickle.load(open(f'{data_dir}/{file_name.upper()}_sheet_dict.pkl',"rb"))

    for i in inFile:
        if ((i[0:6].rstrip()=="ENDMDL") or (i[0:6].rstrip()=='TER' and i[21].rstrip()==chain)):
            break
        if (i[0:6].rstrip()=="MODEL" and int(i[10:14].rstrip())>1):
            break

        if(i[0:4].rstrip())=="ATOM"and(i[13:15].rstrip())=="CA"and(i[16]=='A' or i[16]==' ') and i[21:22].strip()==chain and i[17:20].strip()!= "UNK" :
            #print (i)

            if i[22:27].strip().isdigit() == False: # check if seq number is '123AB' or '123'
                if not ((i[22:27].strip().startswith("-"))==True and (i[22:27].strip()[1:].isdigit())==True):
                    continue

            aminoAcidName[counter]=int(aminoAcidDict[i[17:20].strip()])
            xCord[counter]=(float(i[30:38].strip()))
            yCord[counter]=(float(i[38:46].strip()))
            zCord[counter]=(float(i[46:54].strip()))
            seq_number[counter]=str(i[22:27].strip())
            counter+=1

    protLen=len(yCord)
    print(protLen)
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

                    ##getting the positions of AminoAcids in sequence
                    position0 = str(list(seq_number.values())[indexOf0])
                    position1 = str(list(seq_number.values())[indexOf1])
                    position2 = str(list(seq_number.values())[indexOf2])

                    aacd0 = list(aminoAcidDict.keys())[list(aminoAcidDict.values()).index(aminoAcidName[indexOf0])]
                    aacd1 = list(aminoAcidDict.keys())[list(aminoAcidDict.values()).index(aminoAcidName[indexOf1])]
                    aacd2 = list(aminoAcidDict.keys())[list(aminoAcidDict.values()).index(aminoAcidName[indexOf2])]

                    # GET HELIX or SHEET info
                    seqList = [str(position0)+chain, str(position1)+chain, str(position2)+chain]
                    typeList = []
                    for s in seqList:
                        if s in helix_dict:
                            typeList.append(helix_dict[s])
                        elif s in sheet_dict:
                            typeList.append(sheet_dict[s])
                        else:
                            typeList.append('NONE_0')
                    type_ = determine_type(typeList) # func call determine helix/sheet/none combination type for each triangle

                    x0 = str(xCord.get(indexOf0))
                    y0 = str(yCord.get(indexOf0))
                    z0 = str(zCord.get(indexOf0))

                    x1 = str(xCord.get(indexOf1))
                    y1 = str(yCord.get(indexOf1))
                    z1 = str(zCord.get(indexOf1))

                    x2 = str(xCord.get(indexOf2))
                    y2 = str(yCord.get(indexOf2))
                    z2 = str(zCord.get(indexOf2))

                    # GENERATE 3D key
                    key3D = dLen*dTheta*(numOfLabels**2)*(aminoAcidName[indexOf0]-1)+\
                            dLen*dTheta*(numOfLabels)*(aminoAcidName[indexOf1]-1)+\
                            dLen*dTheta*(aminoAcidName[indexOf2]-1)+\
                            dTheta*(classL1-1)+\
                            (classT1-1)

                    if mirror_image and Theta1 > 90:
                            key3D = (-1) * key3D

                    ## 3D ke-freq file generation with type of helix/sheet
                    if key3D in filesDict3D:
                        index = int(type_.split("_")[0])-1
                        types[index] += 1
                        filesDict3D[key3D] = types
                    else:
                        types = [0] * 18
                        index = int(type_.split("_")[0])-1
                        types[index] += 1
                        filesDict3D[key3D]=types

                    # WRITE LINE TO triplet file
                    if output_option in ['both', 'triplets'] and fileTriplets:
                        line = (str(key3D)+"\t"+\
                                str(aacd0)+"\t"+str(position0)+"\t"+str(aacd1)+"\t"+str(position1)+"\t"+str(aacd2)+"\t"+str(position2)+"\t"+\
                                str(classT1)+"\t"+str(Theta)+"\t"+str(classL1)+"\t"+str(maxDist)+"\t"+\
                                x0+"\t"+y0+"\t"+z0+"\t"+x1+"\t"+y1+"\t"+z1+"\t"+x2+"\t"+y2+"\t"+z2+"\t"+\
                                str(Theta1)+"\t"+str(s3)+"\t"+type_+"\n")
                        fileTriplets.writelines(line)

    ## Write lines in key-freq file with helix sheet type
    if output_option in ['both', 'keys'] and fileKey3D:
        for value_ in filesDict3D:
            types_ = ("\t").join([str(x) for x in filesDict3D[value_]])
            fileKey3D.writelines([str(value_),'\t', types_,'\n'])

    fileKey3D.close()
    fileTriplets.close()
    os.remove(f'{data_dir}/{file_name.upper()}_helix_dict.pkl')
    os.remove(f'{data_dir}/{file_name.upper()}_sheet_dict.pkl')

# Main function to handle input and output
def SSETSR(data_dir, input_files, chain=None, output_option='both', output_subdir='lexicographic', aa_grouping=False, mirror_image=False, size_filter=10000):
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
        delayed(generate_keys_and_triplets_ss_info)(data_dir, file_name.upper(), chain_dict.get(file_name.upper(), chain), output_subdir, output_option, aa_grouping, mirror_image, size_filter)
        for file_name in input_files
    )

'''
Example Usage:

data_dir = "Dataset"
input_files = ["1GTA", "1gtb", "1lbe"]
chain = ["A", "A", "A"]
output_option = "both"
PDB_DL(input_files)
SSETSR(data_dir, input_files, chain=chain, output_option=output_option, aa_grouping=False, mirror_image=False, size_filter=10000)

OR

data_dir = "Dataset"
csv_file = "sample_details.csv"
PDB_DL(csv_file) 
SSETSR(data_dir, csv_file, output_option="keys")
'''