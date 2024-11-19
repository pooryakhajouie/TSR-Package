import pandas as pd
import numpy as np
import glob
import os
import time
import multiprocessing
from joblib import Parallel, delayed
import pickle
import math
from matplotlib import pyplot as plt
from PIL import Image

def KeyToImage(input_folder, extension, csv_location, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    num_cores = multiprocessing.cpu_count()

    # Read grid CSV file for mapping
    map_ = pd.read_csv(csv_location, header=0)
    map_['combined'] = ['_'.join(x) for x in np.sort(map_[['TripletName1', 'TripletName2', 'TripletName3']], axis=1)]
    map_ = map_.apply(lambda x: x.astype(str).str.upper())

    # Get the list of triplet files
    filelist = [x.split("/")[-1].split(".")[0] for x in glob.glob(os.path.join(input_folder, '*' + extension))]
    print("Processing files:", filelist)

    # Function to process each chunk in parallel
    def create_image_chunk(file, j, chunksize, n, triplet_data, name, angle, distance):
        tmpimage = np.zeros((35 * 35, 44 * 29))
        start_idx = j * chunksize
        stop_idx = min(n, (j + 1) * chunksize)

        for i in range(start_idx, stop_idx):
            triplets = name.iloc[i].values.tolist()
            triplets_ = ['HIS' if t in ['HIS', 'HIE', 'HID'] else t for t in triplets]

            if any(type(t) != str or len(t.strip()) < 3 for t in triplets_):
                continue

            triplets_ = '_'.join(sorted(map(str, triplets_))).strip()
            idx = map_[map_['combined'] == triplets_]
            if idx.empty:
                continue

            block_row = int(idx['RowIndex'].iloc[0])
            block_col = int(idx['ColIndex'].iloc[0])

            if math.isnan(distance.iloc[i, 0]) or math.isnan(angle.iloc[i, 0]):
                continue

            idx_i = (block_row - 1) * 35 + int(distance.iloc[i, 0]) - 1
            idx_j = (block_col - 1) * 29 + int(angle.iloc[i, 0]) - 1
            tmpimage[idx_i, idx_j] += 1.0

        pickle.dump(tmpimage, open(os.path.join(input_folder, f'{file}_{j}_tmpimg.p'), 'wb'))

    # Process each file
    for file in filelist:
        if os.path.exists(os.path.join(output_folder, f'{file}.image')):
            continue

        triplet_data = pd.read_csv(os.path.join(input_folder, f'{file}{extension}'), delimiter="\t", names=range(20))
        if triplet_data.empty:
            print(f"{file} is empty, skipping...")
            continue

        name = triplet_data[[1, 3, 5]]
        angle = triplet_data[[7]]
        distance = triplet_data[[9]]
        image = np.zeros((35 * 35, 44 * 29))
        n = len(triplet_data)

        chunksize = 100000
        noofchunks = math.ceil(n / chunksize)
        Parallel(n_jobs=num_cores)(delayed(create_image_chunk)(file, j, chunksize, n, triplet_data, name, angle, distance) for j in range(noofchunks))

        # Combine the chunks into the final image
        for j in range(noofchunks):
            tmpimage = pickle.load(open(os.path.join(input_folder, f'{file}_{j}_tmpimg.p'), 'rb'))
            image += tmpimage
            os.remove(os.path.join(input_folder, f'{file}_{j}_tmpimg.p'))

        np.savetxt(os.path.join(output_folder, f'{file}.image'), image, fmt='%1.2f', delimiter=' ')

        # Convert the .image file to a .png
        image = np.loadtxt(os.path.join(output_folder, f'{file}.image'), delimiter=' ')
        image = np.log(image + 1)
        plt.clf()
        plt.imshow(image, aspect='auto', interpolation='none')
        plt.axis('off')
        plt.savefig(os.path.join(output_folder, f'{file}.png'), bbox_inches='tight', pad_inches=0)
        print(f"Processed {file}")

    print("All files processed.")

'''
Sample usage:
input_folder = "sample_data/"  # Replace with the actual path to your input folder
extension = ".triplets_29_35"              # File extension for the triplet files
csv_location = "35by44grid.csv"   # Replace with the path to your 35by44grid.csv file
output_folder = "sample_data/"  # Replace with the desired output folder path

# Call the function to process the triplets and generate the .image and .png files
KeyToImage(input_folder, extension, csv_location, output_folder)

'''