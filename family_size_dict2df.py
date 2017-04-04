#!/usr/bin/env python3

import os
import collections
import pickle
import numpy as np
import pandas as pd

outfile_path = '/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/projects/CAPP/CellLineDilutions_v2/downsample_x10_x8dilutions/family_size_df/'

file_path = '/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/projects/CAPP/CellLineDilutions_v2/downsample_x10_x8dilutions/family_size'
os.chdir(file_path)
file_list = os.listdir(file_path)

for file in file_list:
    tag_fam_file = open(file, 'rb')
    tag_dict = pickle.load(tag_fam_file)

    # tags_per_fam_size = collections.Counter([i for i in tag_dict.values()])
    # lst_fam_per_read = list(tags_per_fam_size.items())  # convert to list

    # === Convert dictionary to dataframe ===

    tag_df = pd.DataFrame(list(tag_dict.items()), columns=['tag_ID', 'family_size'])
    tag_df_summary = tag_df.join(tag_df['tag_ID'].str.split('_', expand =True))
    tag_df_summary.columns = ['tag_ID', 'family_size','barcode', 'R1chr', 'R1start', 'R2chr', 'R2start', 'R1cigar', 'R2cigar', 'strand', 'orientation', 'RG']

    tag_df_summary.to_csv(outfile_path + file, index=None, sep='\t', mode='a')



