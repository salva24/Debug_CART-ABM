# -----------------------------------------------------------------------------
# Copyright (C) 2025 Salvador de la Torre Gonzalez
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Paths to the data files
csv_data_dir = './final_data_mine/'
num_csv_files_mine=1
dat_data_dir = './final_data_theirs/'
days = 30.1

# Read all CSV files from your model
csv_dataframes = []
for i in range(1, num_csv_files_mine + 1):  # Assuming you have final_data1.csv to final_data10.csv
    csv_file_path = os.path.join(csv_data_dir, f'final_data{i}.csv')
    if os.path.exists(csv_file_path):
        df_csv = pd.read_csv(csv_file_path)
        csv_dataframes.append(df_csv)
        print(f"Loaded {csv_file_path}")
    else:
        print(f"Warning: {csv_file_path} not found")

print(f"Successfully loaded {len(csv_dataframes)} CSV files")

# Updated column names for the .dat file
column_names = ['#tiempo', 'volumen', 'volumen2', 'radio', 'celulas tumorales', 
                'dead_cancer_cells', 'todas las celulas', 'cart_alive']

# Read all 10 .dat files
dat_dataframes = []
for i in range(1, 11):
    dat_file_path = os.path.join(dat_data_dir, f'DatosFinales{i}.dat')
    if os.path.exists(dat_file_path):
        df_dat = pd.read_csv(dat_file_path, sep=' ', names=column_names, skiprows=1)
        df_dat['#tiempo'] = pd.to_numeric(df_dat['#tiempo'])
        df_dat['total_days'] = df_dat['#tiempo'] / 1440  # 1440 minutes in a day
        dat_dataframes.append(df_dat)
        print(f"Loaded {dat_file_path}")
    else:
        print(f"Warning: {dat_file_path} not found")

print(f"Successfully loaded {len(dat_dataframes)} DAT files")

# Process your CSV data to calculate mean and std
csv_time_points = None
csv_mean_cells = None
csv_std_cells = None

if csv_dataframes:
    # Filter data by days and find common time points
    filtered_csv_dfs = [df[df['total_days'] <= days] for df in csv_dataframes]
    
    # Use the first file as reference for time points
    csv_time_points = filtered_csv_dfs[0]['total_days'].values
    
    # Collect cell counts for each time point across all runs
    all_csv_cell_counts = []
    for df_csv in filtered_csv_dfs:
        # Interpolate to common time points if needed
        cell_counts = np.interp(csv_time_points, df_csv['total_days'], df_csv['num_tumor_cells'])
        all_csv_cell_counts.append(cell_counts)
    
    # Convert to numpy array for easier computation
    all_csv_cell_counts = np.array(all_csv_cell_counts)
    
    # Calculate mean and standard deviation
    csv_mean_cells = np.mean(all_csv_cell_counts, axis=0)
    csv_std_cells = np.std(all_csv_cell_counts, axis=0)

# Process the Nature paper data to calculate mean and std
dat_time_points = None
dat_mean_cells = None
dat_std_cells = None

if dat_dataframes:
    # Find common time points (use the first file as reference)
    dat_time_points = dat_dataframes[0][dat_dataframes[0]['total_days'] <= days]['total_days'].values
    
    # Collect cell counts for each time point across all runs
    all_dat_cell_counts = []
    for df_dat in dat_dataframes:
        filtered_df = df_dat[df_dat['total_days'] <= days]
        # Interpolate to common time points if needed
        cell_counts = np.interp(dat_time_points, filtered_df['total_days'], filtered_df['celulas tumorales'])
        all_dat_cell_counts.append(cell_counts)
    
    # Convert to numpy array for easier computation
    all_dat_cell_counts = np.array(all_dat_cell_counts)
    
    # Calculate mean and standard deviation
    dat_mean_cells = np.mean(all_dat_cell_counts, axis=0)
    dat_std_cells = np.std(all_dat_cell_counts, axis=0)

# Create the plot
fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Number of Cells')

# Plot your CARTopiaX model data
if csv_dataframes:
    # Plot mean line
    ax1.plot(csv_time_points, csv_mean_cells, 
             color='red', linestyle='-', linewidth=2,
             label='CARTopiaX Model (Mean)')
    
    # Plot standard deviation band
    sd_multiplier = 3
    ax1.fill_between(csv_time_points, 
                     csv_mean_cells - sd_multiplier * csv_std_cells, 
                     csv_mean_cells + sd_multiplier * csv_std_cells,
                     color='red', alpha=0.3, 
                     label=f'CARTopiaX Model')

# Plot Nature paper data if available
if dat_dataframes:
    # Plot mean line
    ax1.plot(dat_time_points, dat_mean_cells, 
             color='blue', linestyle='-', linewidth=2,
             label='Nature Paper Model (Mean)')
    
    # Plot standard deviation band
    sd_multiplier = 3
    ax1.fill_between(dat_time_points, 
                     dat_mean_cells - sd_multiplier * dat_std_cells, 
                     dat_mean_cells + sd_multiplier * dat_std_cells,
                     color='blue', alpha=0.3, 
                     label=f'Nature Paper Model')

ax1.legend(loc='best')
ax1.grid(True, linestyle='--', alpha=0.7)

plt.title('Tumor Cell Count over Time: CARTopiaX vs Nature Paper Model')
plt.tight_layout()

# plt.savefig('./graphs_only_one_dosage_version8_9.png', dpi=300, bbox_inches='tight')
plt.show()
print("Plot saved")