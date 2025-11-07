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
num_csv_files_mine = 3
paper_data_dir = './final_data_theirs/'
num_csv_files_paper = 10
days = 30.1
sd_multiplier = 3

# Read all CSV files from your model
csv_dataframes = []
for i in range(1, num_csv_files_mine + 1):
    csv_file_path = os.path.join(csv_data_dir, f'final_data{i}.csv')
    if os.path.exists(csv_file_path):
        df_csv = pd.read_csv(csv_file_path)
        csv_dataframes.append(df_csv)
        print(f"Loaded {csv_file_path}")
    else:
        print(f"Warning: {csv_file_path} not found")

print(f"Successfully loaded {len(csv_dataframes)} CSV files from your model")

# Read all CSV files from the Nature paper model
paper_dataframes = []
for i in range(1, num_csv_files_paper + 1):
    paper_file_path = os.path.join(paper_data_dir, f'datos_finales{i}.csv')
    if os.path.exists(paper_file_path):
        df_paper = pd.read_csv(paper_file_path)
        paper_dataframes.append(df_paper)
        print(f"Loaded {paper_file_path}")
    else:
        print(f"Warning: {paper_file_path} not found")

print(f"Successfully loaded {len(paper_dataframes)} CSV files from the Nature paper model")

# Process your CSV data to calculate mean and std
csv_time_points = None
csv_mean_cells = None
csv_std_cells = None

if csv_dataframes:
    filtered_csv_dfs = [df[df['total_days'] <= days] for df in csv_dataframes]
    csv_time_points = filtered_csv_dfs[0]['total_days'].values
    all_csv_cell_counts = []
    for df_csv in filtered_csv_dfs:
        cell_counts = np.interp(csv_time_points, df_csv['total_days'], df_csv['num_tumor_cells'])
        all_csv_cell_counts.append(cell_counts)
    all_csv_cell_counts = np.array(all_csv_cell_counts)
    csv_mean_cells = np.mean(all_csv_cell_counts, axis=0)
    csv_std_cells = np.std(all_csv_cell_counts, axis=0)

# Process the Nature paper data to calculate mean and std
paper_time_points = None
paper_mean_cells = None
paper_std_cells = None

if paper_dataframes:
    filtered_paper_dfs = [df[df['total_days'] <= days] for df in paper_dataframes]
    paper_time_points = filtered_paper_dfs[0]['total_days'].values
    all_paper_cell_counts = []
    for df_paper in filtered_paper_dfs:
        cell_counts = np.interp(paper_time_points, df_paper['total_days'], df_paper['num_tumor_cells'])
        all_paper_cell_counts.append(cell_counts)
    all_paper_cell_counts = np.array(all_paper_cell_counts)
    paper_mean_cells = np.mean(all_paper_cell_counts, axis=0)
    paper_std_cells = np.std(all_paper_cell_counts, axis=0)

# Create the plot
fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Number of Cells')

# Plot your model data (tumor cells)
if csv_dataframes:
    ax1.plot(csv_time_points, csv_mean_cells, 
             color='red', linestyle='-', linewidth=2,
             label='CARTopiaX Tumor Cells (Mean)')
    ax1.fill_between(csv_time_points, 
                     csv_mean_cells - sd_multiplier * csv_std_cells, 
                     csv_mean_cells + sd_multiplier * csv_std_cells,
                     color='red', alpha=0.2)
    # Plot alive CART cells
    all_cart_counts = []
    for df_csv in filtered_csv_dfs:
        cart_counts = np.interp(csv_time_points, df_csv['total_days'], df_csv['num_alive_cart'])
        all_cart_counts.append(cart_counts)
    all_cart_counts = np.array(all_cart_counts)
    mean_cart = np.mean(all_cart_counts, axis=0)
    std_cart = np.std(all_cart_counts, axis=0)
    ax1.plot(csv_time_points, mean_cart, color='green', linestyle='-', linewidth=2, label='CARTopiaX Alive CART Cells (Mean)')
    ax1.fill_between(csv_time_points, mean_cart - sd_multiplier * std_cart, mean_cart + sd_multiplier * std_cart,
                     color='green', alpha=0.2)

# Plot Nature paper data (tumor cells)
if paper_dataframes:
    ax1.plot(paper_time_points, paper_mean_cells, 
             color='red', linestyle='--', linewidth=2,
             label='Nature Paper Tumor Cells (Mean)')
    ax1.fill_between(paper_time_points, 
                     paper_mean_cells - sd_multiplier * paper_std_cells, 
                     paper_mean_cells + sd_multiplier * paper_std_cells,
                     color='red', alpha=0.2)
    # Plot alive CART cells
    all_cart_counts_paper = []
    for df_paper in filtered_paper_dfs:
        cart_counts_paper = np.interp(paper_time_points, df_paper['total_days'], df_paper['num_alive_cart'])
        all_cart_counts_paper.append(cart_counts_paper)
    all_cart_counts_paper = np.array(all_cart_counts_paper)
    mean_cart_paper = np.mean(all_cart_counts_paper, axis=0)
    std_cart_paper = np.std(all_cart_counts_paper, axis=0)
    ax1.plot(paper_time_points, mean_cart_paper, color='green', linestyle='--', linewidth=2, label='Nature Paper Alive CART Cells (Mean)')
    ax1.fill_between(paper_time_points, mean_cart_paper - sd_multiplier * std_cart_paper, mean_cart_paper + sd_multiplier * std_cart_paper,
                     color='green', alpha=0.2)

ax1.legend(loc='best')
ax1.grid(True, linestyle='--', alpha=0.7)

plt.title('Tumor & Alive CART Cells Over Time: CARTopiaX vs Nature Paper Model')
plt.tight_layout()

# plt.savefig('./graphs_free_initial_1_1_cart_2_dosages_num_cells.png', dpi=300, bbox_inches='tight')
plt.show()
# print("Plot saved")

cell_types = [
    ('tumor_cells_type1', 'Type 1'),
    ('tumor_cells_type2', 'Type 2'),
    ('tumor_cells_type3', 'Type 3'),
    ('tumor_cells_type4', 'Type 4'),
    ('tumor_cells_type5_dead', 'Type 5 (Dead)')
]

colors_mine = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']
colors_paper = ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6']

fig, ax = plt.subplots(figsize=(12, 7))

for idx, (col, label) in enumerate(cell_types):
    # Your model
    if csv_dataframes:
        all_type_percent = []
        for df_csv in filtered_csv_dfs:
            # Calculate percentage for each time point
            percent = np.interp(csv_time_points, df_csv['total_days'],
                                100 * df_csv[col] / df_csv['num_tumor_cells'])
            all_type_percent.append(percent)
        all_type_percent = np.array(all_type_percent)
        mean_percent = np.mean(all_type_percent, axis=0)
        std_percent = np.std(all_type_percent, axis=0)
        ax.plot(csv_time_points, mean_percent, color=colors_mine[idx], label=f'CARTopiaX {label}')
        ax.fill_between(csv_time_points, mean_percent - std_percent, mean_percent + std_percent,
                        color=colors_mine[idx], alpha=0.2)

    # Paper model
    if paper_dataframes:
        all_type_percent_paper = []
        for df_paper in filtered_paper_dfs:
            percent_paper = np.interp(paper_time_points, df_paper['total_days'],
                                      100 * df_paper[col] / df_paper['num_tumor_cells'])
            all_type_percent_paper.append(percent_paper)
        all_type_percent_paper = np.array(all_type_percent_paper)
        mean_percent_paper = np.mean(all_type_percent_paper, axis=0)
        std_percent_paper = np.std(all_type_percent_paper, axis=0)
        ax.plot(paper_time_points, mean_percent_paper, color=colors_paper[idx], linestyle='--', label=f'Paper {label}')
        ax.fill_between(paper_time_points, mean_percent_paper - std_percent_paper, mean_percent_paper + std_percent_paper,
                        color=colors_paper[idx], alpha=0.2)

ax.set_xlabel('Time (days)')
ax.set_ylabel('Percentage of Tumor Cells (%)')
ax.set_title('Percentage of Each Tumor Cell Type Over Time')
ax.legend(loc='upper left', fontsize=9)
ax.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
# plt.savefig('./graphs_free_initial_1_1_cart_2_dosages_type_of_cells_percentage.png', dpi=300, bbox_inches='tight')
plt.show()