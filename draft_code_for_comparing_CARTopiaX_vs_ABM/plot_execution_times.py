import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV file
df = pd.read_csv('./timesABM.csv')

# Function to parse time strings (e.g., "517m13,098s" -> seconds)
def parse_time(time_str):
    parts = time_str.replace(',', '.').split('m')
    minutes = float(parts[0])
    seconds = float(parts[1].replace('s', ''))
    return minutes * 60 + seconds

# Prepare data for plotting
data_to_plot = []
positions = []
colors = []

dose_categories = ['0doses', '1dose', '2doses']
sources = ['CARTopiaX', 'Nature']

for i, dose in enumerate(dose_categories):
    for j, source in enumerate(sources):
        col_name = f"{dose}_{source}"
        times = df[col_name].apply(parse_time).values
        data_to_plot.append(times)
        positions.append(i)  # Same position for same dose
        colors.append('royalblue' if source == 'CARTopiaX' else 'crimson')

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Separate data by source for connecting lines
cartopia_means = []
nature_means = []
cartopia_positions = []
nature_positions = []

# Plot min, max, and mean for each dataset
for i, (data, pos, color) in enumerate(zip(data_to_plot, positions, colors)):
    min_val = np.min(data)
    max_val = np.max(data)
    mean_val = np.mean(data)
    
    # Store means for connecting lines
    if color == 'royalblue':
        cartopia_means.append(mean_val)
        cartopia_positions.append(pos)
    else:
        nature_means.append(mean_val)
        nature_positions.append(pos)
    
    # Vertical line from min to max
    ax.plot([pos, pos], [min_val, max_val], color=color, linewidth=1.5, alpha=0.7)
    
    # Min and max horizontal bars
    ax.plot([pos-0.02, pos+0.02], [min_val, min_val], color=color, linewidth=1.5)
    ax.plot([pos-0.02, pos+0.02], [max_val, max_val], color=color, linewidth=1.5)
    
    # Mean point
    ax.plot(pos, mean_val, 'o', color=color, markersize=8, 
            markeredgecolor='black', markeredgewidth=0.5, zorder=3)

# Connect mean points
ax.plot(cartopia_positions, cartopia_means, '-', color='royalblue', linewidth=2, alpha=0.6, zorder=2)
ax.plot(nature_positions, nature_means, '-', color='crimson', linewidth=2, alpha=0.6, zorder=2)

# Set x-axis labels
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['0 doses', '1 dose', '2 doses'])

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='royalblue', alpha=0.7, label='CARTopiaX'),
                   Patch(facecolor='crimson', alpha=0.7, label='Nature Paper')]
ax.legend(handles=legend_elements, loc='upper right')

# Labels and title
ax.set_xlabel('Number of Doses', fontsize=12)
ax.set_ylabel('Execution Time (hours)', fontsize=12)
ax.set_title('Execution Time Comparison: CARTopiaX vs Nature Paper', fontsize=14)

# Set y-axis to start at 0
ax.set_ylim(bottom=0)

# Convert y-axis to hours
y_ticks = ax.get_yticks()
ax.set_yticklabels([f'{tick/3600:.1f}' for tick in y_ticks])

ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('execution_times_violin.png', dpi=300, bbox_inches='tight')
plt.show()