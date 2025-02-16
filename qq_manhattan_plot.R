import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Adjust the file path as necessary
file_path = 'Metal_test1.tbl'

# Read the data, using '\s+' as the separator to handle one or more spaces
df = pd.read_csv(file_path, sep='\s+')

# Parse the MarkerName to extract chromosome, position, reference allele, and alternative allele
df[['Chromosome', 'Position', 'Ref', 'Alt']] = df['MarkerName'].str.split(':', expand=True)
df['Chromosome'] = pd.to_numeric(df['Chromosome'], errors='coerce').astype(str)  # Keep as string for categorical plotting
df['Position'] = pd.to_numeric(df['Position'], errors='coerce')

# Calculate -log10(P-value)
df['minusLog10Pvalue'] = -np.log10(df['P-value'])

# Clean and prepare data
df.dropna(subset=['Chromosome', 'Position', 'minusLog10Pvalue'], inplace=True)

# Sort by chromosome and position for the Manhattan plot
df.sort_values(by=['Chromosome', 'Position'], inplace=True)

# Map chromosomes to numeric values for plotting
chromosome_map = {chrom: i for i, chrom in enumerate(df['Chromosome'].unique(), start=1)}
df['ChromNumeric'] = df['Chromosome'].map(chromosome_map)

# QQ Plot
plt.figure(figsize=(8, 6))
observed = np.sort(df['P-value'])
expected = np.arange(1, len(observed) + 1) / len(observed)
plt.plot(-np.log10(expected), -np.log10(observed), 'o')
plt.plot([0, max(-np.log10(expected))], [0, max(-np.log10(expected))], 'k--')
plt.xlabel('Expected -log10(p-value)')
plt.ylabel('Observed -log10(p-value)')
plt.title('QQ Plot')
plt.show()

# Manhattan Plot with Chromosomes on the x-axis
plt.figure(figsize=(12, 6))
colors = ['skyblue', 'salmon']  

# Significance threshold for the Manhattan plot (e.g., 5e-8)
significance_threshold = 5e-8
threshold_line = -np.log10(significance_threshold)

# Plot each chromosome's data points
for chrom, group in df.groupby('Chromosome'):
    plt.scatter(group['ChromNumeric'], group['minusLog10Pvalue'],
                color=colors[int(chromosome_map[chrom]) % len(colors)],
                alpha=0.6, edgecolor='k', linewidth=0.1, label=f'Chr {chrom}')

# Add a horizontal line for the significance threshold
plt.axhline(y=threshold_line, color='red', linestyle='--', label=f'Significance threshold (-log10({significance_threshold}))')

# Adjust the x-axis to show chromosome labels instead of numeric values
plt.xticks(ticks=list(chromosome_map.values()), labels=list(chromosome_map.keys()), rotation=45)

plt.xlabel('Chromosome')
plt.ylabel('-log10(P-value)')
plt.title('Manhattan Plot')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Chromosome")
plt.tight_layout()
plt.savefig('mnth.png')
plt.show()
