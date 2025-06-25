import pandas as pd
import numpy as np
from semopy import Model, Optimizer, gather_statistics
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.patches import FancyArrowPatch
np.set_printoptions(suppress=True)

# ----------------------------
# 1. Data Preparation
# ----------------------------
# Load data (replace with your file path)
data = pd.read_csv("analyse_data.csv")

# Select target variables (example variables, adjust according to actual data)
bottom = False
if bottom:
    env_factors = ['NDVI_mean', 'WET_mean', 'NDBSI_mean', 'LST_mean', 'RSEI_mean', 'SI_mean']
else:
    env_factors = ['NDVI_std', 'WET_std', 'NDBSI_std', 'LST_std', 'RSEI_std', 'SI_std']
biodiversity = ['Shannon', 'Fisher_alpha', 'Simpson', 'Pielou_evenness']
df = data[env_factors + biodiversity].dropna()

# Standardize data (SEM recommends standardization)
df_std = (df - df.mean()) / df.std()

# Check data
print("First 5 rows of data:\n", df_std.head())
print("\nDescriptive statistics:\n", df_std.describe())

# ----------------------------
# 2. Build SEM Model
# ----------------------------
# Define model syntax (similar to lavaan syntax)
if bottom:
    model_spec = """
    # Structural model (retaining unidirectional paths)
    RSEI_mean ~ NDVI_mean + NDBSI_mean + WET_mean + LST_mean
    SI_mean ~ NDVI_mean + NDBSI_mean + WET_mean + LST_mean
    Shannon ~ RSEI_mean + SI_mean
    Fisher_alpha ~ RSEI_mean + SI_mean
    Simpson ~ RSEI_mean + SI_mean
    Pielou_evenness ~ RSEI_mean + SI_mean

    # Residual correlations between environmental variables (alternative to bidirectional relationships)
    NDVI_mean ~~ NDBSI_mean
    NDVI_mean ~~ WET_mean
    NDVI_mean ~~ LST_mean
    NDBSI_mean ~~ WET_mean
    NDBSI_mean ~~ LST_mean
    WET_mean ~~ LST_mean
    RSEI_mean ~~ SI_mean

    # Residual correlations between biodiversity indices
    Shannon ~~ Fisher_alpha
    Shannon ~~ Simpson
    Shannon ~~ Pielou_evenness
    Fisher_alpha ~~ Simpson
    Fisher_alpha ~~ Pielou_evenness
    Simpson ~~ Pielou_evenness
    """
else:
    model_spec = """
    # Structural model
    RSEI_std ~ NDVI_std + NDBSI_std + WET_std + LST_std
    SI_std ~ NDVI_std + NDBSI_std + WET_std + LST_std
    Shannon ~ RSEI_std + SI_std
    Fisher_alpha ~ RSEI_std + SI_std
    Simpson ~ RSEI_std + SI_std
    Pielou_evenness ~ RSEI_std + SI_std

    # Residual correlations between environmental variables (alternative to bidirectional relationships)
    NDVI_std ~~ NDBSI_std
    NDVI_std ~~ WET_std
    NDVI_std ~~ LST_std
    NDBSI_std ~~ WET_std
    NDBSI_std ~~ LST_std
    WET_std ~~ LST_std
    RSEI_std ~~ SI_std

    # Residual correlations between biodiversity indices
    Shannon ~~ Fisher_alpha
    Shannon ~~ Simpson
    Shannon ~~ Pielou_evenness
    Fisher_alpha ~~ Simpson
    Fisher_alpha ~~ Pielou_evenness
    Simpson ~~ Pielou_evenness
    """

# Create model
model = Model(model_spec)
model.fit(df_std)  # Bind data

# Parameter estimation (removed 'ML' parameter)
opt = Optimizer(model)
opt.optimize()  # Use default optimization method

# ----------------------------
# 3. Model Evaluation
# ----------------------------
# Get fit indices
stats = gather_statistics(model)
print("\nModel fit indices:")
print(f"Chi-square (χ²): {stats.chi2[0]:.3f}, df: {stats.dof}")
print(f"CFI: {stats.cfi:.3f}, RMSEA: {stats.rmsea:.3f}")

# Check parameter estimation results (standardized coefficients)
results = model.inspect(std_est=True)
print("\nParameter estimation results (standardized coefficients):\n", results)

import graphviz

# ----------------------------
# Configuration parameters (modify as needed)
# ----------------------------
output_format = "png"  # Options: pdf, svg, png, etc.
output_filename = "sem_graph"  # Output file name (without extension)

# ----------------------------
# Data preparation
# ----------------------------
path_data = results[results['op'] == '~'].copy()
cov_data = results[results['op'] == '~~'].copy() if '~~' in results['op'].values else None

# Normalize edge width (for visualizing weights)
path_data['scaled_width'] = np.interp(
    np.abs(path_data['Est. Std']), 
    (0, np.abs(path_data['Est. Std']).max()), 
    (1, 5)
)

# ----------------------------
# Create Graphviz directed graph
# ----------------------------
dot = graphviz.Digraph(
    engine='dot',  # Use dot layout engine (hierarchical layout)
    graph_attr={
        'dpi': '300',
        'rankdir': 'TD',     # Top-to-bottom layout
        'splines': 'true',   # Use curved connections
        'overlap': 'false',  # Prevent node overlap
        'nodesep': '0.5',    # Horizontal spacing between nodes
        'ranksep': '0.8',    # Vertical spacing between layers
        'fontname': 'Arial',
    },
    node_attr={
        'shape': 'box',
        'style': 'rounded,filled',
        'fontname': 'Arial',
        'fontsize': '20',
    }
)

# ----------------------------
# Dynamically determine node types (based on the `bottom` parameter)
# ----------------------------
if bottom:
    suffix = '_mean'
else:
    suffix = '_std'
special_nodes = ['RSEI', 'SI']  # Special node prefixes

# Classify nodes
all_nodes = set(path_data['lval']).union(set(path_data['rval']))
env_nodes = [n for n in all_nodes if n.endswith(suffix) and not any(n.startswith(s) for s in special_nodes)]
bio_nodes = [n for n in all_nodes if not n.endswith(suffix)]
special_nodes = [f"{s}{suffix}" for s in special_nodes]

# ----------------------------
# Add nodes (automatic color classification)
# ----------------------------
# Environmental factor nodes (left side)
for node in env_nodes:
    dot.node(
        node, 
        label=node.replace(suffix, ''),  # Remove suffix for display
        color='#76B7B2', 
        fillcolor='#76B7B280',  # Semi-transparent fill
        gradientangle='270'
    )

# Special nodes (RSEI/SI, centered)
for node in special_nodes:
    dot.node(
        node, 
        label=node.replace(suffix, ''), 
        color='#F28E2B', 
        fillcolor='#F28E2B80'
    )

# Biodiversity nodes (right side)
for node in bio_nodes:
    dot.node(
        node, 
        label=node, 
        color='#59A14F', 
        fillcolor='#59A14F80'
    )

# ----------------------------
# Add paths (with weights and direction)
# ----------------------------
for _, row in path_data.iterrows():
    color = '#4E79A7' if row['Est. Std'] >= 0 else '#E15759'
    dot.edge(
        row['rval'], 
        row['lval'],
        label=f" r={row['Est. Std']:.2f} p={row['p-value']:.2f}",  # Display standardized coefficients
        color=color,
        penwidth=str(row['scaled_width']),  # Line width reflects effect size
        arrowhead='normal',
        arrowsize='0.8',
        fontsize='20',
        fontcolor=color
    )

# ----------------------------
# Add residual correlations (dashed lines)
# ----------------------------
if cov_data is not None:
    for _, row in cov_data.iterrows():
        if row['lval'] != row['rval']:  # Exclude variance terms
            dot.edge(
                row['lval'], 
                row['rval'],
                label=f" r={row['Est. Std']:.2f} p={row['p-value']:.2f}",
                color='gray60',
                style='dashed',
                dir='none',  # No direction
                penwidth='1.2'
            )

# ----------------------------
# Render and save the graph
# ----------------------------
dot.render(
    filename=output_filename, 
    format=output_format,
    cleanup=True  # Automatically clean up temporary files
)

print(f"Graph saved as {output_filename}.{output_format}")