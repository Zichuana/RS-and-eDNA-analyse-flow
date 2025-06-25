import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import warnings
import networkx as nx
from jinja2 import Template
import os
from tqdm import tqdm

warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'Times New Roman',
    'axes.unicode_minus': False,
    'figure.max_open_warning': 0,
    'font.size':26
})

# =============================================
# 1. 数据加载与预处理
# =============================================
def load_and_preprocess(species_path, env_path, 
                        min_samples=0.7, 
                        min_total_abundance=0):
    """
    Function for loading and preprocessing large-scale species data
    Includes species filtering functionality (based on occurrence frequency and total abundance)

    Parameters:
        species_path: Path to species abundance table
        env_path: Path to environmental factors table
        min_samples: Minimum sample occurrence proportion for species (default 10%)
        min_total_abundance: Minimum total abundance across all samples for species (default 10)

    Returns:
        species_filtered: Filtered species abundance table
        env_scaled: Standardized environmental factors table
    """
    print("\n[Data loading and preprocessing]")
    
    # 加载数据
    species = pd.read_csv(species_path, index_col=0)
    env = pd.read_csv(env_path, index_col=0)
    
    # 检查数据
    print(f"Original species data shape: {species.shape} (Sample × Species)")
    print(f"Environmental data shape: {env.shape} (Sample × Variable)")
    
    # 样本对齐
    common_samples = species.index.intersection(env.index)
    species = species.loc[common_samples]
    env = env.loc[common_samples]
    print(f"The number of samples after alignment: {len(common_samples)}")
    
    # Species Filtering
    #1. Filter by occurrence frequency (at least in X% of the samples)
    # if 0 < min_samples < 1:
    #     min_samples = int(len(species) * min_samples)
    # keep_species = (species > 0).sum(axis=0) >= min_samples

    #2. Filter by total abundance (where the total abundance in all samples is greater than the threshold)
    # keep_species &= species.sum(axis=0) >= min_total_abundance
    # species_filtered = species.loc[:, keep_species]
    # print(f" Number of species after filtration: {species_filteres.shape [1]}"
    # f"(Filter out the species {species.shape[1]- species_filteres.shape [1]})")
    
    # Standardization of environmental variables
    env_scaled = pd.DataFrame(
        StandardScaler().fit_transform(env),
        index=env.index,
        columns=env.columns
    )
    print(env_scaled)
    
    return species, env_scaled


import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import statsmodels.api as sm

def analyze_dataset(species, env, n_jobs=-1):
    """Support multi-species parallel analysis and incorporate hierarchical FDR correction"""
    print("\n[Statistical analysis]")
    if len(species) < 10:
        print(f"Warning: The sample size is too small(n={len(species)}), Caution is needed when interpreting the results!")
    
    # 1. Cluster and group environmental variables (based on Spearman correlation distance)
    env_corr = env.corr(method='spearman').values
    env_dist = 1 - np.abs(env_corr) 
    env_linkage = linkage(env_dist, method='average') 
    env_groups = fcluster(env_linkage, t=0.5, criterion='distance')  
    
    # Store the grouping information in the dictionary (optional: Print the grouping results)
    env_group_dict = {col: f"EnvGroup_{g}" for col, g in zip(env.columns, env_groups)}
    print("\n环境变量分组结果：")
    for group in set(env_groups):
        group_vars = [col for col, g in env_group_dict.items() if g == f"EnvGroup_{group}"]
        print(f"Group {group}: {group_vars}")
    
    all_results = []
    
    # 3. Species-by-species analysis
    for species_col in tqdm(species.columns, desc="Analyze the progress"):
        species_var = species[species_col]
        results = []
        
        # Analyze each environmental variable one by one
        for env_col in env.columns:
            # Rank correlation coefficient
            corr, pval = spearmanr(species_var, env[env_col])
            
            # Linear regression
            X = sm.add_constant(env[env_col])
            model = sm.OLS(species_var, X).fit()
            
            results.append({
                'Species': species_col,
                'Environment': env_col,
                'Spearman_r': corr,
                'Spearman_p': pval,
                'OLS_coef': model.params[1],
                'OLS_p': model.pvalues[1],
                'R_squared': model.rsquared,
                'EnvGroup': env_group_dict[env_col]  # Record the grouping of environment variables
            })
        
        all_results.extend(results)
    
    result_df = pd.DataFrame(all_results)
    
    # 5. Hierarchical FDR correction (Grouped by species + environmental variables)
    result_df['Spearman_q'] = np.nan
    for species_col in species.columns:
        for env_group in result_df['EnvGroup'].unique():
            mask = (result_df['Species'] == species_col) & (result_df['EnvGroup'] == env_group)
            pvals = result_df.loc[mask, 'Spearman_p'].values
            if len(pvals) > 0:  
                _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
                result_df.loc[mask, 'Spearman_q'] = qvals
    
    return result_df.sort_values(['Species', 'Spearman_p'])

def plot_robust_network(result_df, min_abs_corr=0, figsize=(15, 12), r_threshold=0):
    """Fix type errors and enhance the associated network diagram"""
    import networkx as nx
    
    # Ensure input is of numeric type
    if not pd.api.types.is_numeric_dtype(result_df['Spearman_r']):
        result_df['Spearman_r'] = pd.to_numeric(result_df['Spearman_r'], errors='coerce')
    
    # Create graph object

    # edges = result_df[
    #     (result_df['Spearman_q'] < 0.05) & 
    #     (abs(result_df['Spearman_r']) > r_threshold)
    # ]
    edges = result_df[
        result_df['Spearman_q'] < 0.05
    ]
    print(edges)

    if edges.empty:
        print(f"No significant associations with |r| > {r_threshold}")
        return
    G = nx.Graph()
    
    # Add nodes (distinguish between species and environment)
    species_nodes = edges['Species'].unique()
    env_nodes = edges['Environment'].unique()
    
    G.add_nodes_from(species_nodes, node_type='species', color='#4E79A7', size=800)
    G.add_nodes_from(env_nodes, node_type='environment', color='#E15759', size=1200)
    
    # Add edges (fix type errors)
    for _, row in edges.iterrows():
        try:
            r = float(row['Spearman_r'])
            if abs(r) >= min_abs_corr:
                G.add_edge(
                    row['Species'], 
                    row['Environment'],
                    weight=abs(r),
                    color='red' if r > 0 else 'blue',
                    style='solid' if r > 0 else 'dashed',
                    alpha=min(0.9, 0.3 + abs(r)*0.7)  # Transparency proportional to correlation strength
                )
        except (TypeError, ValueError) as e:
            print(f"Skipping invalid data: {row['Species']}-{row['Environment']}, Error: {e}")
    
    # Layout optimization
    # pos = nx.spring_layout(G, k=1.2, seed=42, iterations=100)
        # Use force-directed layout with optimized parameters
    pos = nx.spring_layout(
        G, 
        k=2,          # Node repulsion coefficient (reduce to decrease overlap)
        seed=42,        # Fix random seed
        iterations=100  # Increase number of iterations
    )
    
    # Drawing settings
    plt.figure(figsize=figsize)
    
    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=species_nodes,
        node_color=[G.nodes[n]['color'] for n in species_nodes],
        node_size=[G.nodes[n]['size'] for n in species_nodes],
        alpha=1
    )
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=env_nodes,
        node_color=[G.nodes[n]['color'] for n in env_nodes],
        node_size=[G.nodes[n]['size'] for n in env_nodes],
        alpha=1,
        node_shape='s'
    )
    
    # Draw edges
    for edge in G.edges(data=True):
        nx.draw_networkx_edges(
            G, pos,
            edgelist=[(edge[0], edge[1])],
            width=edge[2]['weight']*7,
            alpha=edge[2]['alpha'],
            edge_color=edge[2]['color'],
            style=edge[2]['style']
        )
    
    # Label settings (automatic line breaks)
    labels = {}
    for node in G.nodes():
        if len(node) > 10:
            labels[node] = '\n'.join([node[:8], node[8:]])
        else:
            labels[node] = node
    
    nx.draw_networkx_labels(
        G, pos, labels,
        font_size=18,
        font_family='Times New Roman',
        bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=1)
    )
    
    # Add professional legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Species',
              markerfacecolor='#4E79A7', markersize=18),
        Line2D([0], [0], marker='s', color='w', label='Environment',
              markerfacecolor='#E15759', markersize=18),
        Line2D([0], [0], color='red', lw=3, label='Positive'),
        Line2D([0], [0], color='blue', lw=3, linestyle='--', label='Negative')
    ]
    # plt.legend(handles=legend_elements, loc='upper right', framealpha=1)
    plt.legend(
    handles=legend_elements,
    loc='center left',          # Position changed to top-left
    framealpha=1,             # Opacity
    fontsize=18,              # Font size
    title='Legend',           # Legend title (optional)
    title_fontsize='18',      # Legend title font size (optional)
    bbox_to_anchor=(0, 1)     # Adjust position (optional, to prevent overlap)
    )
    plt.title(
        f"Species-Environment Association Network\n"
        f"(Correlation threshold: q <= 0.05)",
        pad=20, fontsize=20
    )
    
    plt.axis('off')
    plt.tight_layout()
    
    # Save the figure
    os.makedirs('results', exist_ok=True)
    plt.savefig(
        'results/robust_network.png',
        dpi=300,
        bbox_inches='tight',
        transparent=True
    )
    plt.close()
    print("Saved the repaired network diagram: results/robust_network.png")

def plot_all_associations(result_df, figsize=(20, 16), font_scale=1.2):
    """
    Display a heatmap of all species-environment associations
    Parameters:
        result_df - Dataframe containing all association results
        figsize - Figure size (width, height)
        font_scale - Font scaling factor
    """
    sns.set(font_scale=font_scale)
    
    # Create matrix data (including all associations)
    matrix = result_df.pivot(index='Species', 
                            columns='Environment', 
                            values='Spearman_r')
    
    # Sort by correlation strength
    row_order = result_df.groupby('Species')['Spearman_r'].mean().sort_values().index
    col_order = result_df.groupby('Environment')['Spearman_r'].mean().sort_values().index
    
    plt.figure(figsize=figsize)
    
    # Draw heatmap (add p-value star markers)
    ax = sns.heatmap(
        matrix.loc[row_order, col_order],
        annot=True, 
        cmap='coolwarm', 
        center=0,
        fmt=".2f",
        linewidths=1,
        cbar_kws={'label': 'Spearman Correlation (r)'},
        annot_kws={'size': font_scale*12}
    )
    ax.set_xlabel('Environment', fontsize=font_scale*16)  # x-axis label
    ax.set_ylabel('Species', fontsize=font_scale*16)      # y-axis label
    
    # **Modify color bar label font size**
    cbar = ax.collections[0].colorbar  # Get color bar object
    cbar.set_label('Spearman Correlation (r)', fontsize=font_scale*16)
    # Add significance markers (directly add stars in heatmap cells)
    for i, species in enumerate(row_order):
        for j, env in enumerate(col_order):
            row = result_df[(result_df['Species']==species) & 
                          (result_df['Environment']==env)]
            if not row.empty:
                pval = row['Spearman_p'].values[0]
                qval = row['Spearman_q'].values[0]
                
                if qval < 0.01:
                    ax.text(j+0.77, i+0.5, "***", 
                           ha='center', va='center', color='black', fontsize=font_scale*14, fontweight='bold')
                elif qval < 0.05:
                    ax.text(j+0.84, i+0.5, "**", 
                           ha='center', va='center', color='black', fontsize=font_scale*14, fontweight='bold')
                    # print("pval",pval)
                    print("qval",qval)
                elif qval < 0.1:
                    ax.text(j+0.77, i+0.5, "*", 
                           ha='center', va='center', color='black', fontsize=font_scale*14, fontweight='bold')
    
    plt.title("All Species-Environment Associations\n" +
             "(***: q < 0.01, **: q < 0.05, *: q < 0.1)", 
             pad=20, fontsize=font_scale*18, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=18)
    plt.yticks(rotation=0, fontsize=18)
    plt.tight_layout()
    
    # Save the figure
    os.makedirs('results', exist_ok=True)
    plt.savefig('results/all_associations_heatmap.png', 
               dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved all association heatmap: results/all_associations_heatmap.png")



def plot_key_scatters(species, env, result_df, n_plots=6):
    """Plot key scatter plots"""
    sig_df = result_df[
        (result_df['Spearman_q'] < 0.05) & 
        (abs(result_df['Spearman_r']) > 0.5)
    ].sort_values('Spearman_r', ascending=False)
    
    if sig_df.empty:
        print("No significant associations with |r| > 0.5")
        return
    
    # Create output directory
    os.makedirs('results/scatter_plots', exist_ok=True)
    
    # Plot the top N most significant relationships
    for i, (_, row) in enumerate(sig_df.head(n_plots).iterrows()):
        plt.figure(figsize=(6, 5))
        sns.regplot(
            x=env[row['Environment']], 
            y=species[row['Species']],
            scatter_kws={'s': 50, 'alpha': 0.7, 'color': 'steelblue'},
            line_kws={'color': 'crimson', 'lw': 1.5}
        )
        
        plt.title(
            f"{row['Species']} vs {row['Environment']}\n"
            f"Spearman r = {row['Spearman_r']:.2f} (q = {row['Spearman_q']:.3f})",
            pad=12, fontsize=18
        )
        plt.xlabel(row['Environment'], fontsize=18)
        plt.ylabel(row['Species'], fontsize=18)
        plt.tight_layout()
        
        # Save the figure
        fname = f"results/scatter_plots/{row['Species']}_vs_{row['Environment']}.png"
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        plt.close()

# =============================================
# 4. Report Generation Module
# =============================================
def generate_report(result_df, species, env):
    """Generate an HTML analysis report"""
    print("\n[Generating Analysis Report]")
    
    # Create result directories
    os.makedirs('results', exist_ok=True)
    os.makedirs('results/scatter_plots', exist_ok=True)  # Ensure scatter plot directory exists
    
    # Calculate summary statistics
    total_species = species.shape[1]
    total_envs = env.shape[1]
    sig_associations = (result_df['Spearman_q'] < 0.05).sum()
    
    # Get top 10 results
    top_results = result_df[
        (result_df['Spearman_q'] < 0.05) & 
        (abs(result_df['Spearman_r']) > 0.3)
    ].sort_values('Spearman_r', ascending=False).head(10)
    
    # Prepare data for the strongest correlations
    top_pos = ""
    top_pos_r = 0
    top_neg = ""
    top_neg_r = 0
    
    if not top_results.empty:
        top_pos = f"{top_results.iloc[0]['Species']} ~ {top_results.iloc[0]['Environment']}"
        top_pos_r = top_results.iloc[0]['Spearman_r']
        
        neg_results = result_df[result_df['Spearman_r'] < 0]
        if not neg_results.empty:
            top_neg_row = neg_results.sort_values('Spearman_r').iloc[0]
            top_neg = f"{top_neg_row['Species']} ~ {top_neg_row['Environment']}"
            top_neg_r = top_neg_row['Spearman_r']
    
    # HTML template (fix syntax errors)
    html_template = """
<!DOCTYPE html>
<html>
<head>
    <title>Species-environment Association Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; line-height: 1.6; margin: 0 auto; max-width: 1200px; padding: 20px }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px }
        h2 { color: #2980b9; margin-top: 30px }
        .summary { background-color: #f8f9fa; padding: 15px; border-radius: 5px }
        table { width: 100%; border-collapse: collapse; margin: 20px 0 }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd }
        th { background-color: #3498db; color: white }
        tr:nth-child(even) { background-color: #f2f2f2 }
        img { max-width: 100%; height: auto; display: block; margin: 20px auto }
        .highlight { font-weight: bold; color: #e74c3c }
    </style>
</head>
<body>
    <h1>Species-environment Association Analysis Report</h1>
    
    <div class="summary">
        <h2>Analysis Summary</h2>
        <p>- Analyze the number of species: <span class="highlight">{{ total_species }}</span></p>
        <p>- The number of environmental variables: <span class="highlight">{{ total_envs }}</span></p>
        <p>- Number of significant associations(q &lt; 0.05): <span class="highlight">{{ sig_associations }}</span></p>
        <p>- Strongest positive correlation: <span class="highlight">{{ top_pos }}</span> (r = {{ "%.2f"|format(top_pos_r) }})</p>
        <p>- Strongest negative correlation: <span class="highlight">{{ top_neg }}</span> (r = {{ "%.2f"|format(top_neg_r) }})</p>
    </div>
    
    <h2>Top 10 Significant association</h2>
    {{ top_table }}
    
    <h2>Associated heat map</h2>
    <img src="all_associations_heatmap.png" alt="Associated heat map">
    
    <h2>Associated network</h2>
    <img src="robust_network.png" alt="Associated network figure">
    
    <h2>Scatter plot of key relationships</h2>
    {% for img in scatter_imgs %}
    <div style="margin-bottom: 30px;">
        <img src="{{ img }}" alt="Scatter plot of key relationships">
    </div>
    {% endfor %}
</body>
</html>
    """
    
    # Prepare scatter plot list
    scatter_imgs = []
    if os.path.exists('results/scatter_plots'):
        scatter_imgs = [
            f"scatter_plots/{f}" for f in os.listdir('results/scatter_plots') 
            if f.endswith('.png')
        ][:6]

    # Render the report
    from jinja2 import Template
    template = Template(html_template)

    # Ensure top_results has data
    if top_results.empty:
        top_table = "<p>There were no significant associated results</p>"
    else:
        top_table = top_results.to_html(index=False)

    report = template.render(
        total_species=total_species,
        total_envs=total_envs,
        sig_associations=sig_associations,
        top_pos=top_pos,
        top_pos_r=top_pos_r,
        top_neg=top_neg,
        top_neg_r=top_neg_r,
        top_table=top_table,
        scatter_imgs=scatter_imgs
    )

    # Save the report
    with open('results/report.html', 'w', encoding='utf-8') as f:
        f.write(report)

# =============================================
# Main Program
# =============================================
if __name__ == "__main__":
    # Data loading
    species_data, env_data = load_and_preprocess(
        species_path="level-2-59.csv",
        env_path="environmental_data_500m.csv"
    )
    
    # Statistical analysis
    analysis_results = analyze_dataset(species_data, env_data)
    
    # Save complete results
    os.makedirs('results', exist_ok=True)
    analysis_results.to_csv("results/full_analysis_results.csv", index=False)
    
    # Visualization
    plot_all_associations(analysis_results)
    # plot_network(analysis_results)
    plot_robust_network(analysis_results)
    plot_key_scatters(species_data, env_data, analysis_results)
    
    # Generate report
    generate_report(analysis_results, species_data, env_data)
    
    print("\n[Analysis Complete]")
    print(f"Results saved as: {os.path.abspath('results')}")