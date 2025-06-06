import pandas as pd
import numpy as np
from scipy.stats import entropy
import skbio.diversity.alpha as alpha

def calculate_diversity_metrics(abundance_df, rarefy_depth=None):
    """
    Improved version of diversity metrics calculation, supports count data
    
    Parameters:
        abundance_df: DataFrame, rows as samples, columns as species, values as counts
        rarefy_depth: int or None, if specified, will perform rarefaction
    """
    # Data preprocessing
    df = abundance_df.copy()
    if rarefy_depth:
        df = rarefy_data(df, rarefy_depth)
    
    # Calculate total reads
    total_counts = df.sum(axis=1)
    if (total_counts == 0).any():
        raise ValueError("Some samples have zero total counts, please check the data")
    
    # Calculate relative abundance
    rel_abundance = df.div(total_counts, axis=0)
    
    # Initialize result DataFrame
    results = pd.DataFrame(index=df.index)
    
    # 1. Richness
    results['Richness'] = (df > 0).sum(axis=1)
    
    # 2. Shannon index
    results['Shannon'] = rel_abundance.apply(
        lambda x: -sum(p * np.log(p) for p in x[x > 0]), axis=1)
    
    # 3. Simpson index
    results['Simpson'] = 1 - (rel_abundance**2).sum(axis=1)
    
    # 4. Pielou evenness
    with np.errstate(divide='ignore', invalid='ignore'):
        results['Pielou_evenness'] = results['Shannon'] / np.log(results['Richness'])
    results['Pielou_evenness'].replace([np.inf, -np.inf], np.nan, inplace=True)
    
    # 5-7. Count-based metrics
    results['Chao1'] = df.apply(lambda x: alpha.chao1(x[x > 0].astype(int)), axis=1)
    results['ACE'] = df.apply(lambda x: alpha.ace(x[x > 0].astype(int)), axis=1)
    results['Fisher_alpha'] = df.apply(lambda x: alpha.fisher_alpha(x[x > 0].astype(int)), axis=1)
    
    return results

def rarefy_data(df, depth=None):
    """Rarefy data to specified sequencing depth"""
    if depth is None:
        depth = df.sum(axis=1).min()
    return df.apply(lambda x: np.random.multinomial(depth, x/x.sum()), axis=1)

if __name__ == "__main__":
    input_file = "level-2-59.csv"
    output_file = "diversity_metrics_improved.csv"
    
    try:
        # Read data
        data = pd.read_csv(input_file, index_col=0)
        print(f"Successfully read data, containing {len(data)} samples, {len(data.columns)} species")
        
        # Calculate diversity metrics (optional rarefaction)
        diversity_df = calculate_diversity_metrics(data)  # Without rarefaction
        # diversity_df = calculate_diversity_metrics(data, rarefy_depth=10000)  # Rarefy to 10,000 reads
        
        # Save results
        diversity_df.to_csv(output_file)
        print(f"\nResults saved to {output_file}")
        print("\nResults preview:")
        print(diversity_df.head())
        
    except Exception as e:
        print(f"Error: {str(e)}")