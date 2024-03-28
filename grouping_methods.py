import numpy as np
import pandas as pd
import itertools

import seaborn as sns
import matplotlib.pyplot as plt 
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering 

from sklearn.metrics import mean_squared_error
from sklearn.decomposition import NMF 
import matplotlib.pyplot as plt 


def get_combinatorial_grouping(matrix_prisma):
    matrix_prisma = matrix_prisma.melt(id_vars = ['gene'], var_name = 'site', value_name = 'intensity')
    matrix_prisma = matrix_prisma.loc[matrix_prisma['intensity'] > 0.9]
    matrix_prisma.reset_index(inplace = True, drop = True)
    matrix_prisma

    # Generate the dataframe of binding sites for each gene
    site_per_gene = matrix_prisma.groupby('gene')['site'].agg(';'.join).reset_index()
    site_per_gene
    # site_per_gene.to_csv...

    # Generate the dataframe of all combinations of binding sites for each gene
    site_comb_per_gene = pd.DataFrame({'site': [], 'gene': []})
    def expand_site(row):
        ##Returns a df of unique binding-site subsets for each protein
        ## e.g..  "ABCC1" which has binding sites: "P74"  "P171" "P227" would result in:
        #        combo protein
        # 1 P171, P227   ABCC1
        # 2  P171, P74   ABCC1
        # 3  P227, P74   ABCC1
        # 4       P171   ABCC1
        # 5       P227   ABCC1
        # 6        P74   ABCC1
        sites = row['site']
        gene = row['gene']
        parts = sites.split(';')
        df_part = pd.DataFrame({'site': [parts], 'gene': gene})
        for l in range(1, len(parts)):
            for sub in itertools.combinations(parts, l):
                sub = list(sub)
                df_sub = pd.DataFrame({'site': [sub], 'gene': gene})
                df_part = pd.concat([df_part, df_sub])
        return df_part

    for r in site_per_gene.apply(axis = 1, func = expand_site):
        site_comb_per_gene = pd.concat([site_comb_per_gene, r])

    site_comb_per_gene

    site_comb_per_gene['site'] = site_comb_per_gene['site'].apply(','.join)
    site_comb_per_gene = site_comb_per_gene.sort_values('site', ascending = True)
    site_comb_per_gene.reset_index(inplace = True, drop = True)

    # Generate the dataframe of all combinations of binding sites for grouped genes
    site_comb_gene_group = site_comb_per_gene.groupby('site')['gene'].agg(';'.join).reset_index()
    site_comb_gene_group

    # Add a column with the count of genes for each group of sites
    site_comb_gene_group['gene'] = site_comb_gene_group['gene'].str.split(';')
    site_comb_gene_group['gene_count'] = site_comb_gene_group['gene'].str.len()
    site_comb_gene_group['gene'] = site_comb_gene_group['gene'].apply(','.join)
    site_comb_gene_group

    # Filter for sites with more than one protein binding
    site_comb_gene_group = site_comb_gene_group[site_comb_gene_group['gene_count'] > 1]

    # Add a column of grouping index
    site_comb_gene_group.index = np.arange(1, len(site_comb_gene_group) + 1)
    site_comb_gene_group = site_comb_gene_group.reset_index(names = 'group')
    site_comb_gene_group = site_comb_gene_group.rename(columns = {'gene':'proteins', 'gene_count':'protein_count'})
    return site_comb_gene_group


def get_hierachical_grouping(matrix_prisma):
    matrix_prisma = matrix_prisma.set_index(matrix_prisma['gene'])
    matrix_prisma = matrix_prisma.drop(columns = 'gene')
    matrix_prisma

    sns.clustermap(matrix_prisma)

    # Using the dendrogram to find the optimal numbers of clusters
    dendrogram = sch.dendrogram(sch.linkage(matrix_prisma, method  = "ward", metric = 'euclidean'))
    # The Ward method is a method that attempts to minimize the within-cluster variance,
    # which is similar to K-means to minimize the within-cluster sum of square (wcss)
    plt.title('Dendrogram')
    plt.xlabel('Genes')
    plt.ylabel('Euclidean distances')
    plt.show()
    # Then find the longest vertical distance we can without crossing any horizontal lines,
    # and count the lines on the diagram and figure out how many clusters are best
    # which should be 4 here
    # AgglomerativeClustering(n_clusters = 4...

    # Fitting hierarchical clustering 
    # with Agglomerative Hierarchical Clustering algorithm using Euclidean distance and ward method for our algorithm class  
    hc = AgglomerativeClustering(n_clusters = None, distance_threshold = 3, linkage ='ward')
    hc

    # fit the hierarchical clustering algorithm to the dataset while creating the 
    # clusters vector that tells for each gene which cluster the gene belongs to
    y_hc = hc.fit_predict(matrix_prisma)

    len(y_hc)
    y_hc
    max(y_hc)

    # Generate the table of genes with their cluster groups
    matrix_prisma.insert(0, 'group', y_hc + 1)
    matrix_prisma

    gene_cluster = matrix_prisma[['group']].sort_values('group')
    gene_cluster = gene_cluster.reset_index(names = 'gene')
    gene_cluster = gene_cluster.groupby('group')['gene'].agg(','.join).reset_index()
    gene_cluster['gene_count'] = gene_cluster['gene'].str.count(',') + 1

    # Filter for sites with more than one protein binding
    gene_cluster = gene_cluster[gene_cluster['gene_count'] > 1]

    gene_cluster = gene_cluster.rename(columns = {'gene':'proteins', 'gene_count':'protein_count'})
    return gene_cluster


def get_nmf_grouping(matrix_prisma):
    print(matrix_prisma)
    matrix_prisma.set_index(matrix_prisma['gene'], inplace = True)
    matrix_prisma.drop(columns = 'gene', inplace = True)
    matrix_prisma

    # Calculate the optimal rank (n of components) of matrix V 
    # (when cophenetic score hits first local minimum)
    def calculate_rank(matrix_V):
        # Calculate the Frobenius Norm for the matrix V
        norm_matrix = np.linalg.norm(matrix_V, ord = 'fro')
        # Define the benchmark for the root mean squared error (RMSE)
        benchmark = norm_matrix * 0.0001
        # Iterate through ranks > 3 to find the optimal rank
        rank = 3
        while True:
            # Initialise the model:
            model = NMF(n_components = rank, init = 'random', random_state = 0, max_iter = 1000)
            W = model.fit_transform(matrix_V)
            H = model.components_
            V = W @ H
            # Calculate RMSE of each-iteration NMF model and the matrix V
            RMSE = np.sqrt(mean_squared_error(matrix_V, V))
            # Return the rank of the NMF model with a small enough RMSE and the model
            if RMSE < benchmark: # if only use this condition, the resulting rank is 220.
                return rank
            # Otherwise increase the rank by 1 and re-calculate
            else:
                rank += 1

    # optimal_rank = calculate_rank(matrix_prisma)
    # optimal_rank


    # Hardcode the resulting optimal rank so don't have to wait for the calculation each time
    optimal_rank  = 220


    # Decompose the matrix using NMF
    # default random_state = 0?
    nmf_model = NMF(n_components = optimal_rank)
    nmf_model.fit(matrix_prisma)

    # Get the coefficient matrix H
    # it represents how much each site belongs to the signature group
    H = pd.DataFrame(nmf_model.components_, columns = matrix_prisma.columns)
    H.index = H.index + 1
    H
    sns.heatmap(H)
    sns.clustermap(H) # A line
    # Normalise matrix H by the columns - for each site
    H.max()
    H_normalized = (H - H.min()) / (H.max() - H.min())
    H_normalized
    # H_normalized.to_csv('/home/ruoci/Desktop/PRISMA papers/code and data/ready_parts/out/H_normalized.csv')
    H_normalized.max()
    # Roughly visulise the normalised coefficient matrix H
    sns.heatmap(H_normalized)
    sns.clustermap(H_normalized)
    # Filtering for the sites with top 1-3 non-zero intensity for each signature group
    H_filtered = H_normalized[H_normalized > 0]
    H_filtered
    H_filtered.count()
    plt.boxplot(H_filtered.count())

    H_filtered = H_filtered.apply(pd.Series.nlargest, n = 3, axis = 0)
    H_filtered.count()
    plt.boxplot(H_filtered.count())

    # Generate the table with top sites, their signature group and their intensity
    top_sites = pd.DataFrame(H_filtered)
    top_sites['group'] = top_sites.index
    top_sites = top_sites.melt(id_vars = 'group', var_name = 'site')
    top_sites = top_sites.dropna()
    top_sites.reset_index(inplace = True, drop = True)
    top_sites = top_sites[['group', 'site']]
    top_sites

    sig_sites = top_sites.groupby('group')['site'].agg(';'.join).reset_index()
    sig_sites

    # Get the feature matrix W
    # it represents how much each gene contributes to each signature
    W = pd.DataFrame(nmf_model.fit_transform(matrix_prisma), index = matrix_prisma.index)
    W.columns = W.columns + 1
    W
    W.max()
    max(W.max())
    sns.heatmap(W)
    sns.clustermap(W)
    # Normalise the feature matrix W by the rows - for each gene
    W_transposed = W.T
    W_transposed
    W_normalized = (W_transposed - W_transposed.min()) / (W_transposed.max() - W_transposed.min())
    W_normalized
    W_normalized.mean()
    # Roughly visulise the normalised feature matrix W
    sns.heatmap(W_normalized)
    sns.clustermap(W_normalized)
    # Filtering for the genes with top 1-3 non-zero intensity for each signature group
    W_filtered = W_normalized[W_normalized > 0]
    W_filtered
    W_filtered.count()
    plt.boxplot(W_filtered.count())

    W_filtered = W_filtered.apply(pd.Series.nlargest, n = 3, axis = 0)
    W_filtered
    W_filtered.count()
    plt.boxplot(W_filtered.count())


    # Generate the table with top genes, their signature group and their intensity
    top_genes = pd.DataFrame(W_filtered)
    top_genes['group'] = top_genes.index
    top_genes = top_genes.melt(id_vars = 'group', var_name = 'gene')
    top_genes = top_genes.dropna()
    top_genes.reset_index(inplace = True, drop = True)
    top_genes = top_genes[['group', 'gene']]
    top_genes

    sig_genes = top_genes.groupby('group')['gene'].agg(','.join).reset_index()
    sig_genes

    # Generate a combined table of the signature group with the top sites and top genes
    df_top = sig_sites.merge(sig_genes)
    df_top['gene_count'] = df_top['gene'].str.count(',') + 1
    df_top = df_top.rename(columns = {'gene':'proteins', 'gene_count':'protein_count'})
    df_top

    # Filter for sites with more than one protein binding
    df_top = df_top[df_top['protein_count'] > 1]

    return(df_top)