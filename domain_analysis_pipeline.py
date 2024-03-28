
import pandas as pd
from scipy.stats import fisher_exact
import numpy as np

def preprocess_interpro(interpro_raw, id_map):
    interpro_anno = interpro_raw.iloc[:,0:3]
    interpro_anno.columns = ['UniProt ID', 'Interpro ID', 'domain']
    interpro_anno.dropna(inplace = True)
    id_map = id_map.iloc[:,1:3]
    id_map.columns = ['proteins', 'UniProt ID']
    interpro_anno = interpro_anno.merge(id_map, on='UniProt ID', how='left').dropna()
    interpro_anno = interpro_anno.drop_duplicates()
    return interpro_anno

def domain_analysis(grouped_protein, domain_anno):

    # Load the dataframe of lists of proteins with groups of binding sites
    grouped_protein['proteins'] = grouped_protein['proteins'].str.split(',')

    # Generate a long table of groups of binding sites for each protein
    grouped_protein = grouped_protein.explode('proteins')

    # Merge the protein binding site table with the domain annotation table
    protein_domains = grouped_protein.merge(domain_anno,  how='left')
    protein_domains = protein_domains.drop_duplicates()
    protein_domains['group'] = protein_domains['group'].apply(str)
    protein_domains

    # Summerise the domains within each binding site group
    group_domains = pd.DataFrame()
    group_domains['domain_protein_in_group_count'] = protein_domains.groupby(['group', 'protein_count', 'domain'])['proteins'].count()
    group_domains = group_domains.reset_index()
    group_domains = group_domains.rename(columns = {'protein_count': 'group_protein_count'})
    group_domains

    # Calculate the ratio of protein counts with each domain within each binding site group
    group_domains['coverage_of_domain_proteins_in_group'] = group_domains['domain_protein_in_group_count'] / group_domains['group_protein_count']
    group_domains['group'] = group_domains['group'].apply(int)
    group_domains = group_domains.sort_values(['group', 'domain_protein_in_group_count'], ascending = True)
    group_domains['group'] = group_domains['group'].apply(str)
    group_domains

    # Calculate the statistics of each type of domain across binding site groups
    domain_summary= pd.DataFrame()
    domain_summary['total_domain_protein_count'] = group_domains.groupby('domain')['domain_protein_in_group_count'].sum()
    domain_summary['groups'] =  group_domains.groupby('domain')['group'].apply('; '.join)
    domain_summary['group_count'] =  group_domains.groupby('domain')['domain_protein_in_group_count'].count()
    domain_summary['average_domain_protein_count_per_group'] = domain_summary['total_domain_protein_count'] / domain_summary['group_count']
    domain_summary['mean_coverage'] = group_domains.groupby('domain')['coverage_of_domain_proteins_in_group'].mean()
    domain_summary['std_coverage'] = group_domains.groupby('domain')['coverage_of_domain_proteins_in_group'].std().fillna(0)
    domain_summary.sort_values(['mean_coverage', 'group_count', ], ascending = False)
    domain_summary = domain_summary.reset_index()
    domain_summary

    # Add the column of representation_of representation_of_total_domain_proteins
    group_domains = group_domains.merge(domain_summary, left_on = 'domain', right_on = 'domain', how = 'left') 
    group_domains = group_domains[['group', 'group_protein_count', 'domain', 'total_domain_protein_count', 'domain_protein_in_group_count', 'coverage_of_domain_proteins_in_group']]
    group_domains['representation_of_total_domain_proteins'] = group_domains['domain_protein_in_group_count'] / group_domains['total_domain_protein_count']
    group_domains

    return (protein_domains, group_domains, domain_summary)


def domain_enrichment(grouped_protein, domain_anno):
    
    # Generate a long table of groups of binding sites for each protein
    grouped_protein['proteins'] = grouped_protein['proteins'].str.split(',')
    grouped_protein = grouped_protein.explode('proteins')
    grouped_protein

    # Merge the protein binding site table with the domain annotation table
    protein_domains = grouped_protein.merge(domain_anno, how = 'inner')
    protein_domains

    protein_domains = protein_domains[['group', 'proteins', 'protein_count', 'Interpro ID', 'domain']]
    # test_domain = test_domain.groupby([['group', 'Interpro ID']])['proteins'].agg(';'.join).reset_index()
    protein_domains = protein_domains.drop_duplicates()

    # Calculate the total number of proteins 
    total_protein_count = len(grouped_protein['proteins'].unique())
    total_protein_count 

    # Generate the table of domain summary table for all proteins included in the PRISMA exp
    across_group = protein_domains[['proteins', 'Interpro ID']]
    across_group = protein_domains.groupby('Interpro ID')['proteins'].agg(';'.join).reset_index()
    across_group['domain_count'] = across_group['proteins'].str.count(';') + 1
    across_group


    # Generate a table with the domains from PRISMA and useful info including domain counts in each group
    test_domain = pd.DataFrame()
    for group in protein_domains['group'].unique():
        each_group = pd.DataFrame()
        each_group = protein_domains.loc[protein_domains['group'] == group]
        each_group = each_group[['proteins', 'Interpro ID', 'protein_count']]
        each_group = each_group.groupby(['Interpro ID', 'protein_count'])['proteins'].agg(';'.join).reset_index()
        each_group['group'] = group
        each_group['domain_count'] = each_group['proteins'].str.count(';') + 1
        # print(f'group protein count:{group_protein_count}')
        # print(f'the group table: {each_group}')
        test_domain = pd.concat([test_domain, each_group], axis = 0)

    test_domain

    domain_enrichment_results = pd.DataFrame()

    def domain_enrichment_ingroup(row):
        tested_domain = row['Interpro ID']
        # print(tested_domain)
        group_domain_count = int(row['domain_count'])
        group_protein_count = int(row['protein_count'])
        total_domain_count = int(across_group.loc[across_group['Interpro ID'] == tested_domain]['domain_count'].iloc[0])
        other_group_domain_count = total_domain_count - group_domain_count

        test_table = np.array([[group_domain_count, group_protein_count - group_domain_count],
                            [other_group_domain_count, total_protein_count - group_protein_count - other_group_domain_count]])
        oddsr, p = fisher_exact(table = test_table, alternative = 'two-sided')
        # print(p)

        result = pd.DataFrame()
        result['Interpro ID'] = [tested_domain]
        result['group'] = row['group']
        result['proteins'] = row['proteins']
        result['p_value'] = [p]
        result['odds_ratio'] = [oddsr]
        # The odds ratio indicates that the odds of getting the domain in a group 
        # is how many times that of not getting the domain in that group.
        result['domain_in_group'] = [group_domain_count]
        result['total_domain'] = [total_domain_count]
        result['protein_in_group'] = [group_protein_count]
        result['total_protein'] = [total_protein_count]

        return result

    for r in test_domain.apply(axis = 1, func = domain_enrichment_ingroup):
        domain_enrichment_results = pd.concat([domain_enrichment_results, r])

    domain_enrichment_results.sort_values(by = ['p_value', 'odds_ratio'], axis = 0,
                                        ascending = [True, False], inplace = True,
                                        ignore_index = True)

    domain_enrichment_results
    protein_domains[['Interpro ID', 'domain']]

    domain_enrichment_results = domain_enrichment_results.merge(protein_domains[['Interpro ID', 'domain']], on='Interpro ID').drop_duplicates()
    domain_enrichment_results
    # Filter for the p_value < 0.01
    domain_enrichment_results = domain_enrichment_results[domain_enrichment_results['p_value'] < 0.01]

    ## Drop duplicates that only differ by Interpro ID and keep the first occurrence
    # domain_enrichment_results = domain_enrichment_results.astype('str')

    # domain_enrichment_results.drop_duplicates(subset = ['group', 'proteins', 'p_value', 'odds_ratio', 'domain_in_group',
    #                                                      'total_domain', 'protein_in_group', 'total_protein'],
    #                                                      keep = 'first', inplace = True)
    # domain_enrichment_results
    # can't do this because can't differentiate redundant domain variants and truly different domains

    # domain_enrichment_results.to_csv(domain_enrichment_path, index = False)

    domain_enrichment_results = domain_enrichment_results.astype('str')
    grouped_domain_enrichment = domain_enrichment_results.groupby('group')[['domain', 'Interpro ID']].agg('//'.join).reset_index()
    grouped_domain_enrichment['domain_count'] = grouped_domain_enrichment['domain'].str.count('//') + 1
    grouped_domain_enrichment

    # grouped_domain_enrichment.to_csv(grouped_domain_enrichment_path, index = False)

    return (domain_enrichment_results, grouped_domain_enrichment)



