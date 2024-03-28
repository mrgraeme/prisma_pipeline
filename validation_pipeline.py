

import pandas as pd
import itertools
from functools import reduce


def validate_groups(grouped_protein, df_corum, df_corr, df_string_link, df_string_info):

    # Generate protein pairs from the grouping file
    grouped_protein

    protein_pair = pd.DataFrame({'group': [], 'source': [], 'target': []})
    def get_protein_pair(row):
        df_pair = pd.DataFrame({'group': [], 'source': [], 'target': []})
        proteins = row['proteins'].split(',')
        if len(proteins) > 1:
            df_pair = pd.DataFrame(list(itertools.combinations(proteins, 2)))
            df_pair.columns = ['source', 'target']
            df_pair['group'] = str(row['group'])
        return df_pair

    for r in grouped_protein.apply(axis = 1, func = get_protein_pair):
        protein_pair = pd.concat([protein_pair, r])

    protein_pair

    # Sort the protein pairs alphabetically and drop duplicates
    def sort_protein_pair(row):
        temp_source = min(row['source'], row['target'])
        temp_target = max(row['source'], row['target'])
        row['source'] = temp_source
        row['target'] = temp_target
        return row

    protein_pair = protein_pair.apply(axis = 1, func = sort_protein_pair)
    protein_pair = protein_pair.drop_duplicates()
    protein_pair


    # CORUM (manually curated mammalian protein complexes) for interaction validation:
    df_corum = df_corum.dropna()
    df_corum
    # Assign the weighted score of 3 for CORUM interactions
    df_corum['corum'] = 3
    df_corum = df_corum[['source', 'target', 'corum']]
    # Make sure all gene names in columns are in uppercase
    df_corum['source'] = df_corum['source'].apply(lambda x: x.upper())
    df_corum['target'] = df_corum['target'].apply(lambda x: x.upper())
    # Sort the pairs alphabetically and drop duplicates
    df_corum = df_corum.apply(axis = 1, func = sort_protein_pair)
    df_corum = df_corum.drop_duplicates()
    df_corum


    # Paired correlation data based on protein abundance across all CCLE cell lines for interaction validation:
    df_corr = df_corr.dropna()
    df_corr
    # Filter for pairs with pearson correlation > 0.5
    df_corr = df_corr[df_corr['pearsons_correlaton'] > 0.5]
    # Assign the weighted score of 1 for correlated interactions
    df_corr['corr'] = 1
    df_corr = df_corr[['source', 'target', 'corr']]
    # Make sure all gene names in columns are in uppercase
    df_corr['source'] = df_corr['source'].apply(lambda x: x.upper())
    df_corr['target'] = df_corr['target'].apply(lambda x: x.upper())
    # Sort the pairs alphabetically and drop duplicates
    df_corr = df_corr.apply(axis = 1, func = sort_protein_pair)
    df_corr = df_corr.drop_duplicates()
    df_corr


    # String (functional protein association networks) for interaction validation:
    # Load the detailed table for protein network data in Homo Sapiens from string
    df_string_link = df_string_link.dropna()
    df_string_link
    # Select and rename the columns needed
    df_string_link = df_string_link[['protein1', 'protein2', 'experimental']]
    df_string_link = df_string_link.rename(columns = {'protein1': 'source', 'protein2': 'target'})
    # Filter for pairs with experimental > 150
    df_string_link = df_string_link[df_string_link['experimental'] > 150]
    # Sort the pairs alphabetically and drop duplicates
    df_string_link = df_string_link.apply(axis = 1, func = sort_protein_pair)
    df_string_link = df_string_link.drop_duplicates()
    df_string_link

    # Load the accessory table for list of STRING proteins incl. their display namesnetwork data in Homo Sapiens from string
    df_string_info = df_string_info.dropna()
    # Select and rename the columns needed
    df_string_info = df_string_info[['#string_protein_id', 'preferred_name']]
    df_string_info = df_string_info.rename(columns={'#string_protein_id': 'id', 'preferred_name': 'name'})
    # Make sure all gene names in columns are in uppercase
    df_string_info['name'] = df_string_info['name'].apply(lambda x: x.upper())
    df_string_info

    # Convert the id in the link table into its corresponding name in the info table
    df_string = df_string_link.merge(df_string_info, how = 'left', left_on = 'source', right_on = 'id')
    df_string = df_string[['name', 'target', 'experimental' ]]
    df_string = df_string.rename(columns = {'name': 'source'})
    df_string = df_string.merge(df_string_info, how = 'left', left_on = 'target', right_on = 'id')
    df_string = df_string[['source', 'name', 'experimental' ]]
    df_string = df_string.rename(columns = {'name': 'target'})
    df_string

    # Filter for pairs with experimental > 400
    df_string_high = df_string[df_string['experimental'] > 400]
    df_string_high

    # Assign the weighted score of 1 for string experimental score > 150
    df_string['string_exp_150'] = 1
    df_string = df_string[['source', 'target', 'string_exp_150']]
    df_string

    # Assign the weighted score of 2 for string experimental score > 400
    df_string_high['string_exp_400'] = 2
    df_string_high = df_string_high[['source', 'target', 'string_exp_400']]
    df_string_high

    # Merge the PRISMA pairs with the interaction validations together:
    protein_pair_validation = protein_pair.merge(df_corum, on = ['source', 'target'], how = 'left')
    protein_pair_validation = protein_pair_validation.merge(df_corr, on = ['source', 'target'], how = 'left')
    protein_pair_validation = protein_pair_validation.merge(df_string_high, on = ['source', 'target'], how = 'left')
    protein_pair_validation = protein_pair_validation.merge(df_string, on = ['source', 'target'], how = 'left')
    protein_pair_validation = protein_pair_validation.fillna(0)
    protein_pair_validation

    # Calculate the overall score for each pair
    protein_pair_validation['overall_score']  = protein_pair_validation[['corum', 'corr', 'string_exp_400', 'string_exp_150']].apply(axis = 1, func = max)
    protein_pair_validation

    # Calculate the average score for each group
    group_average_validation = pd.DataFrame()
    group_average_validation['protein_pair_count'] = protein_pair_validation.groupby('group')[['source']].count()
    group_average_validation['mean_score'] = protein_pair_validation.groupby('group')['overall_score'].mean()
    group_average_validation['std_score'] = protein_pair_validation.groupby('group')['overall_score'].std()
    group_average_validation['corum_coverage'] = protein_pair_validation.groupby('group')['corum'].sum() / 3 / group_average_validation['protein_pair_count']
    group_average_validation['string_exp_400_coverage'] = protein_pair_validation.groupby('group')['string_exp_400'].sum() / 2 / group_average_validation['protein_pair_count']
    group_average_validation['string_exp_150_coverage'] = protein_pair_validation.groupby('group')['string_exp_150'].sum() / group_average_validation['protein_pair_count']
    group_average_validation['corr_coverage'] = protein_pair_validation.groupby('group')['corr'].sum() / group_average_validation['protein_pair_count']
    group_average_validation = group_average_validation.fillna(0)

    group_average_validation = group_average_validation.reset_index()
    group_average_validation['group'] = group_average_validation['group'].apply(int)
    group_average_validation

    return (protein_pair_validation, group_average_validation)

