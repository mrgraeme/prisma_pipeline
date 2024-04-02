import pandas as pd
import os


## Define relative and absolute paths to read and write files
project_data_dir = './project_data/'
global_data_dir = './global_data/'

# Read in the ELM SLiM prediction result for ARID1A
ARID1A_pred = pd.read_csv('%s/elm_arid1a_predicted_slims.csv' % global_data_dir)
ARID1A_pred.dropna(inplace = True)
ARID1A_pred

# Remove the position column suffix of author notation
ARID1A_pred['Positions'] = ARID1A_pred['Positions'].str.removesuffix(' [A]')

# Split the position column into start and end position columns
ARID1A_pred[['Start', 'End']] = ARID1A_pred['Positions'].str.split('-', expand = True)
ARID1A_pred.drop(columns = ['Positions'], inplace = True)
ARID1A_pred

# Create a dataframe for PRISMA site P1-P228
PRISMA_site = [f'P{i + 1}' for i in range(228)]
df_site = pd.DataFrame(PRISMA_site, columns = ['Site'])
df_site['Start'] = [i * 10 for i in range(228)]
df_site['End'] = [i * 10 + 20 for i in range(228)]
df_site


pred_slim_site = pd.DataFrame()
def find_slim_site(row):
    row_site = []
    start_pos = int(row['Start'])
    for ind in df_site.index:
        start_prisma = int(df_site['Start'][ind])
        end_prisma = int(df_site['End'][ind])
        site_prisma = df_site['Site'][ind]
        if (start_pos > start_prisma) & (start_pos < end_prisma):
            row_site.append(site_prisma)
    result = pd.DataFrame({'Elm Name': [row['Elm Name']],
                           'Instances': [row['Instances']],
                           'Start': [row['Start']],
                           'End': [row['End']]})
    result['Site'] = [row_site]
    return result

for r in ARID1A_pred.apply(axis = 1, func = find_slim_site):
    pred_slim_site = pd.concat([pred_slim_site, r])

pred_slim_site

pred_slim_site_path = '%s/elm_arid1a_predicted_slims_with_sites.csv' % project_data_dir
pred_slim_site.to_csv(pred_slim_site_path, index = False)


### Then match the SLiM sites with the group from different grouping methods
# Prompt for the grouping method prefix input 
method_prefix = input('Enter the grouping method used: '.strip()) or 'combinatorial'

grouped_protein_path = project_data_dir + 'ARID1A_' + method_prefix + '_grouped_protein.csv'
domain_enrichment_path = project_data_dir + 'ARID1A_' + method_prefix + '_grouped_domain_enrichment.csv'
grouped_enrichment_path = project_data_dir + 'ARID1A_' + method_prefix + '_grouped_func_enrichment.csv'

predicted_groups_path = project_data_dir + 'ARID1A_' + method_prefix + '_predicted_groups.csv'
predicted_groups_with_enrichments_path = project_data_dir + 'ARID1A_' + method_prefix + '_predicted_groups_with_enrichments.csv'

grouped_protein = pd.read_csv(grouped_protein_path)
grouped_protein['site'] = grouped_protein['site'].apply(lambda x: x.split(','))
grouped_protein

predicted_groups = pd.DataFrame()
def find_slim_group(row):
    result = pd.DataFrame()
    sites = row['Site']
    for slim_site in sites:
        for ind in grouped_protein.index:
            if slim_site in grouped_protein['site'][ind]:
                result['group'] = [grouped_protein['group'][ind]]
                result['group_site'] = [grouped_protein['site'][ind]]
                result['slim_site'] = [slim_site]
                result['Elm_name'] = row['Elm Name']
                result['proteins'] = grouped_protein['proteins'][ind] 
                result['protein_count'] = grouped_protein['protein_count'][ind]
    return result

for r in pred_slim_site.apply(axis = 1, func = find_slim_group):
    predicted_groups = pd.concat([predicted_groups, r])

predicted_groups.reset_index(inplace = True, drop = True)
predicted_groups['group_site'] = predicted_groups['group_site'].apply(lambda x: ','.join(x))
predicted_groups.drop_duplicates(inplace = True)
predicted_groups

predicted_groups.to_csv(predicted_groups_path, index = False)


### Add functional enrichment terms to each group
grouped_enrichment = pd.read_csv(grouped_enrichment_path)
grouped_enrichment

predicted_groups_func = predicted_groups.merge(grouped_enrichment, how = 'left')
predicted_groups_func


### Add domain enrichment info to each group
grouped_domain_enrichment = pd.read_csv(domain_enrichment_path)

# Merge with groups and keeps the domain name and interpro ID
predicted_groups_with_enrichments = predicted_groups_func.merge(grouped_domain_enrichment, how = 'left')
predicted_groups_with_enrichments

predicted_groups_with_enrichments.to_csv(predicted_groups_with_enrichments_path, index = False)

