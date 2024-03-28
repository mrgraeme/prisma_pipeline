
import pandas as pd

prisma_protein_label  = 'ARID1A'
grouping_method = 'combinatorial'

# Define relative and absolute paths to read and write files
project_data_dir = '~/work/data/prisma_analysis/arid1a'
global_data_dir = '~/work/data/prisma_analysis/data'
                                                     
# Read in the raw PRISMA file and preprocess
from PRISMA_preprocessing import preprocess_prisma
prisma_raw = pd.read_csv('%s/arid1a_interactions_raw.csv' % project_data_dir, header=1)
prisma = preprocess_prisma(prisma_raw)
prisma.to_csv('%s/arid1a_interactions_processed.csv' % project_data_dir, index = False)
prisma

### Perform grouping - Uncomment method required
from grouping_methods import get_combinatorial_grouping, get_hierachical_grouping, get_nmf_grouping
grouping = get_combinatorial_grouping(prisma)
# grouping = get_hierachical_grouping(prisma)
# grouping = get_nmf_grouping(prisma)
grouping

### Domain analysis
from domain_analysis_pipeline import preprocess_interpro, domain_analysis, domain_enrichment
interpro_raw = pd.read_csv('%s/protein2ipr_human.dat' % global_data_dir, delimiter= '\t')
id_map = pd.read_csv('%s/uniprot_map.csv' % global_data_dir)
interpro_anno = preprocess_interpro(interpro_raw, id_map)

domains_results = domain_analysis(grouping.copy(), interpro_anno)
protein_domains = domains_results[0]
group_domains = domains_results[1]
group_domains

domain_enrichment = domain_enrichment(grouping.copy(), interpro_anno)
domain_enrichment[0]
# domain_enrichment[0].to_csv('~/work/data/prisma_analysis/arid1a/domain_enrichment.csv', index=False)
# domain_enrichment[1].to_csv('~/work/data/prisma_analysis/arid1a/grouped_domain_enrichment.csv', index=False)

from validation_pipeline import validate_groups
corum = pd.read_csv('%s/corum_pairs_complex.csv' % global_data_dir)  ### need to preprocess raw corum
corr = pd.read_csv('%s/protein_corr_pearson.csv' % global_data_dir)
string = pd.read_csv('%s/9606.protein.links.detailed.v12.0.txt' % global_data_dir, sep = ' ')
string_info = pd.read_csv('%s/9606.protein.info.v12.0.txt' % global_data_dir, sep = '\t')

validation = validate_groups(grouping.copy(), corum, corr, string, string_info)
protein_pair_validation = validation[0]
protein_pair_validation


