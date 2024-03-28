
import pandas as pd

prisma_protein_label  = 'ARID1A'
grouping_method = 'combinatorial'

## Define relative and absolute paths to read and write files
project_data_dir = './project_data/'
global_data_dir = './global_data/'
                                                     
## Read in the raw PRISMA file and preprocess
from PRISMA_preprocessing import preprocess_prisma
print('Loading and processing PRISMA data...')
prisma_raw = pd.read_csv('%s/arid1a_interactions_raw.csv' % project_data_dir, header=1)
prisma = preprocess_prisma(prisma_raw)
prisma.to_csv('%s/arid1a_interactions_processed.csv' % project_data_dir, index = False)
prisma

### Perform grouping - Uncomment method required
print('Grouping data...')
from grouping_methods import get_combinatorial_grouping, get_hierachical_grouping, get_nmf_grouping
grouping = get_combinatorial_grouping(prisma)
# grouping = get_hierachical_grouping(prisma)
# grouping = get_nmf_grouping(prisma)
grouping

### Domain analysis
from domain_analysis_pipeline import preprocess_interpro, domain_analysis, domain_enrichment
## Human protein / domain pairs extracted from protein2ipr.dat - download from https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz
interpro_raw = pd.read_csv('%s/protein2ipr_human.dat' % global_data_dir, delimiter= '\t')
## ID map downloaded from ensembl.org/biomart extracting: Gene name, UniProtKB/Swiss-Prot ID
id_map = pd.read_csv('%s/uniprot_map.csv' % global_data_dir)
interpro_anno = preprocess_interpro(interpro_raw, id_map)

print('Performing group domain analysis...')
domains_results = domain_analysis(grouping.copy(), interpro_anno)
protein_domains = domains_results[0]
group_domains = domains_results[1]
domain_summary = domains_results[2]
protein_domains.to_csv('%s/%s_%s_protein_domains.csv' % (project_data_dir, prisma_protein_label, grouping_method), index=False)
group_domains.to_csv('%s/%s_%s_group_domains.csv' % (project_data_dir, prisma_protein_label, grouping_method), index=False)
domain_summary.to_csv('%s/%s_%s_domain_summary.csv' % (project_data_dir, prisma_protein_label, grouping_method), index=False)

print('Performing group domain enrichment analysis...')
domain_enrichment_results = domain_enrichment(grouping.copy(), interpro_anno)
protein_domain_enrichment = domain_enrichment_results[0]
grouped_domain_enrichment = domain_enrichment_results[1]
protein_domain_enrichment.to_csv('%s/%s_%s_protein_domain_enrichment.csv' % (project_data_dir, prisma_protein_label, grouping_method), index=False)
grouped_domain_enrichment.to_csv('%s/%s_%s_grouped_domain_enrichment.csv' % (project_data_dir, prisma_protein_label, grouping_method), index=False)

### Group protein pair interaction validation
print('Loading external interaction data for validation...')
from validation_pipeline import validate_groups
## Load external interaction data
## Corum data downloaded from https://mips.helmholtz-muenchen.de/corum/ and preprocessed 
corum = pd.read_csv('%s/corum_pairs_complex.csv' % global_data_dir)  ### need to preprocess raw corum
## Pearson correlation of normalised CCLE protein abundance data https://gygi.hms.harvard.edu/publications/ccle.html
corr = pd.read_csv('%s/protein_corr_pearson.csv' % global_data_dir)
## STRING data downlaoded from https://stringdb-downloads.org/download/protein.links.detailed.v12.0.txt.gz
string = pd.read_csv('%s/9606.protein.links.detailed.v12.0.txt' % global_data_dir, sep = ' ')
## STRING annotation - https://stringdb-downloads.org/download/protein.info.v12.0.txt.gz
string_info = pd.read_csv('%s/9606.protein.info.v12.0.txt' % global_data_dir, sep = '\t')

print('Performing group protein pair interaction validation...')
validation = validate_groups(grouping.copy(), corum, corr, string, string_info)
protein_pair_validation = validation[0]
group_average_validation = validation[1]

protein_pair_validation.to_csv('%s/%s_%s_protein_pair_validation.csv' % (project_data_dir, prisma_protein_label, grouping_method), index=False)
group_average_validation.to_csv('%s/%s_%s_group_average_validation.csv' % (project_data_dir, prisma_protein_label, grouping_method), index=False)

