protein_domain_counts = protein_domains.merge(group_domains[['group', 'domain', 'domain_protein_in_group_count']], on=['group', 'domain'], how='left').dropna()
protein_domain_counts.sort_values('domain_protein_in_group_count')
protein_domain_counts
shared_domain_count = protein_domain_counts.groupby(['group', 'proteins'])['domain_protein_in_group_count'].aggregate(sum)
shared_domain_count = pd.DataFrame(shared_domain_count).reset_index()
shared_domain_count.sort_values(['group', 'domain_protein_in_group_count'], ascending=False)
shared_domain_count

protein_singlet_validation = pd.concat([protein_pair_validation.drop(columns=['target']).rename(columns={'source':'proteins'}), protein_pair_validation.drop(columns=['source']).rename(columns={'target':'proteins'})])
protein_singlet_validation
interaction_count = pd.DataFrame(protein_singlet_validation.groupby(['group', 'proteins'])['overall_score'].aggregate(sum)).reset_index()
interaction_count

group_protein_scores = shared_domain_count.merge(interaction_count, on=['group', 'proteins'])
group_protein_scores.columns = ['group', 'protein', 'domains_shared', 'known_interaction_score']


group_protein_scores = group_protein_scores.sort_values(['group', 'domains_shared'])

group_protein_scores