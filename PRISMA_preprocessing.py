
import os
import pandas as pd
import re

def preprocess_prisma(df_prisma):

    # Rename the column containing gene names
    def rename_gene_column(df):
        """
        Renames columns in a pandas DataFrame that contain the substring 'Gene' at any part of its name into 'gene'.
        Raises a warning if there is more than one column that fits the criteria.
        """
        gene_column = [col for col in df.columns if 'gene' in col.lower()]
        if len(gene_column) > 1:
            print(f"Warning: More than one column contains the substring 'Gene': {gene_column}")
        if gene_column:
            df = df.rename(columns = {gene_column[0]: 'gene'})
        return df
    
    def extract_site(match):
        result = re.search(f'{site_pattern}', match)
        if result:
            site = result.group(0)
            return site
        else:
            return match

    df_prisma = rename_gene_column(df_prisma)

    # Rename the columns containing the interaction signal intensities
    site_pattern = r'P\d+'

    df_prisma.rename(columns = {col: extract_site(col) for col in df_prisma.columns}, inplace = True)
    df_prisma

    # Parsing the dataframe to a subset with gene names and intensity per site
    matrix_prisma = df_prisma.filter(regex = site_pattern, axis = 1)
    matrix_prisma = pd.concat([df_prisma['gene'], matrix_prisma], axis = 1)

    # Drop the NaN gene names
    matrix_prisma.dropna(inplace = True)
    # matrix_prisma = matrix_prisma.dropna(subset = 'gene')

    # Make sure all gene names are in uppercase
    matrix_prisma['gene'] = matrix_prisma['gene'].apply(lambda x: x.upper())

    # Unlist the gene names
    # Choose to only keep the first gene when unlisting homologs to aviod inflating
    matrix_prisma['gene'] = matrix_prisma['gene'].str.split(',|;').str[0]
    # df_prisma = df_prisma.explode('gene')

    # Remove duplicated gene names
    matrix_prisma.drop_duplicates(subset = 'gene', keep = 'first', inplace = True)

    # Return processed PRISMA matrix for downstream analyses
    return matrix_prisma


