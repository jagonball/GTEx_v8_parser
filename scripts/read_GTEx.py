'''
Read and search GTEx Analysis V8 file
'''
import os
from types import new_class
import pandas as pd
import re

__author__ = "Johnathan Lin <jagonball@gmail.com>"
__email__ = "jagonball@gmail.com"

def main():
    ### The target gene names
    genes = ['CCDC158', 'MKL1', 'DLG1', 'BCAS3', 'TBX2']
    # Create a DataFrame to store matched rows
    df_new = pd.DataFrame()
    ### Save GTEx Analysis file names into a list
    file_list = os.listdir('.//data//gene_tpm//')
    #print(file_list)
    ### Whether to remove data or not
    remove_data = False
    # Loop over each file
    for file in file_list:
        print(f'Current file: {file}')
        file_path = f'.//data//gene_tpm//{file}'
        df_new = get_gene_match(file_path, genes, df_new, remove_data)
        print(f'{df_new.shape} df_new')
    # Output df_new to file
    df_new.to_csv('./analysis/df_new.txt', sep='\t', index=False)

#    ### Separate each gene into a DataFrame
#    for i in genes:
#        print(i)
#        df_sep = gene_separator(df_new, i)
#        print(df_sep)
    

def file_importer(file_path):
    df = pd.read_csv(file_path, sep='\t', header=0, skiprows=2)
    return df


def gene_search(df, genes, source, remove_data = False):
    df_match = df.query(f'Description == {genes}')
    # Add source and average column
    df_match.insert(loc = 3, column = 'Source', value = source)
    df_match.insert(loc = 4, column = 'Average', value = df_match.iloc[:, 5:].mean(axis=1))
    if remove_data == True:
        df_match = df_match.iloc[:, 0:5]
    return df_match


def get_gene_match(file_path, genes, df_new, remove_data):
    ### Read GTEx file into DataFrame df_raw
    df_raw = file_importer(file_path)
    #print(df_raw.head())
    print(f'{df_raw.shape} df_raw')
    ### Rename columns from sample ID to subject ID
    # Slice the column names with regular expression
    reg = r'^([\w]+-[\w]+)'
    new_columns = [re.match(reg, x).group() for x in df_raw.columns[3:]]
    colnames = ['id', 'Name', 'Description']
    colnames += new_columns # Full new column names
    df_raw.columns = colnames
    ### Search DataFrame with list of gene names
    ### Return a new DataFrame with selected genes
    # Retrieve source name, turn list into string
    source_name = re.findall('v8_(.*).gct', file_path)\
                                    [0].replace('_', ' ')
    df_match = gene_search(df_raw, genes, source_name,
                           remove_data=remove_data)
    print(f'{df_match.shape} df_match')
    # append to df_new
    df_new = pd.concat([df_new, df_match], axis=0, join='outer')
    return df_new


#def gene_separator(df, gene_name):
#    df_sep = df[df['Description'] == gene_name]
#    return df_sep

if __name__ == '__main__':
    main()
    