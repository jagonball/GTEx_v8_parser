'''
Perform data filter, analysis then generate plot. 
'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_ind
import os
import re
#from statannot import add_stat_annotation

def main():
    ### Read file
    file_path = './/analysis//df_new.txt'
    df = pd.read_csv(file_path, sep='\t', header=0)
    #print(df.head())
    ### The target gene names
    genes = ['BCAS3', 'TBX2', 'CCDC158', 'MKL1', 'DLG1']
    remove_0 = True
    # Set plot style
    plt.style.use('ggplot')
    #sns.set_theme(style="ticks")
    ### Bar chart ###
#    bar_chart(df_plot, gene)
    

#    ### Generate boxplot for each gene in genes ###
#    ### compared to all other tissues
#    for gene in genes:
#        # Subset DataFrame to gene for plot drawing
#        df_plot = tissue_filter(df, gene, remove_0=remove_0)       
#        # Transpose df for boxplot input
#        df_plot = df_plot.set_index('Source')
#        df_plot = df_plot.iloc[:, 4:].T
#        print(df_plot.head())
#        # Remove outlier data
#        if gene == 'CCDC158':
#            df_plot = df_plot.drop(columns='testis')
#        boxplot(df_plot, gene, remove_0=remove_0)


#    ### Perform target tissue versus other tissue and generate boxplot ###
#    ### Now with 2 target tissues
#    # Filter df with gene and tissue
#    for gene in genes:
#        tissue = 'kidney cortex'
#        tissue_2 = 'kidney medulla'
#        other_tissue = 'whole blood'
#        ### DataFrame for gene
#        # Subset DataFrame to gene
#        df_plot = df[df['Description'] == gene]
#        # Remove 0
#        if remove_0:
#            df_plot = df_plot.replace(0, np.nan)
#        # Remove outlier data
#        if gene == 'CCDC158':
#            df_plot = df_plot[df_plot['Source'] != 'testis']
#        # Save to file
#        df_plot.to_csv(f'.//analysis//filtered_data//df_plot_{gene}_{tissue}.txt',
#                        sep='\t', index=False)
#        # Tissue and tissue_2 rows into arrays for boxplot
#        array_1, array_other, array_2 = tissue_arrays_for_ttest(df_plot, tissue,
#                                                   tissue_2=tissue_2,
#                                                   other=other_tissue,
#                                                   dropna=True)
#        # Print array size
#        print(f'array_1: {array_1.shape}')
#        print(f'array_other: {array_other.shape}')
#        print(f'array_2: {array_2.shape}')
#        # Perform t-test
#        t_test_result, ttest_type = t_test(array_1, array_other)
#        p_value = t_test_result[1]
#        p_value = '{:.3g}'.format(p_value)
#        print(p_value)
#        t_test_result_2, ttest_type = t_test(array_2, array_other)
#        p_value_2 = t_test_result_2[1]
#        p_value_2 = '{:.3g}'.format(p_value_2)
#        print(p_value_2)
#        # Create a list of arrays
#        arrays = [array_1, array_other, array_2]
#        # Custom boxplot color
#        #box_color = ['red', 'g', 'blue', 'g']
#        # Generate boxplot
#        boxplot(arrays, gene, tissue=tissue, other=other_tissue, tissue_2=tissue_2,
#                p_value=p_value, p_value_2=p_value_2, remove_0=remove_0)
#
#
#    ### Perform t-test between target tissue against all other tissues
#    tissue = 'kidney cortex'
#    print(f'Target tissue: {tissue}')
#    for gene in genes:
#        print(f'Current gene: {gene}')
#        # Create lists for new data
#        sample_size = []
#        test_type = []
#        statistic = []
#        p_value = []
#        # Subset DataFrame to gene
#        df_gene = df[df['Description'] == gene]
#        #print(f'{df_gene.shape} df_gene')
#        df_gene = df_gene.reset_index(drop=True)
#        ### Prepare data for t-test
#        # Target tissue array
#        array_1 = tissue_array(df_gene, tissue)
#        # Remove nan
#        array_1 = array_1[~np.isnan(array_1)]
#        #print(array_1)
#        # Get other tissue array then perform t-test
#        df_gene_other = df_gene[df_gene['Source'] != tissue]
#        for other_tissue in df_gene_other['Source']:
#            print(f'Current tissue: {other_tissue}')
#            array_2 = tissue_array(df_gene_other, other_tissue)
#            # Remove nan
#            array_2 = array_2[~np.isnan(array_2)]
#            # Print array size
#            print(f'Target tissue size: {array_1.shape}')
#            print(f' Other tissue size: {array_2.shape}')
#            # Add to list sample_size
#            sample_size.append(array_2.shape[0])
#            # Perform t-test
#            ttest_result, ttest_type = t_test(array_1, array_2)
#            # Add results to lists statistic and p_value
#            test_type.append(ttest_type)
#            statistic.append(ttest_result[0])
#            p_value.append(ttest_result[1])
#        # Subset df_gene_other for t-test output
#        df_ttest = df_gene_other.iloc[:, 0:4]
#        #print(len(sample_size))
#        #print(len(test_type))
#        #print(len(statistic))
#        #print(len(p_value))
#        df_ttest['Sample size'] = sample_size
#        df_ttest['Test type'] = test_type
#        df_ttest['Statistic'] = statistic
#        df_ttest['P value'] = p_value
#        # Add target tissue row into output
#        row_target_tissue = df_gene[df_gene['Source'] == tissue].iloc[:, 0:4]
#        #print(row_target_tissue)
#        row_target_tissue['Sample size'] = array_1.shape[0]
#        df_ttest = pd.concat([df_ttest, row_target_tissue])
#        df_ttest = df_ttest.sort_index()
#        df_ttest.to_csv(f'.//analysis//t_test/df_ttest_{gene}_{tissue}.txt',
#                        sep='\t', index=False)


    ### Heatmap ###
    # Create a DataFrame to store p-values
    df_pvalue = pd.DataFrame()
    # File list of t-test output
    file_list = os.listdir('.//analysis//t_test//')
    # Loop over each file
    for file in file_list:
        print(f'Current file: {file}')
        df_pvalue = combine_pvalue(file, df_pvalue)
        print(f'{df_pvalue.shape} df_pvalue')
    # Remove index name
    df_pvalue.index.name = None
    #print(df_pvalue.columns)
    # Rearrange column for heatmap
    df_pvalue = df_pvalue.iloc[:, [0, 8, 2, 6, 4, 1, 9, 3, 7, 5]]
    #print(df_pvalue.columns)
    # Save to file
    df_pvalue.to_csv(f'.//analysis//df_pvalue.txt',
                        sep='\t', index=True)
    heatmap(df_pvalue, annotate=True)


#    ### Histogram ###
#    file_path = './/data//gene_tpm//gene_tpm_2017-06-05_v8_kidney_cortex.gct'
#    df_hist = file_importer(file_path)
#    #print(df_hist.head())
#    # Remove 0
#    df_hist.iloc[:, 3:] = df_hist.iloc[:, 3:].replace(0, np.nan)
#    # Add median column
#    df_hist.insert(loc = 3, column = 'Median', value = df_hist.iloc[:, 3:].median(axis=1))
#    print(df_hist.shape)
#    # Remove rows with no value
#    df_hist = df_hist[df_hist['Median'].notna()]
#    df_hist.insert(loc=4, column = 'Median_log10', value = np.log10(df_hist['Median']))
#    print(df_hist.shape)
#    # Save to file
#    df_hist.to_csv(f'.//analysis//df_hist.txt',
#                        sep='\t', index=False)
#    # Retrieve target gene values and save into a dictionary
#    gene_value = {}
#    for gene in genes:
#        value = df_hist[df_hist['Description'] == gene]['Median_log10']
#        gene_value[gene] = value.values[0]
#    print(gene_value)
#    # Generate histogram for column
#    column = 'Median_log10'
#    histogram(df_hist, column)



def boxplot(df, gene, tissue='all', other='others', tissue_2=None,
            p_value=None, p_value_2=None, remove_0=False, palette=None):
    ### Rule for figsize ###
    if tissue == 'all':
        fig, ax1 = plt.subplots(figsize=(10, 5))
    else:
        fig, ax1 = plt.subplots(figsize=(6, 4))
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    sns.boxplot(data = df,
                width = 0.3,
                palette = palette,
                fliersize = 2,
                linewidth = 0.5,
                ax = ax1)
    ax1.set_title(gene,
                  fontsize = 16,
                  fontweight = 'bold')
    # Set the axes labels
    if tissue == 'all':
        plt.xticks(fontsize=8, rotation=45, ha='right', rotation_mode='anchor')
    else:
        plt.xticks(fontsize=16)
        x_labels = [tissue, other, tissue_2]
        ax1.set_xticklabels(x_labels)
    ax1.set_ylabel('gene TPM',
                   fontsize = 14)
                   #fontweight = 'bold')
    # Annotate p_value
    if p_value:
        bbox = dict(boxstyle ="round", fc ="0.8")
        ax1.annotate(text = f'      p value\nkc.wb. {p_value}\nkm.wb. {p_value_2}',
                     xy = (0.76, 0.82),
                     xycoords = 'axes fraction',
                     fontstyle = 'italic',
                     bbox = bbox)
    #    ax1.annotate(text = f'p value\n{p_value_2}',
    #                 xy = (0.60, -0.15),
    #                 xycoords = 'axes fraction',
    #                 fontstyle = 'italic',
    #                 alpha = 0.5)
    plt.tight_layout() ### Rescale the fig size to fit the data
    ### Rule for file name ###
    if tissue == 'all':
        if remove_0:
            plt.savefig(f'.//analysis//plot//boxplot_{gene}_{tissue}_rm0.png', dpi=900)
        else:
            plt.savefig(f'.//analysis//plot//boxplot_{gene}_{tissue}.png', dpi=900)
    elif tissue_2:
        if remove_0:
            plt.savefig(f'.//analysis//plot//boxplot_{gene}_{tissue}&{tissue_2}_{other}_rm0.png', dpi=900)
        else:
            plt.savefig(f'.//analysis//plot//boxplot_{gene}_{tissue}&{tissue_2}_{other}.png', dpi=900)
    else:
        if remove_0:
            plt.savefig(f'.//analysis//plot//boxplot_{gene}_{tissue}_{other}_rm0.png', dpi=900)
        else:
            plt.savefig(f'.//analysis//plot//boxplot_{gene}_{tissue}_{other}.png', dpi=900)
    plt.show()
    

def bar_chart(df, gene):
    # Create a figure
    fig, ax = plt.subplots()
    tissue = df['Source']
    counts = df['Average']    
    ax.bar(tissue, counts)
    ax.set_ylabel('Average')
    ax.set_title(f'{gene}')
    plt.show()


def tissue_filter(df, gene, tissue='all',
                  savefile=True, remove_0=False):
    """Filter DataFrame with gene, 
    and remove columns with no value on target tissue row.

    Args:
        df (DataFrame): Input DataFrame.
        gene (str): Target gene.
        tissue (str, optional): Target tissue. Defaults to 'all'.
        savefile (bool, optional): Save filtered DataFrame to file.\
                                   Defaults to True.
        remove_0 (bool, optional): Remove 0 from DataFrame.\
                                   Defaults to False.
    """
    ### Filter df with gene and tissue
    # Subset DataFrame to gene
    df_gene = df[df['Description'] == gene]
    if tissue != 'all':
        #print(f'{df_gene.shape} df_gene')
        # Select the tissue row and check which columns have value
        col_bool = df_gene[df_gene['Source'] == tissue].notna().iloc[0,:]
        # Use the above boolean to get the column names to keep
        col_to_keep = df_gene.columns[col_bool]
        df_gene = df_gene.loc[:, col_to_keep]
        #print(f'{df_gene.shape} df_gene')
#        # Add new column 'boxplot' for different category
#        df_gene.insert(loc = 5, column = 'boxplot', value = 'others')
#        # Get the target tissue index then replace the boxplot value
#        tissue_index = df_gene.index[df_gene['Source'] == tissue].tolist()
#        df_gene.at[tissue_index[0], 'boxplot'] = tissue
    if remove_0:
        df_gene = df_gene.replace(0, np.nan)
    if savefile:
        df_gene.to_csv(f'.//analysis//filtered_data//df_gene_{gene}_{tissue}.txt',
                        sep='\t', index=False)
    return df_gene


def tissue_array(df, tissue):
    # Filter df with tissue on column 'Source'
    df_temp = df[df['Source'] == tissue]
    #print(df_temp)
    # Flatten 2D array to 1D with ".ravel()"
    array = df_temp.iloc[:, 5:].to_numpy(copy=True).ravel()
    #print(array)
    return array


def tissue_arrays_for_ttest(df, tissue, tissue_2=None,
                            other='all', dropna=False):
    ## Target tissue ##
    array_1 = tissue_array(df, tissue)
    if tissue_2:
        array_2 = tissue_array(df, tissue_2)
    ## Other tissues ##
    if other == 'all':
        array_other = multiple_tissues_array(df, tissue)
    else:
        array_other = tissue_array(df, other)
    if dropna:
        # Remove nan
        array_1 = array_1[~np.isnan(array_1)]
        array_2 = array_2[~np.isnan(array_2)]
        array_other = array_other[~np.isnan(array_other)]
    return array_1, array_other, array_2


def multiple_tissues_array(df, tissue):
    df_temp = df[df['Source'] != tissue]
    #print(df_temp)
    df_temp = df_temp.iloc[:, 5:]
    # Create a list for other tissues, turn into an array later.
    array = []
    # Iterate over rows with Index and Series
    for i, s in df_temp.iterrows():
        # Save series to temporary list
        list_temp = s.dropna().to_list()
        #print(list_temp)
        array += list_temp
    array = np.array(array)
    #print(array)
    return array


def array_remove_0(array):
    # Remove "0" from array
    array = array[array != 0]
    return array


def t_test(array_1, array_2):
    # Check variance of both data groups
    var_1 = np.var(array_1)
    var_2 = np.var(array_2)
    print(f'array_1 variance: {var_1}\n'\
          f'array_2 variance: {var_2}')
    # Variance ratio
    var_ratio = var_1/var_2
    print(f'variance ratio: {var_ratio}')
    # If variance ratio is greater than 4, perform Welch’s t-test
    if var_ratio > 4 or var_ratio < 0.25:
        print('Welch’s t-test')
        ttest_type = "Welch's t-test"
        ttest_result = ttest_ind(array_1, array_2, equal_var=False)
    else:
        print('Two-sample t-test')
        ttest_type = 'Two-sample t-test'
        ttest_result = ttest_ind(array_1, array_2)
    print(ttest_result)
    return ttest_result, ttest_type


def file_importer(file_path):
    df = pd.read_csv(file_path, sep='\t', header=0, skiprows=2)
    return df


def combine_pvalue(file, df_pvalue):
    ### Read file into DataFrame df
    file_path = f'.//analysis//t_test//{file}'
    df = pd.read_csv(file_path, sep='\t', header=0)
    #print(df.head())
    # Set index to 'Source' then get Series 'P_value'
    df = df.set_index('Source')
    series = df['P value']
    # Rename series based on file name
    source_name = re.findall('ttest_(.*).txt', file)[0]
    series =series.rename(source_name)
    # Add series to df_pvalue as new column
    df_pvalue[source_name] = series
    return df_pvalue


def heatmap(df, annotate=False):
    fig, ax1 = plt.subplots()
    fig.set_figwidth(12)
    fig.set_figheight(16)
    sns.heatmap(df,
                vmin=0,
                vmax=0.1,
                cmap='coolwarm',
                center=0.05,
                cbar=True,
                cbar_kws={'shrink':0.4,
                          'label':'p value'},
                linewidths=0.5,
                linecolor='grey',
                xticklabels=True,
                yticklabels=True,
                ax=ax1
                )
    plt.xticks(fontsize=12, rotation=45, ha="right",
               rotation_mode="anchor")
    plt.yticks(fontsize=10)
    # Annotate p-values
    if annotate:
        # Loop over data dimensions and create text annotations.
        for y in range(df.shape[0]):
            for x in range(df.shape[1]):
                plt.text(x + 0.5, y + 0.5, '{:.3g}'.format(df.iloc[y, x]),
                         horizontalalignment='center',
                         verticalalignment='center',
                         color="#afaf00")
        fig.tight_layout() ### Rescale the fig size to fit the data
        plt.savefig(f'.//analysis//plot//heatmap_annotate.png', dpi=900)
    else:
        fig.tight_layout() ### Rescale the fig size to fit the data
        plt.savefig(f'.//analysis//plot//heatmap.png', dpi=900)
    plt.show()


def histogram(df, x):
    fig, ax1 = plt.subplots()
    sns.histplot(df, x=x, ax=ax1)
    plt.tight_layout() ### Rescale the fig size to fit the data
    plt.show()


if __name__ == '__main__':
    main()