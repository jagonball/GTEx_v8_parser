'''
Generate plot
'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_ind

def main():
    ### Read file
    file_path = './/analysis//df_new.txt'
    df = pd.read_csv(file_path, sep='\t', header=0)
    #print(df.head())
    ### The target gene names
    genes = ['CCDC158', 'BCAS3', 'MKL1', 'DLG1']
    # Set plot style
    plt.style.use('ggplot')
    #sns.set_theme(style="ticks")
    # Bar chart
#    bar_chart(df_plot, gene)

    remove_0 = False
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

#    ### Filter df with gene and tissue ###
#    for gene in genes:
#        tissue = 'kidney medulla'
#        other_tissue = 'whole blood'
#        df_plot = tissue_filter(df, gene, tissue=tissue, 
#                                savefile=True, remove_0=remove_0)
#        # Remove outlier data
#        if gene == 'CCDC158':
#            df_plot = df_plot[df_plot['Source'] != 'testis']
#        # Turn rows into arrays for boxplot
#        array_1, array_2 = tissue_arrays_for_ttest(df_plot, tissue)#, other=other_tissue)
#        # Print array size
#        print(f'array_1: {array_1.shape}')
#        print(f'array_2: {array_2.shape}')
#        # Perform t-test
#        t_test_result = t_test(array_1, array_2)
#        p_value = t_test_result[1]
#        p_value = '{:.3g}'.format(p_value)
#        print(p_value)
#        # Create a list of arrays
#        arrays = [array_1, array_2]
#        # Generate boxplot
#        boxplot(arrays, gene, tissue=tissue, #other=other_tissue,
#                p_value=p_value, remove_0=remove_0)


    ### Perform t-test between target tissue against all other tissues
    for gene in genes:
        tissue = 'kidney cortex'
        # Create lists for new data
        sample_size = []
        statistic = []
        p_value = []
        # Subset DataFrame to gene
        df_gene = df[df['Description'] == gene]
        #print(f'{df_gene.shape} df_gene')
        # Target tissue array
        array_1 = tissue_array(df_gene, tissue)
        # Remove nan
        array_1 = array_1[~np.isnan(array_1)]
        print(array_1)
        # Get other tissue array then perform t-test
        df_gene_other = df_gene[df_gene['Source'] != tissue]
        for other_tissue in df_gene_other['Source']:
            array_2 = tissue_array(df_gene_other, other_tissue)
            # Remove nan
            array_2 = array_2[~np.isnan(array_2)]
            # Print array size
            print(f'array_1: {array_1.shape}')
            print(f'array_2: {array_2.shape}')
            # Add to list sample_size
            sample_size.append(array_2.shape[0])
            # Perform t-test
            t_test_result = t_test(array_1, array_2)
            # Add results to lists statistic and p_value
            statistic.append(t_test_result[0])
            p_value.append(t_test_result[1])
        # Subset df_gee_other for t-test output
        df_ttest = df_gene_other.iloc[:, 0:4]
        #print(len(sample_size))
        #print(len(statistic))
        #print(len(p_value))
        df_ttest['Sample_size'] = sample_size
        df_ttest['Statistic'] = statistic
        df_ttest['P_value'] = p_value
        df_ttest.to_csv(f'.//analysis//t_test/df_ttest_{gene}_{tissue}.txt',
                        sep='\t', index=False)
        
        
#    ### Histogram ###
#    file_path = './/data//gene_tpm//gene_tpm_2017-06-05_v8_kidney_medulla.gct'
#    df = file_importer(file_path)
#    #print(df.head())
#    fig, ax1 = plt.subplots()
#    df_rm0 = df.iloc[:, 3:].replace(0, np.nan)
#    df_log = np.log10(df_rm0)
#    #print(df_log.head())
#    # Generate histogram for each column
#    for column in df_log.columns:
#        histogram(df_log, column)



def boxplot(df, gene, tissue='all', other='others',
            p_value=None, remove_0=False):
    ### Rule for figsize ###
    if tissue == 'all':
        fig, ax1 = plt.subplots(figsize=(10, 5))
    else:
        fig, ax1 = plt.subplots(figsize=(6, 5))
    #fig.canvas.manager.set_window_title('A Boxplot Example')
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    sns.boxplot(data=df,
                width=0.8,
                fliersize=1,
                linewidth=0.5,
                ax=ax1)

    ## Add in points to show each observation
    #sns.stripplot(data=df,
    #              size=1,
    #              color=".3",
    #              linewidth=0)

    #bp= df.boxplot(ax=ax1,
    #               fontsize=8,
    #               rot=45)
    #bp = ax1.boxplot(df, notch=False, sym='+', vert=True, whis=1.5)
    #plt.setp(bp['boxes'], color='black')
    #plt.setp(bp['whiskers'], color='black')
    #plt.setp(bp['fliers'], color='red', marker='+')

    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
#    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#                   alpha=0.5)
    ax1.set(#axisbelow=True,  # Hide the grid behind plot objects
            title=f'{gene}',
            xlabel=None,
            ylabel='gene TPM')
    
    # Set the axes labels
    if tissue == 'all':
        plt.xticks(fontsize=8, rotation=45, ha='right', rotation_mode='anchor')
    else:
        plt.xticks(fontsize=12)
        x_labels = [tissue, other]
        ax1.set_xticklabels(x_labels)
    # Annotate p_value
    if p_value:
        ax1.annotate(text = f'p value: {p_value}',
                     xy = (0.75, 1.02),
                     xycoords = 'axes fraction',
                     fontstyle = 'italic')
#    num_boxes = len(df.columns)
#    ax1.set_xlim(0.5, num_boxes + 0.5)
    #top = 40
    #bottom = -5
    #ax1.set_ylim(bottom, top)
    #ax1.set_xticks(range(num_boxes), labels=df.columns, rotation=45, fontsize=8)
    #ax1.set_xticklabels(df.columns,
    #                    rotation=45, fontsize=8)
    plt.tight_layout() ### Rescale the fig size to fit the data
    ### Rule for file name ###
    if tissue == 'all':
        if remove_0:
            plt.savefig(f'.//analysis//plot//boxplot_{gene}_{tissue}_rm0.png')
        else:
            plt.savefig(f'.//analysis//plot//boxplot_{gene}_{tissue}.png')
    else:
        if remove_0:
            plt.savefig(f'.//analysis//plot//boxplot_{gene}_{tissue}_{other}_rm0.png')
        else:
            plt.savefig(f'.//analysis//plot//boxplot_{gene}_{tissue}_{other}.png')
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


def tissue_arrays_for_ttest(df, tissue, other='all'):
    ## Target tissue ##
    array_1 = tissue_array(df, tissue)
    ## Other tissues ##
    if other == 'all':
        array_2 = multiple_tissues_array(df, tissue)
    else:
        array_2 = tissue_array(df, other)
    return array_1, array_2


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
        t_test_result = ttest_ind(array_1, array_2, equal_var=False)
    else:
        print('Two-sample t-test')
        t_test_result = ttest_ind(array_1, array_2)
    print(t_test_result)
    return t_test_result


def file_importer(file_path):
    df = pd.read_csv(file_path, sep='\t', header=0, skiprows=2)
    return df


def histogram(df, x):
    sns.histplot(df, x=x)
    plt.tight_layout() ### Rescale the fig size to fit the data
    plt.show()


if __name__ == '__main__':
    main()