'''
Generate plot
'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    ### Read file
    file_path = './/analysis//df_new.txt'
    df = pd.read_csv(file_path, sep='\t', header=0)
    #print(df.head())
    ### The target gene names
    genes = ['CCDC158', 'BCAS3', 'MKL1', 'DLG1']
    gene = 'CCDC158'

    # Set plot style
    plt.style.use('ggplot')
    #sns.set_theme(style="ticks")
    # Bar chart
#    bar_chart(df_plot, gene)

#    ### Generate boxplot for each gene in genes
#    for gene in genes:
#        # Subset DataFrame to gene for plot drawing
#        df_plot = tissue_filter(df, gene)
#        # Transpose df for boxplot input
#        df_plot = df_plot.set_index('Source')
#        df_plot = df_plot.iloc[:, 4:].T
#        #print(df_plot.head())
#        df_plot = df_plot.drop(columns='testis')
#        boxplot(df_plot, gene)

    ### Filter df with gene and tissue
    df_plot = tissue_filter(df, gene, tissue='kidney cortex')
    # Transpose df for boxplot input
    df_plot = df_plot.set_index('Source')
    df_plot = df_plot.iloc[:, 4:]
    print(df_plot.head())


def boxplot(df, gene):
    fig, ax1 = plt.subplots(figsize=(10, 4.5))
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
    plt.xticks(fontsize=8, rotation=45, ha='right', rotation_mode='anchor')
#    num_boxes = len(df.columns)
#    ax1.set_xlim(0.5, num_boxes + 0.5)
    #top = 40
    #bottom = -5
    #ax1.set_ylim(bottom, top)
    #ax1.set_xticks(range(num_boxes), labels=df.columns, rotation=45, fontsize=8)
    #ax1.set_xticklabels(df.columns,
    #                    rotation=45, fontsize=8)
    plt.tight_layout() ### Rescale the fig size to fit the data
    plt.savefig(f'./analysis/boxplot_{gene}.png')
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


def tissue_filter(df, gene, tissue='all'):
    ### Filter df with gene and tissue
    # Subset DataFrame to gene for plot drawing
    df_plot = df[df['Description'] == gene]
    if tissue != 'all':
        #print(f'{df_plot.shape} df_plot')
        # Select the tissue row and check which columns have value
        col_bool = df_plot[df_plot['Source'] == tissue].notna().iloc[0,:]
        # Use the above boolean to get the column names to keep
        col_to_keep = df_plot.columns[col_bool]
        df_plot = df_plot.loc[:, col_to_keep]
        #print(f'{df_plot.shape} df_plot')
        # Add new column 'boxplot' for different category
        df_plot.insert(loc = 5, column = 'boxplot', value = 'others')
        tissue_index = df_plot.index[df_plot['Source'] == tissue].tolist()
        # Replace the target tissue value
        df_plot.at[tissue_index[0], 'boxplot'] = tissue
    else:
        pass
    df_plot.to_csv(f'.//analysis//plot_data//df_plot_{gene}_{tissue}.txt', sep='\t', index=False)
    return df_plot



if __name__ == '__main__':
    main()