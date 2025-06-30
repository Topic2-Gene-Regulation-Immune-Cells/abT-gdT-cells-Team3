def set_user(user):
    global who_this
    who_this = user

# import of packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats
import sklearn.manifold as sklm
import os
from sklearn.decomposition import PCA

# data_clean_up: vorläufig
def call_data_clean(p_threshold=None, qc_thresholds=None, normalization=None):
    ''' 
    loads data set;
    p_threshold: 0 < val > 5;
    qc_thresholds: input of qc_threshold_metrics;
    normalization: "none" or log2 automatically
    '''
    
    if who_this == "Laila":
        ATAC_seq_path = "/Users/laila/github/bioinfo25/Raw_data/ATAC-seq_called-peaks_ImmGenATAC18_AllOCRsInfo.csv"
        RNA_seq_path = "/Users/laila/github/bioinfo25/Raw_data/RNA-seq_mmc2.csv"
        Transcription_exons_path = "/Users/laila/github/bioinfo25/Raw_data/Transkrips-exon_refFlat.txt"
        Cell_population_qc_path = "/Users/laila/github/bioinfo25/Raw_data/cell-populations_qc-matrices_mmc1.xlsx"
        Voluntary_path = "/Users/laila/github/bioinfo25/Raw_data/voluntary_ImmGenATAC18_AllTFmotifsInOCRs.txt"
    elif who_this == "Kaja":
        ATAC_seq_path = "D:\\Uni\\4. Semester\\Bioinfo\\datasets\\processed atac seq data and called peaks.csv"
        RNA_seq_path = "D:\\Uni\\4. Semester\\Bioinfo\\datasets\\processed rna seq data.csv"
        Transcription_exons_path = "D:\\Uni\\4. Semester\\Bioinfo\\datasets\\refFlat.txt"
        Cell_population_qc_path = "D:\\Uni\\4. Semester\\Bioinfo\\datasets\\summary of immune cell populationsprofiled by atac.xlsx"
        Voluntary_path = "D:\\Uni\\4. Semester\\Bioinfo\\datasets\\summary of immune cell populationsprofiled by atac.xlsx"
    elif who_this == "Pia":
        ATAC_seq_path = "/Users/piakentschke/Documents/Uni/Data Analysis/ImmGenATAC18_AllOCRsInfo.csv"
        RNA_seq_path = "/Users/piakentschke/Documents/Uni/Data Analysis/RNA-seq_data.csv"
        Transcription_exons_path = "/Users/piakentschke/Documents/Uni/Data Analysis/refFlat.txt"
        Cell_population_qc_path = "/Users/piakentschke/Documents/Uni/Data Analysis/mmc1.xlsx"
        Voluntary_path = "/Users/piakentschke/Documents/Uni/Data Analysis/ImmGenATAC18_AllTFmotifsInOCRs.txt"
    elif who_this == "Helen":
        ATAC_seq_path = r"C:\Users\helen\Documents\downloads\ImmGenATAC18_AllOCRsInfo.csv"
        RNA_seq_path = r"C:\Users\helen\Documents\downloads\mmc2.csv" 
        Transcription_exons_path= r"C:\Users\helen\Downloads\datasets\Transkrips-exon_refFlat.txt"
        Cell_population_qc_path = r"C:\Users\helen\Dowloads\datasets\mmc1.xlsx"
        Voluntary_path = r"C:\Users\helen\Downloads\datasets\ImmGenATAC18_AllTFmotifsInOCRs.txt"

    # loading of data sets
    ATAC_seq = pd.read_csv(ATAC_seq_path, keep_default_na=False, header=0, index_col=0)
    RNA_seq = pd.read_csv(RNA_seq_path, header=0, index_col=0)
    QC_metrics = pd.read_excel(Cell_population_qc_path, header=0, index_col=0)
    exons=pd.read_csv(Transcription_exons_path, sep="\t")
    exons.columns = ["Genname", "ID", "Chromosom", "Strand", "Transcription Start", "Transcription End", "CDS_Start", "CDS_End", 
    "ExonCount", "ExonStarts", "ExonEnds"]

    # remove NAs and inf values
    ## ATAC-seq
    inf_rows = ATAC_seq[ATAC_seq['_-log10_bestPvalue'].isin([np.inf, -np.inf])]
    ATAC_seq = ATAC_seq[ATAC_seq['_-log10_bestPvalue'] != np.inf]
    ## RNA-seq
    RNA_NA = RNA_seq.isna().sum().sum()
    RNA_val = (RNA_seq.iloc[:, :] < 5).sum().sum()
    RNA_seq = RNA_seq[(RNA_seq >= 5).all(axis=1)]
    RNA_seq = RNA_seq.dropna(axis=0, how='any')
    QC_metrics = QC_metrics.dropna(axis=0, how='any')
    ATAC_seq = ATAC_seq.dropna(axis=0, how='any')
    
    # lisiting for sub sets and transponing of ATAC
    list_ATAC_stem_Tc_Bc = list(ATAC_seq.loc[:,'LTHSC.34-.BM': 'MPP4.135+.BM'])
    list_ATAC_diff_Tc_all= list(ATAC_seq.loc[:,'preT.DN1.Th':'Tgd.Sp'])
    list_ATAC_diff_Tc_pre_ab_act = list(ATAC_seq.loc[:,'preT.DN1.Th':'NKT.Sp.LPS.3d'])
    list_ATAC_diff_Tc_gd = list(ATAC_seq.loc[:,'Tgd.g2+d17.24a+.Th':'Tgd.Sp'])
    list_ATAC_Tc_all = list_ATAC_stem_Tc_Bc + list_ATAC_diff_Tc_all
    list_ATAC_Tc_all = list(dict.fromkeys(list_ATAC_Tc_all))

    ATAC_seq_T = ATAC_seq.T
    ATAC_seq_only_scores = ATAC_seq.loc[:,'LTHSC.34-.BM':]

    # normalization
    if normalization is None:
        # log2 normalization
        ATAC_seq_only_scores_norm = np.log2(ATAC_seq_only_scores)
        ATAC_seq_only_head = ATAC_seq.loc[:,:'genes.within.100Kb']
        ATAC_seq = pd.concat([ATAC_seq_only_head, ATAC_seq_only_scores_norm], axis=1)
        ATAC_seq_T = ATAC_seq.T

        RNA_seq = np.log2(RNA_seq)
        RNA_seq_T = RNA_seq.T

    elif normalization == "none":
        # no normalization
        pass

    #thresholds
    if p_threshold is not None:
        ATAC_seq = ATAC_seq[ATAC_seq["_-log10_bestPvalue"] >= p_threshold]
    
    #summary

    data = {
        'ATAC_seq': ATAC_seq,
        'ATAC_seq_T': ATAC_seq_T,
        'ATAC_seq_only_scores': ATAC_seq_only_scores,
        'norm': ATAC_seq_only_scores_norm,
        'RNA_seq': RNA_seq,
        'QC_metrics': QC_metrics,
        'exons': exons,
        'list_ATAC_stem_Tc_Bc': list_ATAC_stem_Tc_Bc,
        'list_ATAC_diff_Tc_all': list_ATAC_diff_Tc_all,
        'list_ATAC_diff_Tc_pre_ab_act': list_ATAC_diff_Tc_pre_ab_act,
        'list_ATAC_diff_Tc_gd': list_ATAC_diff_Tc_gd,
        'list_ATAC_Tc_all': list_ATAC_Tc_all
    }
    return data

# cluster

# UMAP

# t-SNE
def tSNE (df, cols, components, perplexity, rows=None, gini_coloring=None):
    '''
    df: pandas.df [rows x cols];
    cells: list of cell row names;
    components: no. of PCAs < len(rows);
    perplexity: perplexity < len(rows);
    gini_coloring: feature that should be colored, automatically maximun, manually after e.g. specifig cell / TSS
    '''
    
    # preparation of data
    if rows is not None:
        subset_df = df.loc[rows, cols]
    else:
        subset_df = df.loc[:,cols]
    components = min(components, subset_df.shape[0], subset_df.shape[1])
    perplexity = min(perplexity, subset_df.shape[0] - 1)

    # PCA
    pca = PCA(n_components= components)
    pcs = pca.fit_transform(subset_df)

    # t-SNE
    tsne = sklm.TSNE(n_components=2, perplexity=perplexity, random_state=42)
    tsne_results = tsne.fit_transform(pcs)
    tsne_df = pd.DataFrame(tsne_results, columns=['tSNE1', 'tSNE2'], index=subset_df.index)

    # gini index
    def gini_index (x):
        x = np.array(x)
        x = np.sort(x)
        n = len(x)
        if np.mean(x) == 0:
            return 0.0
        diff_sum = np.abs(np.subtract.outer(x, x)).sum()
        return diff_sum / (2 * n**2 * np.mean(x))
    
    gini_scores = subset_df.apply(gini_index, axis=0)
    
    # coloring
    if gini_coloring is None:
        name_gini_coloring = gini_scores.idxmax()
    elif isinstance(gini_coloring, int):
        name_gini_coloring = subset_df.columns[gini_coloring]
    else:
        name_gini_coloring = gini_coloring

    color_values = subset_df[name_gini_coloring].values

    # plot
    plt.figure(figsize=(8,6))
    sc = plt.scatter(tsne_df['tSNE1'], tsne_df['tSNE2'], c=color_values, cmap='viridis', alpha=0.7)
    plt.colorbar(sc, label=f'activity of {name_gini_coloring}')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title(f't-SNE, colored with gini index after peak: {name_gini_coloring}')
    plt.tight_layout()
    plt.show()

    return tsne_df, gini_scores

# correlation_pearson

# linear regression