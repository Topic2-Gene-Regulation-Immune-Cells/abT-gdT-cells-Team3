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
import scanpy as sc

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
        ATAC_seq_path = r"C:\Users\helen\Downloads\datasets\ImmGenATAC18_AllOCRsInfo.csv"
        RNA_seq_path = r"C:\Users\helen\Downloads\datasets\mmc2.csv" 
        Transcription_exons_path= r"C:\Users\helen\Downloads\datasets\Transkrips-exon_refFlat.txt"
        Cell_population_qc_path = r"C:\Users\helen\Downloads\datasets\mmc1.xlsx"
        Voluntary_path = r"C:\Users\helen\Downloads\datasets\ImmGenATAC18_AllTFmotifsInOCRs.txt"
        
    # loading of data sets
    ATAC_seq = pd.read_csv(ATAC_seq_path, header=0)
    ATAC_seq = ATAC_seq.set_index('ImmGenATAC1219.peakID')
    ATAC_copy = ATAC_seq.copy()
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

    
    # lisiting for sub sets and transponing of ATAC
    list_ATAC_stem_Tc_Bc = list(ATAC_seq.loc[:,'LTHSC.34-.BM': 'MPP4.135+.BM'])
    list_ATAC_diff_Tc_all= list(ATAC_seq.loc[:,'preT.DN1.Th':'Tgd.Sp'])
    list_ATAC_diff_Tc_pre_ab_act = list(ATAC_seq.loc[:,'preT.DN1.Th':'NKT.Sp.LPS.3d'])
    list_ATAC_diff_Tc_gd = list(ATAC_seq.loc[:,'Tgd.g2+d17.24a+.Th':'Tgd.Sp'])
    list_ATAC_diff_Tc_ab = (
    list(ATAC_seq.loc[:, 'preT.DN1.Th':'T.4.Nve.Fem.Sp']) +
    list(ATAC_seq.loc[:, 'Treg.4.FP3+.Nrplo.Co':'T.8.Nve.Sp']) +
    ['NKT.Sp']
)

    list_ATAC_diff_Tc_ab_gd = list_ATAC_diff_Tc_ab + list_ATAC_diff_Tc_gd
    list_ATAC_Tc_all = list_ATAC_stem_Tc_Bc + list_ATAC_diff_Tc_all
    list_ATAC_Tc_all = list(dict.fromkeys(list_ATAC_Tc_all))
    
    
    ab_tc= ['preT.DN1.Th','preT.DN2a.Th', 'preT.DN2b.Th', 'preT.DN3.Th', 'T.DN4.Th', 'T.ISP.Th', 'T.DP.Th', 'T.4.Th', 'T.8.Th', 'T.4.Nve.Sp', 'T.4.Nve.Fem.Sp', 'Treg.4.FP3+.Nrplo.Co', 'Treg.4.25hi.Sp', 'T.8.Nve.Sp', 'NKT.Sp']
    gd_tc=['Tgd.g2+d17.24a+.Th', 'Tgd.g2+d17.LN', 'Tgd.g2+d1.24a+.Th', 'Tgd.g2+d1.LN', 'Tgd.g1.1+d1.24a+.Th','Tgd.g1.1+d1.LN', 'Tgd.Sp']
    ab_gd_tc = ['preT.DN1.Th','preT.DN2a.Th', 'preT.DN2b.Th', 'preT.DN3.Th', 'T.DN4.Th', 'T.ISP.Th', 'T.DP.Th', 'T.4.Th', 'T.8.Th', 'T.4.Nve.Sp', 'T.4.Nve.Fem.Sp', 'Treg.4.FP3+.Nrplo.Co', 'Treg.4.25hi.Sp', 'T.8.Nve.Sp', 'NKT.Sp','Tgd.g2+d17.24a+.Th', 'Tgd.g2+d17.LN', 'Tgd.g2+d1.24a+.Th', 'Tgd.g2+d1.LN', 'Tgd.g1.1+d1.24a+.Th','Tgd.g1.1+d1.LN', 'Tgd.Sp']

    #thresholds
    if p_threshold is not None:
        ATAC_seq = ATAC_seq[ATAC_seq["_-log10_bestPvalue"] >= p_threshold]
        
    
    ATAC_seq_T = ATAC_seq.T
    ATAC_seq_only_scores = ATAC_seq.loc[:,'LTHSC.34-.BM':]
    

    # normalization
    if normalization is None:
        # CPM + log2 normalization
        libsize = ATAC_seq_only_scores.sum(axis=0)
        cpm = ATAC_seq_only_scores.divide(libsize, axis=1) * 1e6
        ATAC_seq_only_scores_norm = np.log2(cpm+1)
        ATAC_seq_only_head = ATAC_seq.loc[:,:'genes.within.100Kb']
        ATAC_seq = pd.concat([ATAC_seq_only_head, ATAC_seq_only_scores_norm], axis=1)
        ATAC_seq_T = ATAC_seq.T

        RNA_seq = np.log2(RNA_seq+1)
        RNA_seq_T = RNA_seq.T
    elif normalization == "noCPM":
        # log2 normalization w/o CPM
        ATAC_seq_only_scores_norm = np.log2(ATAC_seq_only_scores+1)
        ATAC_seq_only_head = ATAC_seq.loc[:,:'genes.within.100Kb']
        ATAC_seq = pd.concat([ATAC_seq_only_head, ATAC_seq_only_scores_norm], axis=1)
        ATAC_seq_T = ATAC_seq.T

        RNA_seq = np.log2(RNA_seq+1)
        RNA_seq_T = RNA_seq.T
    elif normalization == "none":
        # no normalization
        pass

    data = {
        'ATAC_seq': ATAC_seq,
        'ATAC_seq_T': ATAC_seq_T,
        'ATAC_seq_scores_no_norm': ATAC_seq_only_scores,
        'norm_scores': ATAC_seq_only_scores_norm,
        'RNA_seq': RNA_seq,
        'RNA_seq_T': RNA_seq_T, 
        'QC_metrics': QC_metrics,
        'exons': exons,
        'list_ATAC_stem_Tc_Bc': list_ATAC_stem_Tc_Bc,
        'list_ATAC_diff_Tc_all': list_ATAC_diff_Tc_all,
        'list_ATAC_diff_Tc_pre_ab_act': list_ATAC_diff_Tc_pre_ab_act,
        'list_ATAC_diff_Tc_gd': list_ATAC_diff_Tc_gd,
        'list_ATAC_diff_Tc_ab': list_ATAC_diff_Tc_ab,
        'list_ATAC_diff_Tc_ab_gd': list_ATAC_diff_Tc_ab_gd,
        'list_ATAC_Tc_all': list_ATAC_Tc_all,
        'test1': ATAC_copy
        'ab_tc': ab_tc,
        'gd_tc': gd_tc,
        'ab_gd_tc': ab_gd_tc,
    }
    return data

# cluster

# UMAP

# t-SNE
def tSNE (df, cols, components, perplexity, rows=None, gini_coloring=None):
    '''
    no more than 20 000 cols!
    df: pandas.df [rows x cols];
    cols: list of cell row names;
    components: no. of PCAs < len(rows);
    perplexity: perplexity < len(rows);
    gini_coloring: feature that should be colored, automatically maximun, manually after e.g. specifig cell / TSS
    '''
    
    # preparation of data
    df_peaks = df.select_dtypes(include=[np.number])
    if rows is not None:
        subset_df = df_peaks.loc[rows, cols]
    else:
        subset_df = df_peaks.loc[:,cols]
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
    # if gini_coloring is None:
        # name_gini_coloring = gini_scores.idxmax()
    #elif isinstance(gini_coloring, int):
        #name_gini_coloring = subset_df.columns[gini_coloring]
    #else:
       #name_gini_coloring = gini_coloring

    #color_values = subset_df[name_gini_coloring].values

    if gini_coloring is None:
        color_values = 'lightgrey'
        colorbar_label = None
        plot_title = f't-SNE (keine Gini-Färbung)'
        scatter_kwargs = dict(color=color_values, alpha=0.7)
    elif gini_coloring == "TSS":
        # TSS peaks blue, rest gray
        if "TSS" in df.columns:
            mask = df.loc[subset_df.index, "TSS"] != ''
            color_values = np.where(mask, "royalblue", "lightgrey")
            colorbar_label = None
            plot_title = f't-SNE: Peaks in TSS-Nähe (blau)'
            scatter_kwargs = dict(c=color_values, alpha=0.7)
        else:
            raise ValueError('Spalte "TSS" nicht im DataFrame!')
    else:
        # color specific gene red (TSS-col)
        if "TSS" in df.columns:
            mask = df.loc[subset_df.index, "TSS"] == gini_coloring
            color_values = np.where(mask, "crimson", "lightgrey")
            colorbar_label = None
            plot_title = f't-SNE: Peaks mit TSS für {gini_coloring} (rot)'
            scatter_kwargs = dict(c=color_values, alpha=0.7)
        else:
            raise ValueError('Spalte "TSS" nicht im DataFrame!')

    # plot
    plt.figure(figsize=(8,6))
    sc = plt.scatter(tsne_df['tSNE1'], tsne_df['tSNE2'], **scatter_kwargs)
    if colorbar_label is not None:
        plt.colorbar(sc, label=colorbar_label)
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title(plot_title)
    plt.tight_layout()
    plt.show()

    return tsne_df, gini_scores

# correlation_pearson

# linear regression