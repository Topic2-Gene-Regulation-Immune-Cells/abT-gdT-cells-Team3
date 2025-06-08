# data_clean_up: vorläufig
def call_data_clean(p_threshold=None, qc_thresholds=None, normalization=None):
    ''' 
    loads data set;
    p_threshold: 0 < val > 5;
    qc_thresholds: input of qc_threshold_metrics;
    normalization: "none" or log2 automatically
    '''

    # loading of data sets
    ATAC_seq = pd.read_csv(ATAC_seq_path, keep_default_na=False)
    RNA_seq = pd.read_csv(RNA_seq_path)
    QC_metrics = pd.read_excel(Cell_population_qc_path)
    exons=pd.read_csv(Transcription_exons_path, sep="\t")
    exons.columns = ["Genname", "ID", "Chromosom", "Strand", "Transcription Start", "Transcription End", "CDS_Start", "CDS_End", 
    "ExonCount", "ExonStarts", "ExonEnds"]
    
    # lisiting for sub sets and transponing of ATAC
    list_ATAC_stem_Tc_Bc = list(ATAC_seq.loc['LTHSC.34-.BM': 'proB.CLP.BM']),
    list_ATAC_diff_Tc_all= list(ATAC_seq.loc['preT.DN1.Th':'NK.27+11b-.BM']),
    list_ATAC_diff_Tc_pre_ab_act = list(ATAC_seq.loc['preT.DN1.Th':'Tgd.g2+d17.24a+.Th']),
    list_ATAC_diff_Tc_gd = list(ATAC_seq.loc['Tgd.g2+d17.24a+.Th':'NK.27+11b-.BM']),
    list_ATAC_Tc_all = list_ATAC_stem_Tc_Bc + list_ATAC_diff_Tc_all
    
    list_cells = {
        'list_ATAC_stem_Tc_Bc': list_ATAC_stem_Tc_Bc,
        'list_ATAC_diff_Tc_all': list_ATAC_diff_Tc_all,
        'list_ATAC_diff_Tc_pre_ab_act': list_ATAC_diff_Tc_pre_ab_act,
        'list_ATAC_diff_Tc_gd': list_ATAC_diff_Tc_gd,
        'list_ATAC_Tc_all': list_ATAC_Tc_all
    }

    ATAC_seq_T = ATAC_seq.T
    ATAC_seq_only_scores = ATAC_seq.loc['LTHSC.34-.BM':]

    # normalization
    if normalization is None:
        # log2 normalization
        ATAC_seq_only_scores_norm = np.log2(ATAC_seq_only_scores)
        ATAC_seq_only_head = ATAC_seq.loc[:'LTHSC.34-.BM']
        ATAC_seq = pd.concat([ATAC_seq_only_head, ATAC_seq_only_scores_norm])
        ATAC_seq_T = ATAC_seq.T
    elif normalization == "none":
        # no normalization
        pass

    #thresholds

    #summary

    results = {
        'ATAC_seq': ATAC_seq,
        'ATAC_seq_T': ATAC_seq_T,
        'ATAC_seq_only_scores': ATAC_seq_only_scores,
        'list_cells': list_cells,
        'RNA_seq': RNA_seq,
        'QC_metrics': QC_metrics,
        'exons': exons,
    }
    return results

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