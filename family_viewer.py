import pandas as pd
import pyBigWig

CONSERVATION_FILE = 'data/hg38.phastCons7way.atacIntervals.bw'
# CONSERVATION_FILE = '/dsi/gonen-lab/users/roni/useful_files/hg38.phastCons7way.atacIntervals.bw'
# TFBS Files
JASPAR_FILE = 'data/jaspar2024TFBSresults.prq'
HOMER_FILE = 'data/homer_TFBS_results.prq'
# JASPAR_FILE = '/dsi/gonen-lab/users/roni/useful_files/outputs/jaspar2024TFBSresults.prq'
# HOMER_FILE = '/dsi/gonen-lab/users/roni/useful_files/outputs/homer_TFBS_results.prq'
# VARIANT_INHERITENCE_FILE = '/dsi/gonen-lab/shared_files/WGS_on_DSD/data/pipeline_outputs/variants_with_layers/2024-01-11/inheritance/all_vars.csv'
PEAK_FILE = 'data/merged_mATAC_hATAC_0507.bed.gz'
SAMPLE_META_DATA_FILE = 'data/sample_metadata.xlsx'
VARIANT_INHERITENCE_FILE = 'data/all_vars.csv'
# PEAK_FILE = '/dsi/gonen-lab/users/toozig/projects/WT_canidate_enhancer/data/merged_mATAC_hATAC_0507.bed.gz'
# SAMPLE_META_DATA_FILE = '/dsi/gonen-lab/shared_files/WGS_on_DSD/data/read_only/samples/sample_metadata.xlsx'
# VARIANT_INHERITENCE_FILE = '/dsi/gonen-lab/users/roni/DSDncVariants/variants_pipeline/2024-01-11/inheritance/all_vars.csv'


PROBAND_RELATION = 0
UNAFFECTED_SIBILING = 3
HETRO = ["1/0", "0/1", "1|0", "0|1"]
HOMO_REF = ["0/0", "0|0", "", " ", None]
HOMO_ALT = ["1/1", "1|1"]
VAR_DF_COLS_TO_KEEP = ['CHROM',
 'POS',
 'REF',
 'ALT',
 'AF',
 'AF_popmax',
 'TAD',
 'total_probands',
 'healthy_members',]



## dfs

def process_df(df):
    df['strand'] = df.strand.replace({'-':-1,'+':1})
    df['score'] = df['score'].astype(float)
    return df


JASPER_DF = pd.read_parquet(JASPAR_FILE) 
HOMER_DF = pd.read_parquet(HOMER_FILE)
VARIANT_DF = pd.read_csv(VARIANT_INHERITENCE_FILE)
VARIANT_DF.CHROM = VARIANT_DF.CHROM.str.replace('chr','')
jasper_df = process_df(JASPER_DF)
homer_df = process_df(HOMER_DF)
SAMPLE_META_DATA_DF =  pd.read_excel(SAMPLE_META_DATA_FILE)
SAMPLE_META_DATA_DF.family_id = SAMPLE_META_DATA_DF['family_id'].astype(str)
PEAK_DF = pd.read_csv(PEAK_FILE, sep='\t', header=None,  names=['chrom','start', 'end', 'id'])
BW = pyBigWig.open(CONSERVATION_FILE)


SOURCE_DICT = {'JASPAR': JASPER_DF,
                'HOMER': HOMER_DF}


def get_gt_column_name(sample_id):
    """
    Generate a genotype column name for a given sample ID.

    Args:
        sample_id: The ID of the sample.

    Returns:
        A string representing the genotype column name.
    """
    return f'{sample_id}:GT'

def filter_rows_any_var(row, gt_cols):
    """
    Filter rows where any of the provided columns have a variant.

    Args:
        row: The row to be checked.
        gt_cols: The genotype columns to be checked.

    Returns:
        True if any of the columns have a variant, False otherwise.
    """
    return any(row[col] not in HOMO_REF for col in gt_cols)

def filter_rows_all_vars(row, gt_cols):
    """
    Filter rows where all of the provided columns have a variant.

    Args:
        row: The row to be checked.
        gt_cols: The genotype columns to be checked.

    Returns:
        True if all of the columns have a variant, False otherwise.
    """
    return all(row[col] not in HOMO_REF for col in gt_cols)

def filter_not_in_vars(row, gt_cols):
    """
    Filter rows where none of the provided columns have a variant.

    Args:
        row: The row to be checked.
        gt_cols: The genotype columns to be checked.

    Returns:
        True if none of the columns have a variant, False otherwise.
    """
    return all(row[col] in HOMO_REF for col in gt_cols)

def get_family_variants(probands, variants):
    """
    Get variants dataframe with rows where probands have a variant.

    Args:
        probands: The probands to be checked.
        variants: The variants to be checked.

    Returns:
        A dataframe of variants where probands have a variant.
    """
    gt_cols = probands.apply(get_gt_column_name)
    fam_variants = variants[variants.apply(filter_rows_any_var, axis=1, gt_cols=gt_cols)]
    return fam_variants

def two_probands_one_unaffected(variants, probands, unaffected):
    """
    Get variants dataframe with rows where probands have a variant and unaffected sibling doesn't.

    Args:
        variants: The variants to be checked.
        probands: The probands to be checked.
        unaffected: The unaffected to be checked.

    Returns:
        A dataframe of variants where probands have a variant and unaffected sibling doesn't.
    """
    prob_gt_cols = probands.apply(get_gt_column_name)
    prob_variants = variants[variants.apply(filter_rows_all_vars, axis=1, gt_cols=prob_gt_cols)]
    unaffected_gt_cols = unaffected.apply(get_gt_column_name)
    final_variants = prob_variants[prob_variants.apply(filter_not_in_vars, axis=1, gt_cols=unaffected_gt_cols)]
    return final_variants


def get_family_metadata(family_id):
    
    family_df = SAMPLE_META_DATA_DF[SAMPLE_META_DATA_DF['family_id'] == str(family_id)]
    return family_df




def get_variant_mask(variant_df, ids,second_filter = 'all'):
    if not len(ids):
        # return all False of if list is empty
        print('no ids for mask')
        return pd.Series([False] * len(variant_df))
    col_names = [f'{sample_id}:GT' for sample_id in ids]
    mask = ~(variant_df[col_names].isin(HOMO_REF))
    if second_filter == 'all':
        mask = mask.all(axis=1)
    elif second_filter == 'any':
        mask = mask.any(axis=1)
    return mask


def get_filter_variants_by_peak(family_id, peak_id):
    var_df = get_filter_variant_df(family_id)
    result =  var_df[var_df.INTERVAL_ID == peak_id]
    return result

def get_filter_variant_df(family_id):

    family_df = get_family_metadata(family_id)
    #get the probands and unaffecteds
    probands = family_df[family_df['fam_relation'] == PROBAND_RELATION].ID
    unaffected = family_df[family_df['fam_relation'] == UNAFFECTED_SIBILING].ID

    #get the variants of the family
    # fam_variants = get_family_variants(probands, variant_df)


    # create maske of all variants where all prband have and can not be found in unaffected
    mask_probands = get_variant_mask(VARIANT_DF,probands, 'all')
    unaffected_mask = get_variant_mask(VARIANT_DF, unaffected, 'any')
    mask = mask_probands & ~unaffected_mask

    result =  VARIANT_DF[mask]
    return result



#### table functions:
SEGMENT_ID_COL = 14
LENGTH_COL = 18
TAD_COL =11
PROBAND_NAMES_COL = 22
HEALTHY_NAMES_COL = 24

def sum_names(names_col):
    names = set()
    for i in names_col:
        names = names |  set(i.split(';'))
    return len(names)
    
def get_peak_data(peak_id):
    print(f'getting peak data for peak_id: {peak_id}')
    var_df = VARIANT_DF[VARIANT_DF.INTERVAL_ID == peak_id]
    cols = var_df.columns
    columns = ['CHROM'] + cols[SEGMENT_ID_COL:LENGTH_COL].tolist() + cols[[TAD_COL]].tolist() 
    columns = [i for i in columns if i in cols and i != 'median_DP']
    segment_df = var_df[columns].copy()
    segment_df.loc[:,'n_probands'] = sum_names(cols[PROBAND_NAMES_COL])
    segment_df.loc[:,'n_healthy'] = sum_names(cols[HEALTHY_NAMES_COL])
    result = pd.DataFrame(segment_df.iloc[0,: ])
    print (f'shape of segment_df: {result.shape}')
    return result.T

def get_family_list():
    return SAMPLE_META_DATA_DF.family_id.unique().tolist()

def get_source_list():
    return list(SOURCE_DICT.keys())

def get_peak_list(family_id):
    print(f'updating peak list for family_id: {family_id}')
    var_df = get_filter_variant_df(family_id)
    print(f'found {var_df.shape[0]} variants')
    print(f'found {var_df.INTERVAL_ID.unique().shape[0]} peaks')
    return var_df.INTERVAL_ID.unique().tolist()

def get_variant_df(family_id, peak_id, to_filter= True):
    print (f'getting variant_df.family_id: {family_id}, peak_id: {peak_id}')
    var_df = get_filter_variants_by_peak(str(family_id), peak_id)

    family_ids = get_family_metadata(family_id).ID
    gt_cols = [f'{sample_id}:GT' for sample_id in family_ids]
    relevant_cols = VAR_DF_COLS_TO_KEEP + gt_cols
    # drop_cols = ['FILTER','repeatsMasker']
    # relevant_cols = [i for i in relevant_cols if i not in drop_cols]
    var_df = var_df[relevant_cols] if to_filter else var_df
    print(f'shape of var_df: {var_df.shape}')
    return var_df.reset_index(drop=True)
    








def parse_names(names, header):
    string = f"**{header}**\n"
    string += ''.join(['- ' + name + '\n' for name in names.split(';')])
    return string + '\n'

def get_variant_info(var_idx, family_id, peak_id):
    var_df = get_variant_df(family_id, peak_id, to_filter=False)

    variant = var_df.iloc[var_idx]
    result = ''
    for col in ['probands_names', 'healthy_names']:
        result  += parse_names(variant[col], col.split('_')[0].capitalize())
    return result

def get_peak_conservation(peak_id):
    peak = PEAK_DF[PEAK_DF.id == peak_id].to_dict('records')[0]
    scores = BW.values(peak['chrom'], peak['start'], peak['end'])
    return scores



#### TFBS table function

def max_score_per_coord(df):

    # Create a new column 'group' to represent overlapping segments with the same 'strand'
    groups = (df['start'] <= df['end'].shift(1)).tolist()
    group = 0 
    i = 0
    flag = False
    label = []
    while i < len(df):
        cur = groups[i]
       
        if cur:
            if not flag:
                group += 1
                label[-1] = group
                flag = True

            label.append(group)
        else:
            group += 1
            label.append(group)
            flag = False
        i += 1
    df.insert(0,'label',label)
    # Find the row index with the maximum score for each group
    df = df.astype({'score': 'int32'})
    idx_max_score = df.groupby(['label','strand'])['score'].idxmax()

    # Extract the corresponding rows from the DataFrame
    result_df = df.loc[idx_max_score].drop(columns=['label'])

    return result_df

def get_tfbs_filtered_df(source_name, peak_id, threshold=0, one_per_site=False):
    source_df = SOURCE_DICT[source_name]
    peak = PEAK_DF[PEAK_DF.id == peak_id].to_dict('records')[0]
    tfbs_df = source_df[(source_df.chr == peak['chrom']) & (source_df.start >= peak['start']) & (source_df.end <= peak['end'])]
    tfbs_df = tfbs_df[tfbs_df.score >= threshold]
    if one_per_site:
        tfbs_df = max_score_per_coord(tfbs_df)
    return tfbs_df


def get_threshold_min_max(source_name, peak_id):
    source_df = get_tfbs_filtered_df(source_name, peak_id)
    return source_df['score'].min(), source_df['score'].max()

