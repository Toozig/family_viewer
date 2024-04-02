import pandas as pd
import pyBigWig
from viewer_function import make_plot
import time

def timer_decorator(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"Execution time of {func.__name__}: {end_time - start_time} seconds")
        return result
    return wrapper

MAX_SCORE_TFBS = 'maximal score TFBS per site'
SHOW_SEQ = 'Show sequence'

CONSERVATION_FILE = 'data/hg38.phastCons7way.atacIntervals.bw'
# TFBS Files
JASPAR_FILE = 'data/jaspar2024TFBSresults_f350.prq'

HOMER_FILE = 'data/homer_TFBS_results_f7581.prq'


PEAK_FILE = 'data/merged_mATAC_hATAC_0507.bed.gz'
CONFIG_FILE = 'data/config.ini'
SAMPLE_META_DATA_FILE = 'data/sample_metadata.csv'
VARIANT_INHERITENCE_FILE = 'data/all_vars.csv'

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

#### table cols ARGS
SEGMENT_ID_COL = 14
LENGTH_COL = 18
TAD_COL =11
PROBAND_NAMES_COL = 22
HEALTHY_NAMES_COL = 24




## dfs

def process_df(df):
    df['chr'] = df['chr'].str.replace('chr','')
    df['strand'] = df.strand.replace({'-':-1,'+':1})
    df['score'] = df['score'].astype(float)
    return df


BW = pyBigWig.open(CONSERVATION_FILE)


# SOURCE_DICT = {'JASPAR': JASPER_DF,
#                 'HOMER': HOMER_DF}

SOURCE_DICT = {'JASPAR': JASPAR_FILE,
                'HOMER': HOMER_FILE}

SCORE_DEFAULT_THRESHOLD = {'JASPAR': 400,
                           'HOMER': 7}


def filter_chrom_df(source_df, chrom_peak_df):
    chrom_tfbs = source_df[source_df.chr == chrom_peak_df.CHROM.values[0]]
    
    # Create 2D arrays for start and end values
    start = chrom_tfbs['start'].values[:, None]
    end = chrom_tfbs['end'].values[:, None]
    
    # Check for overlaps using broadcasting
    overlaps = (start >= chrom_peak_df['from'].values) & (end <= chrom_peak_df['to'].values)
    
    # Check if any overlap exists for each row
    is_relevant = overlaps.any(axis=1)
    
    # Filter the DataFrame
    f_tfbs = chrom_tfbs[is_relevant]
    
    return f_tfbs




@timer_decorator
def open_family_matadata():
    family_metadata = pd.read_csv(SAMPLE_META_DATA_FILE)
    family_metadata.family_id = family_metadata['family_id'].astype(str)
    return family_metadata

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



class currentState:
    _instance = None
    family_list = None


    @classmethod
    def get_family_list(cls):
        if cls.family_list is None:
            family_metadata =  open_family_matadata()
            cls.family_list = family_metadata.family_id.unique().tolist()
        return cls.family_list


    @classmethod
    def get_source_list(cls):
        return list(SOURCE_DICT.keys())

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(currentState, cls).__new__(cls)
            cls._instance.score_threshold = 400
            cls._instance.family_id = None
            cls._instance.family_df = None
            # init by set_family
            cls._instance.peak_id = None
            cls._instance.source_df = None
            cls._instance.var_df = None
            source_list = cls.get_source_list()
            cls._instance.source = source_list[0]
            # cls._instance.set_source(source_list[0])
            cls._instance.set_family_id('10726')

        return cls._instance

    def set_score_threshold(self, score_threshold):
        self.score_threshold = score_threshold

    def set_source(self, source):
        if source != self.source:
            print("setting source")
            self.source = source
            self.score_threshold = SCORE_DEFAULT_THRESHOLD[source]
            self.__set_source_df(source)

    @timer_decorator
    def __set_source_df(self, source):
        source_df = process_df(pd.read_parquet(SOURCE_DICT[source]))
        print(f'setting source_df, shape: {source_df.shape}')
        cur_peaks = self.var_df[['INTERVAL_ID','CHROM','from','to']].drop_duplicates()
        source_df = cur_peaks.groupby('CHROM').apply(lambda x : filter_chrom_df(source_df, x))
        self.source_df = source_df

    def get_source_df(self):
        if self.source is None:
            print('no source')
            return pd.DataFrame()
        return self.source_df

    def set_family_id(self, family_id):
        """
        Setting the family_id will update the var_df, family_df and peak_id.
        The order of the update is important-
        1. family_df - get the family metadata
        2. var_df - get the variants that are in the family and not in the unaffected #todo - add the option to filter by probands
        3. peak_id - get the peaks that are in the var_df
        """
        cur_family_id = str(family_id)
        if cur_family_id != self.family_id:
            self.family_id = cur_family_id
            self.family_df = self.__set_family_metadata()
            self.__set_full_variant_df()
            self.peak_id = self.get_peak_list()[0]
            self.__set_source_df(self.source)

    def __set_family_metadata(self):
        if self.family_id is None:
            print('no family_id')
            return pd.DataFrame()
        familiy_df = open_family_matadata()
        if self.family_list is None:
            self.family_list = familiy_df.family_id.unique().tolist()
        family_df = familiy_df[familiy_df.family_id == self.family_id]
        return family_df

    def get_family_metadata(self):
        return self.family_df

    def set_peak_id(self, peak_id):
        if peak_id not in self.get_peak_list():
            print(f'peak_id: {peak_id} not in the list')
            return
        self.peak_id = peak_id


    @timer_decorator
    def __set_full_variant_df(self):
        family_df = self.get_family_metadata()
        #get the probands and unaffecteds
        probands = family_df[family_df['fam_relation'] == PROBAND_RELATION].ID
        unaffected = family_df[family_df['fam_relation'] == UNAFFECTED_SIBILING].ID
        # create maske of all variants where all prband have and can not be found in unaffected
        variant_df = pd.read_csv(VARIANT_INHERITENCE_FILE)
        variant_df.CHROM = variant_df.CHROM.str.replace('chr','')
        mask_probands = get_variant_mask(variant_df,probands, 'all')
        unaffected_mask = get_variant_mask(variant_df, unaffected, 'any')
        mask = mask_probands & ~unaffected_mask
        result =  variant_df[mask]
        self.var_df = result

    def get_tfbs_filtered_df(self, one_per_site=True):
        source_df = self.get_source_df()
        peak = self.get_peak_data().to_dict('records')[0]
        tfbs_df = source_df[(source_df.chr == peak['CHROM']) & (source_df.start >= peak['from']) & (source_df.end <= peak['to'])]
        tfbs_df = tfbs_df[tfbs_df.score >= self.score_threshold]
        if one_per_site:
            tfbs_df = max_score_per_coord(tfbs_df)
        return tfbs_df

    
    def get_variant_df(self,  to_filter= True):
        var_df = self.var_df[self.var_df.INTERVAL_ID == self.peak_id].copy()
        print (f'getting variant_df.family_id: {self.family_id}, peak_id: {self.peak_id}')
        if to_filter:
            family_ids = self.get_family_metadata().ID
            gt_cols = [f'{sample_id}:GT' for sample_id in family_ids]
            relevant_cols = VAR_DF_COLS_TO_KEEP + gt_cols
            var_df = var_df[relevant_cols] 
        print(f'shape of var_df: {self.var_df.shape}')
        return var_df.reset_index(drop=True)


    def get_peak_list(self):
        print(f'updating peak list for family_id: {self.family_id}')
        print(f'found {self.var_df.shape[0]} variants')
        print(f'found {self.var_df.INTERVAL_ID.unique().shape[0]} peaks')
        return self.var_df.INTERVAL_ID.unique().tolist()


    def __sum_names(self,names_col):
        names = set()
        for i in names_col:
            names = names |  set(i.split(';'))
        return len(names)



    def get_peak_data(self):
        print(f'getting peak data for peak_id: {self.peak_id}')
        cols = self.var_df.columns
        columns = ['INTERVAL_ID', 'CHROM']
        columns += cols[SEGMENT_ID_COL:LENGTH_COL].tolist() + cols[[TAD_COL]].tolist()
        columns = [i for i in columns if i in cols and i != 'median_DP']
        segment_df = self.get_variant_df(to_filter=False)[columns].copy()
        segment_df.loc[:,'n_probands'] = self.__sum_names(cols[PROBAND_NAMES_COL])
        segment_df.loc[:,'n_healthy'] = self.__sum_names(cols[HEALTHY_NAMES_COL])
        result = pd.DataFrame(segment_df.iloc[0,: ])
        print (f'shape of segment_df: {result.shape}')
        return result.T


    def get_threshold_min_max(self):
        cur_peak_tfbs = self.get_tfbs_filtered_df()
        if cur_peak_tfbs.shape[0] == 0:
            print(f'no tfbs for peak_id: {self.peak_id} of source {self.source}')
            return 0,0
        return cur_peak_tfbs['score'].min(), cur_peak_tfbs['score'].max()


    
    def __parse_names(self,names, header):
        if type(names) != str :
            return ''
        string = f"**{header}**\n"
        string += ''.join(['- ' + name + '\n' for name in names.split(';')])
        return string + '\n'




    def get_variant_info(self,var_idx):
        var_df = self.get_variant_df(to_filter=False)
        variant = var_df.iloc[var_idx]
        result = ''
        for col in ['probands_names', 'healthy_names']:
            result  += self.__parse_names(variant[col], col.split('_')[0].capitalize())
        return result
    
    
    def get_peak_conservation(self, peak_dict=None):
        peak = peak_dict if peak_dict is not None else self.get_peak_data().to_dict('records')[0]
        chrom = peak['CHROM'] if 'chr' in peak['CHROM'] else f'chr{peak["CHROM"]}'
        scores = BW.values(chrom, peak['from'], peak['to'])
        return scores
    
    def get_peak_plot(self, n_lines,checked_box):
        one_per_site = MAX_SCORE_TFBS in checked_box
        tfbs_df = self.get_tfbs_filtered_df(one_per_site)
        var_df = self.get_variant_df()
        peak_dict = self.get_peak_data().to_dict('records')[0]
        conservation_list = self.get_peak_conservation(peak_dict)
        show_seq = SHOW_SEQ in checked_box
        print(f'getting plot')
        return make_plot(CONFIG_FILE, peak_dict)

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
