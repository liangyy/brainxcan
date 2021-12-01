import re
# def get_emp_null(fn):
#     return re.sub('.csv$', '.emp_null.csv', fn)
# def get_perm_null(fn):
#     return re.sub('.csv$', '.perm_null.csv', fn)
def get_null_file(fn, name):
    return re.sub('.csv$', f'.{name}.csv', fn)
def append_null_to_df(df, fn, null_name):
    null_fn = get_null_file(fn, null_name)
    if os.path.exists(null_fn):
        df_null = pd.read_csv(null_fn)
        print(f'Loading {null_name} from {fn}: num. of records = {df_null.shape[0]}')
        df.append(df_null)
def add_null_to_result(df_result, df_null, null_name, output_prefix):
    if len(df_null) > 0:
        df_null = pd.concat(df_null, axis=0).reset_index(drop=True)
        print(f'Forming adjusted z-score using {null_name} as null (n = {df_null.shape[0]})')
        df_null.to_csv(output_prefix + f'.{null_name}.csv', index=False)
        varz_null = np.var(df_null.value)
        zcol, pcol = f'z_adj_{null_name}', f'pval_adj_{null_name}'
        df_result[zcol] = df_result.z_brainxcan / np.sqrt(varz_null)
        df_result[pcol] = z2p(df_result[zcol])  

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='merge_bxcan.py', description='''
        Merge BrainXcan results from dMRI and T1.
    ''')
    parser.add_argument('--dmri', help='''
        Input S-BrainXcan result for dMRI IDPs.
    ''')
    parser.add_argument('--t1', help='''
        Input S-BrainXcan result for T1 IDPs.
    ''')
    parser.add_argument('--idp_meta_file', help='''
        A meta file for annotating IDPs.
    ''')
    parser.add_argument('--output_prefix', help='''
        Output CSV table prefix.
    ''')
    args = parser.parse_args()
    
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    import pandas as pd
    import numpy as np
    from brainxcan.sbxcan.util.misc import z2p
    from brainxcan.sbxcan.run_sbrainxcan import genomic_control
    
    df_null = []
    df_perm = []
    logging.info('Loading S-BrainXcan dMRI.')
    df1 = pd.read_csv(args.dmri)
    df1['modality'] = 'dMRI'
    logging.info('{} IDPs in total.'.format(df1.shape[0]))
    append_null_to_df(df_null, args.dmri, 'emp_null')
    append_null_to_df(df_perm, args.dmri, 'perm_null')
    
    
    logging.info('Loading S-BrainXcan T1.')
    df2 = pd.read_csv(args.t1)
    df2['modality'] = 'T1'
    logging.info('{} IDPs in total.'.format(df2.shape[0]))
    append_null_to_df(df_null, args.t1, 'emp_null')
    append_null_to_df(df_perm, args.t1, 'perm_null')
    
    logging.info('Generating adjusted BrainXcan z-score.')
    df = pd.concat([df1, df2], axis=0)
    df['z_adj_gc'], lambda_gc = genomic_control(df.z_brainxcan)
    df['pval_adj_gc'] = z2p(df.z_adj_gc)
    logging.info(f'GC lambda = {lambda_gc}.')
    
    add_null_to_result(df, df_null, 'emp_null', args.output_prefix)
    add_null_to_result(df, df_perm, 'perm_null', args.output_prefix)
 
    logging.info('Loading the IDP meta file.')
    meta = pd.read_csv(args.idp_meta_file)
    
    logging.info('Saving outputs.')
    df = pd.merge(df, meta.drop(columns=['t1_or_dmri', 'ukb_link']), on='IDP', how='left')
    df.fillna('NA', inplace=True)
    df.sort_values(by='pval').to_csv(args.output_prefix + '.csv', index=False)
    
    logging.info('Done.')
