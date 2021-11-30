import re
def get_emp_null(fn):
    return re.sub('.csv$', '.emp_null.csv', fn)
    
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
    from brainxcan.sbxcan.run_sbrainxcan import p2z, genomic_control
    
    df_null = []
    logging.info('Loading S-BrainXcan dMRI.')
    df1 = pd.read_csv(args.dmri)
    df1['modality'] = 'dMRI'
    logging.info('{} IDPs in total.'.format(df1.shape[0]))
    null_dmri = get_emp_null(args.dmri)
    if os.path.exists(null_dmri):
        null_dmri = pd.read_csv(null_dmri)
        logging.info('Loading empirical nulls from dMRI: nrepeat = {null_dmri.shape[0]}')
        df_null.append(null_dmri)
    
    logging.info('Loading S-BrainXcan T1.')
    df2 = pd.read_csv(args.t1)
    df2['modality'] = 'T1'
    logging.info('{} IDPs in total.'.format(df2.shape[0]))
    null_t1 = get_emp_null(args.t1)
    if os.path.exists(null_t1):
        null_t1 = pd.read_csv(null_t1)
        logging.info('Loading empirical nulls from T1: nrepeat = {null_t1.shape[0]}')
        df_null.append(null_t1)
    
    logging.info('Generating adjusted BrainXcan z-score.')
    df = pd.concat([df1, df2], axis=0)
    df['z_adj_gc'], lambda_gc = genomic_control(df.z_brainxcan)
    logging.info(f'GC lambda = {lambda_gc}.')
    
    if len(df_null) > 0:
        df_null = pd.concat(df_null, axis=0).reset_index(drop=True)
        df_null.to_csv(args.output_prefix + '.null.csv', index=False)
        varz_null = np.var(df_null.value)
        df['z_adj_emp'] = df.z_brainxcan / np.sqrt(varz_null)
    
    logging.info('Loading the IDP meta file.')
    meta = pd.read_csv(args.idp_meta_file)
    
    logging.info('Saving outputs.')
    df = pd.merge(df, meta.drop(columns=['t1_or_dmri', 'ukb_link']), on='IDP', how='left')
    df.fillna('NA', inplace=True)
    df.sort_values(by='pval').to_csv(args.output_prefix + '.csv', index=False)
    
    logging.info('Done.')
