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
    parser.add_argument('--output', help='''
        Output table.
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

    logging.info('Loading S-BrainXcan dMRI.')
    df1 = pd.read_csv(args.dmri)
    df1['modality'] = 'dMRI'
    logging.info('{} IDPs in total.'.format(df1.shape[0]))
    
    logging.info('Loading S-BrainXcan T1.')
    df2 = pd.read_csv(args.t1)
    df2['modality'] = 'T1'
    logging.info('{} IDPs in total.'.format(df2.shape[0]))
    
    logging.info('Saving outputs.')
    df = pd.concat([df1, df2], axis=0)
    df.sort_values(by='pval').to_csv(args.output, index=False)
    
    logging.info('Done.')
