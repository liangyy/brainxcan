if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='grab_signif.py', description='''
        Grab the significant IDPs from S-BrainXcan results.
    ''')
    parser.add_argument('--sbxcan', help='''
        Input S-BrainXcan result.
    ''')
    parser.add_argument('--pval', type=float, help='''
        P-value cutoff to call significance.
    ''')
    parser.add_argument('--max_idps', type=int, help='''
        Maximum number of IDPs to output.
        If there are more than max_idps number of IDPs passing 
        the criteria, the top max_idps will be output.
    ''')
    parser.add_argument('--pval_col', default='pval', help='''
        The name of the p-value column to be used.
    ''')
    parser.add_argument('--outdir', help='''
        Output directory
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

    logging.info('Loading S-BrainXcan.')
    df = pd.read_csv(args.sbxcan)
    logging.info('{} IDPs in total.'.format(df.shape[0]))
    
    logging.info('Applying criteria: p-value cutoff = {} and pval column = {}'.format(args.pval, args.pval_col))
    df = (
        df[ df[args.pval_col] < args.pval ]
            .sort_values(by=args.pval_col, ascending=True)
            .reset_index(drop=True)
    )
    logging.info('{} IDPs pass the criteria.'.format(df.shape[0]))
    
    if df.shape[0] > args.max_idps:
        logging.info('Taking top {} IDPs.'.format(args.max_idps))
        # keep all ties 
        pval_last_to_select = df.iloc[:args.max_idps, :][args.pval_col].values[-1]
        df = df[ df[args.pval_col] <= pval_last_to_select ].reset_index(drop=True)
        logging.info(
            'Select {} IDPs with p-value <= {}'.format(
                df.shape[0], pval_last_to_select
            )
        )
    
    logging.info('Saving outputs.')
    for idp, modality in zip(df.IDP, df.modality):
        tmp = df[ df.IDP == idp ]
        tmp.to_csv('{}/{}_{}.txt'.format(args.outdir, modality.lower(), idp), index=False)

    
    
    logging.info('Done.')
