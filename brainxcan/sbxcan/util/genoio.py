import pandas as pd
import numpy as np
import scipy.stats
from pandas_plink import read_plink1_bin 


def load_genotype_from_bedfile(bedfile, indiv_list, snplist_to_exclude, chromosome=None, load_first_n_samples=None, 
    missing_rate_cutoff=0.5, return_snp=False, standardize=True):
    G = read_plink1_bin(bedfile, verbose=False)
    
    if chromosome is not None:
        chr_str = G.chrom[0].values.tolist()
        if 'chr' in chr_str:
            chromosome = 'chr' + str(chromosome)
        else:
            chromosome = str(chromosome)
        G = G.where(G.chrom == chromosome, drop=True)
    
    df_geno_indiv = pd.DataFrame({'indiv': G.sample.to_series().tolist()})
    df_geno_indiv['idx'] = [ i for i in range(df_geno_indiv.shape[0]) ]
    
    if indiv_list is None:
        indiv_list = G.sample.to_series().tolist()
        if load_first_n_samples is not None:
            indiv_list = indiv_list[:load_first_n_samples]
    df_target_indiv = pd.DataFrame({'indiv': indiv_list})
    df_geno_indiv = pd.merge(df_geno_indiv, df_target_indiv, on='indiv').sort_values(by=['idx'])
    if df_geno_indiv.shape[0] != len(indiv_list):
        raise ValueError('There are input individuals that do not appear in BED file.')
    query_indiv_list = df_geno_indiv.indiv.tolist()
    
    
    snpid = G.variant.variant.to_series().to_list()
    snpid = np.array([ s.split('_')[1] for s in snpid ])
    if return_snp is True:
        a0 = G.variant.a0.to_series().to_numpy()
        a1 = G.variant.a1.to_series().to_numpy()       
        chrom = G.variant.chrom.to_series().to_numpy()    
    
    geno = G.sel(sample=query_indiv_list).values

    # re-order to target indiv_list
    geno = geno[match_y_to_x(np.array(query_indiv_list), np.array(indiv_list)), :]
    
    # filter out unwanted snps
    geno = geno[:, ~np.isin(snpid, snplist_to_exclude)]
    if return_snp is True:
        a0 = a0[~np.isin(snpid, snplist_to_exclude)]
        a1 = a1[~np.isin(snpid, snplist_to_exclude)]
        chrom = chrom[~np.isin(snpid, snplist_to_exclude)]
        
    snpid = snpid[~np.isin(snpid, snplist_to_exclude)]
   
    # filter out genotypes with high missing rate
    missing_rate = np.isnan(geno).mean(axis=0)
    geno = geno[:, missing_rate < missing_rate_cutoff]
    if return_snp is True:
        snpid = snpid[missing_rate < missing_rate_cutoff]
        a0 = a0[missing_rate < missing_rate_cutoff]
        a1 = a1[missing_rate < missing_rate_cutoff]
        chrom = chrom[missing_rate < missing_rate_cutoff]
        
    maf = np.nanmean(geno, axis=0) / 2
    
    # impute genotype missing value
    miss_x, miss_y = np.where(np.isnan(geno))
    geno[(miss_x, miss_y)] = maf[miss_y] * 2
    var_geno = 2 * maf * (1 - maf)
    
    # keep only genotypes with variance != 0
    to_keep = var_geno != 0
    geno = geno[:, to_keep]
    if return_snp is True:
        snpid = snpid[to_keep]
        a0 = a0[to_keep]
        a1 = a1[to_keep]
        chrom = chrom[to_keep]
        
    maf = maf[to_keep]
    var_geno = var_geno[to_keep]
    if standardize is True:
        geno = (geno - 2 * maf) / np.sqrt(var_geno)
    
    if return_snp is True:
        return geno, indiv_list, np.sqrt(var_geno), (snpid.tolist(), a0.tolist(), a1.tolist(), chrom.tolist())
    else:
        return geno, indiv_list, np.sqrt(var_geno)


def match_y_to_x(x, y):
    '''
    x, y are 1d np.array 
    y is a subset of x.
    return idx such that x[idx] = y
    '''
    return np.where(y.reshape(y.size, 1) == x)[1]
