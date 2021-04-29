library(optparse)

option_list <- list(
    make_option(c("-b", "--brainxcan"), type="character", default=NULL,
                help="BrainXcan result (merged CSV).",
                metavar="character"),
    make_option(c("-p", "--output_prefix"), type="character", default=NULL,
                help="Output plot prefix",
                metavar="character"),
    make_option(c("-d", "--datadir"), type="character", default=NULL,
                help="Meta data directory",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

srcpath = Sys.getenv("R_VIS_SRC")
source(paste0(srcpath, '/', 'vis_helper.R'))

library(dplyr)
library(ggplot2)
library(patchwork)
library(oro.nifti)
theme_set(theme_classic(base_size = 12))
options(stringsAsFactors = F)

meta_list = readRDS(paste0(opt$datadir, '/meta_plot.rds'))
tags = names(meta_list)

logging::loginfo('Loading BrainXcan results.')
df = read.csv(opt$brainxcan) %>% mutate(zscore = p2z(pval, bhat))

logging::loginfo('Plotting for all tags.')
p = list(T1 = list(), dMRI = list()) 
for(tag in tags) {
  logging::loginfo(paste0('Tag = ', tag))
  idp_type = meta_list[[tag]]$type
  p[[idp_type]][[length(p[[idp_type]]) + 1]] = vis_by_tag(opt$datadir, tag, df, 'zscore') + ggtitle(meta_list[[tag]]$full_name)
}
for(kk in names(p)) {
  npp = length(p[[kk]])
  if(npp > 0) {
    pp = wrap_plots(p[[kk]], ncol = 1, nrow = npp)
    ggsave(paste0(opt$output_prefix, '.', kk, '.png'), pp, height = 4 * npp, width = 10)
  }
}



