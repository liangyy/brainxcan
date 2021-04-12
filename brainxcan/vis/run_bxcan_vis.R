library(optparse)

option_list <- list(
    make_option(c("-b", "--brainxcan"), type="character", default=NULL,
                help="BrainXcan result (merged CSV).",
                metavar="character"),
    make_option(c("-p", "--output"), type="character", default=NULL,
                help="Output plot filename",
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

tags = names(readRDS(paste0(opt$datadir, '/meta_plot.rds')))

logging::loginfo('Loading BrainXcan results.')
df = read.csv(opt$brainxcan) %>% mutate(zscore = p2z(pval, bhat))

logging::loginfo('Plotting for all tags.')
p = list() 
for(tag in tags) {
  logging::loginfo(paste0('Tag = ', tag))
  p[[length(p) + 1]] = vis_by_tag(opt$datadir, tag, df, 'zscore') + ggtitle(tag)
}
pp = wrap_plots(p, ncol = 1, nrow = length(tags))
ggsave(opt$output, pp, height = 4 * length(tags), width = 10)

