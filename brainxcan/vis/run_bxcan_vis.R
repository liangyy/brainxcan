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
                metavar="character"),
    make_option(c("-r", "--region_vis"), action="store_true",
              help="Also generate interactive region htmls along the way")
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
h = list(T1 = c(), dMRI = c())
for(tag in tags) {
  logging::loginfo(paste0('Tag = ', tag))
  idp_type = meta_list[[tag]]$type
  ptmp = vis_by_tag(opt$datadir, tag, df, 'zscore')
  ptmp = add_title(ptmp, meta_list[[tag]]$full_name, size = 20, hjust = 0.2)
  p[[idp_type]][[length(p[[idp_type]]) + 1]] = ptmp
  h[[idp_type]] = c(h[[idp_type]], ceiling(length(unique(ptmp$data$direction)) / 3))
}
for(kk in names(p)) {
  npp = length(p[[kk]])
  if(npp > 0) {
    pp = wrap_plots(p[[kk]], ncol = 1, nrow = npp, heights = h[[kk]])
    ggsave(paste0(opt$output_prefix, '.', kk, '.pdf'), pp, height = 3 * sum(h[[kk]]), width = 10)
  }
}


if(isTRUE(opt$region_vis)) {
  dirn = paste0(opt$output_prefix, '.region_vis')
  dir.create(dirn)
  logging::loginfo('Plotting region labels.')
  df = read.csv(paste0(opt$datadir, '/idp_meta_data.csv'))
  tbss = F
  # plab = list()
  for(tag in tags) {
    logging::loginfo(paste0('Tag = ', tag))
    title = meta_list[[tag]]$full_name
    save_name = tag
    if(substr(tag, 1, 4) == 'TBSS') {
      if(isTRUE(tbss)) {
        next
      } else {
        tbss = T
        title = 'Diffusion MRI'
        save_name = 'TBSS'
      }
    }
    p = vis_region(opt$datadir, tag, df)
    p = p + ggtitle(title)
    fig = plotly::ggplotly(p)
    htmlwidgets::saveWidget(plotly::as_widget(fig), paste0(dirn, '/label_', save_name, '.html'), selfcontained = T)
  }
}




