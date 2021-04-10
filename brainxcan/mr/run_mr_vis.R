library(optparse)

option_list <- list(
    make_option(c("-i", "--mr_rds"), type="character", default=NULL,
                help="Input MR RDS (output of run_mr_local.R)",
                metavar="character"),
    make_option(c("-b", "--bxcan"), type="character", default=NULL,
                help="BrainXcan result (single line).",
                metavar="character"),
    make_option(c("-p", "--output_plot"), type="character", default=NULL,
                help="Output MR plot.",
                metavar="character"),
    make_option(c("-t", "--output_table"), type="character", default=NULL,
                help="Output MR table.",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

# srcpath = Sys.getenv("R_SRC")
# source(paste0(srcpath, '/', 'mr_helper.R'))
library(dplyr)
library(ggplot2)
library(patchwork)
options(stringsAsFactors = F)
theme_set(theme_classic(base_size = 12))

# logging config
logging::basicConfig(level = opt$log_level)

logging::loginfo('Loading MR RDS and BrainXcan 1-line result.')
mr_res = readRDS(opt$mr_rds)
bxcan = read.csv(opt$bxcan)

logging::loginfo('Plotting.')
p1 = mr_res$idp2pheno$data %>% ggplot() + geom_hline(yintercept = 0, color = 'gray') + 
  geom_vline(xintercept = 0, color = 'gray') + 
  geom_point(aes(x = beta.exposure, y = beta.outcome), alpha = 0.5) + 
  geom_errorbar(aes(x = beta.exposure, ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), alpha = 0.5) +
  geom_errorbarh(aes(y = beta.outcome, xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), alpha = 0.5) +
  xlab('SNP estimated effect in IDP') + ylab('SNP estimated effect in phenotype') + ggtitle('MR: IDP -> Phenotype')
p2 = mr_res$pheno2idp$data %>% ggplot() + geom_hline(yintercept = 0, color = 'gray') + 
  geom_vline(xintercept = 0, color = 'gray') + 
  geom_point(aes(x = beta.exposure, y = beta.outcome), alpha = 0.5) + 
  geom_errorbar(aes(x = beta.exposure, ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), alpha = 0.5) +
  geom_errorbarh(aes(y = beta.outcome, xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), alpha = 0.5) +
  xlab('SNP estimated effect in phenotype') + ylab('SNP estimated effect in IDP') + ggtitle('MR: Phenotype -> IDP')
ggsave(opt$output_plot, p1 + p2, height = 4, width = 8)

logging::loginfo('Collecting summary of results.')
mr_methods = c('Inverse variance weighted', 'Weighted median', 'MR Egger')
df_mr = rbind(
  mr_res$idp2pheno$mr %>% filter(method %in% mr_methods) %>% mutate(direction = 'IDP -> Phenotype'),
  mr_res$pheno2idp$mr %>% filter(method %in% mr_methods) %>% mutate(direction = 'Phenotype -> IDP') 
) %>% select(direction, method, nsnp, b, pval)
bxcan = bxcan %>% mutate(direction = NA, method = 'BrainXcan') %>% rename(nsnp = nsnp_used, b = bhat)
df_res = rbind(
  df_mr, 
  bxcan %>% select(direction, method, nsnp, b, pval)
)
write.table(df_res, opt$output_table, sep = '\t', col = T, row = F, quo = F)

logging::loginfo('Done.')

