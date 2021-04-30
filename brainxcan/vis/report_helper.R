load_color_code_yaml = function(fn) {
  kk = yaml::read_yaml(fn)
  res = c()
  for(n in names(kk)) {
    if(length(kk[n]) > 1) {
      message('More than one color code, will load the first one.')
    }
    res = c(res, kk[[n]][1])
  }
  names(res) = names(kk)
  res
}

plot_bxcan_ordered <- function(xcandf, color_map)
{
  pp = xcandf %>% 
    mutate(kk = reorder(region, zscore ^ 2, FUN = max)) %>% 
    group_by(kk) %>% mutate(zmax = max(zscore), zmin = min(zscore)) %>% ungroup() %>%  
    ggplot() + 
    geom_hline(yintercept = c(-z_thres, z_thres), col = 'gray') + 
    geom_segment(aes(x = kk, xend = kk, y = zmin, yend = zmax), color = 'black', size = 0.1) + 
    geom_point(aes(x = kk, y = zscore, color = subtype)) + 
    coord_flip() +
    scale_color_manual(values = color_map) + 
    ylab('BrainXcan z-score') + 
    theme(axis.title.y = element_blank()) + 
    theme(legend.title = element_blank())
  pp
}

load_mr_sum_stats = function(prefix) {
  kk = Sys.glob(paste0(prefix, '*.MR_sumstats.tsv'))
  df_mr = list()
  for(k in kk) {
    mm = stringr::str_match(basename(k), '(t1_|dmri_)([A-Z-0-9a-z]+).MR_sumstats.tsv')[, 3]
    tmp = read.delim2(k)
    tmp$IDP = mm
    df_mr[[length(df_mr) + 1]] = tmp
  }
  df_mr = do.call(rbind, df_mr)
  df_mr$b = as.numeric(df_mr$b)
  df_mr$pval = as.numeric(df_mr$pval)
  df_mr
}

gen_mr_sign = function(direction, pval) {
  df_code = data.frame(
    direction = c(1, -1),
    sign = c('\u25B2', '\u25BC')
  )
  kk = left_join(
    data.frame(direction = sign(direction), pval = pval),
    df_code, by = 'direction'
  )
  paste0(signif(kk$pval, digits = 3), ' (', kk$sign, ')')
}
