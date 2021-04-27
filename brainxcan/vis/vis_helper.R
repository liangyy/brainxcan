p2z = function(p, b) {
  sign(b) * abs(qnorm(p / 2, lower.tail = T))
}

gen_five_scales = function(sig_cutoff, nominal_cutoff) {
  res = c(
    paste0('z > ', sig_cutoff),
    paste0(nominal_cutoff, ' < z <= ', sig_cutoff),
    paste0(-nominal_cutoff, ' <= z <= ', nominal_cutoff),
    paste0(-sig_cutoff, ' <= z < ', -nominal_cutoff),
    paste0('z < ', - sig_cutoff)
  )
  return(res)
}

zscore2category = function(zval, sig_cutoff = 4, nominal_cutoff = 2, colors = c('#D41159', '#EA88AC', '#FFFFFF', '#8DC2FF', '#1A85FF')) {
  sig_cutoff = abs(sig_cutoff)
  nominal_cutoff = abs(nominal_cutoff)
  scales = gen_five_scales(sig_cutoff, nominal_cutoff)
  oo = rep(scales[3], length(zval))
  oo[zval > sig_cutoff] = scales[1]
  oo[zval < - sig_cutoff] = scales[5]
  oo[zval <= sig_cutoff & zval > nominal_cutoff] = scales[2]
  oo[zval >= - sig_cutoff & zval < - nominal_cutoff] = scales[4]
  color_scale = colors
  names(color_scale) = scales
  return(list(category = oo, color_code = color_scale))
}

vis_by_tag = function(datadir, tag, df, score) {
  vis_bg = readRDS(paste0(datadir, '/bg_img.rds'))
  vis_meta = readRDS(paste0(datadir, '/meta_plot.rds'))
  if(!tag %in% names(vis_meta)) {
    message('Wrong tag.')
    return(NULL)
  } else {
    vis_meta = vis_meta[[tag]]
  }
  df_sub = df %>% filter(IDP %in% vis_meta$IDP)
  vis = readRDS(paste0(datadir, '/', vis_meta$vis_rds))
  
  sub = vis$df_color %>% filter(IDP %in% df_sub$IDP)
  if(sum(duplicated(sub$color_code)) > 0) {
    message('Wrong IDP set.')
    return(NULL)
  }
  
  img = vis$img
  img2 = img
  img2[] = NA
  for(i in 1 : nrow(sub)) {
    img2[img == sub$color_code[i]] = df_sub[[score]][ df_sub$IDP == sub$IDP[i] ] 
  }
  tt = reshape2::melt(img2) %>% filter(value != 0)
  bb = reshape2::melt(vis_bg) %>% filter(value != 0) %>% rename(bg = value)
  bb = left_join(bb, tt, by = c('Var1', 'Var2', 'Var3'))
  mid1 = vis_bg@dim_[2]
  mid2 = vis_bg@dim_[3]
  mid3 = vis_bg@dim_[4]
  d1 = vis_meta$slide_position[1]
  d2 = vis_meta$slide_position[2]
  d3 = vis_meta$slide_position[3]
  tmp = rbind(
    bb %>% filter(Var1 == floor(mid1 / d1)) %>% mutate(direction = 'a1') %>% rename(x = Var2, y = Var3, ref = Var1),
    bb %>% filter(Var2 == floor(mid2 / d2)) %>% mutate(direction = 'a2') %>% rename(x = Var1, y = Var3, ref = Var2),
    bb %>% filter(Var3 == floor(mid3 / d3)) %>% mutate(direction = 'a3') %>% rename(x = Var1, y = Var2, ref = Var3)
  )
  
  if(score == 'zscore') {
    kk = zscore2category(tmp$value)
    tmp$value_category = kk$category
    p = tmp %>% ggplot() + 
      geom_raster(aes(x, y, fill = value_category)) +
      scale_fill_manual(values = kk$color_code, na.value = 'transparent') + 
      # scale_fill_gradient2(name = score, low = 'blue', mid = 'white', high = 'red', midpoint = 0, na.value = 'transparent') +
      geom_tile(aes(x, y, alpha = -bg), fill = "grey20") +
      scale_alpha(range = c(0.2, 0.8)) +
      coord_equal() + facet_grid(.~direction, labeller = label_both) +
      guides(color = guide_legend("-bg"), alpha = FALSE)
  } else {
    p = tmp %>% ggplot() + 
      geom_raster(aes(x, y, fill = value)) +
      scale_fill_gradient2(name = score, low = 'blue', mid = 'white', high = 'red', midpoint = 0, na.value = 'transparent') +
      geom_tile(aes(x, y, alpha = -bg), fill = "grey20") +
      scale_alpha(range = c(0.2, 0.8)) +
      coord_equal() + facet_grid(.~direction, labeller = label_both) +
      guides(color = guide_legend("-bg"), alpha = FALSE)
  }
  
  p + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
}


