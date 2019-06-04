load("bfs.Rda")

library(tidyverse)

bfs.df <- as.data.frame(bfs)
names(bfs.df) <- c("x", "n",  "nobs", "ef.cont","mcmcsteps", "bf.bfp",
                   "bf.brms.bridge", "bf.brms.sd", "bf.brs.sd")

#### post proc ####
bfs.df$n[bfs.df$n==20] <- "N = 20"
bfs.df$n[bfs.df$n==40] <- "N = 40"
bfs.df$n <- as.factor(bfs.df$n)


# all against bfp
bfs.df$bf_brms_bridge_vs_bfp <- bfs.df$bf.brms.bridge - bfs.df$bf.bfp
bfs.df$bf_brms_sd_vs_bfp <- bfs.df$bf.brms.sd - bfs.df$bf.bfp
bfs.df$bf_brs_sd_vs_bfp <- bfs.df$bf.brs.sd - bfs.df$bf.bfp
# all against brms bridge
bfs.df$bf_bfp_vs_brms_bridge <- bfs.df$bf.bfp - bfs.df$bf.brms.bridge
bfs.df$bf_brms_sd_vs_brms_bridge <- bfs.df$bf.brms.sd - bfs.df$bf.brms.bridge
bfs.df$bf_brs_sd_vs_brms_bridge <- bfs.df$bf.brs.sd - bfs.df$bf.brms.bridge
# compare two sd methods
bfs.df$bf_brms_sd_vs_brs_sd <- bfs.df$bf.brms.sd - bfs.df$bf.brs.sd

bfs_long_comp_bfp <- bfs.df %>% 
  select(bf.bfp, bf_brms_bridge_vs_bfp, bf_brms_sd_vs_bfp, bf_brs_sd_vs_bfp) %>%
  gather(Comparison, BF_log_diff, bf_brms_bridge_vs_bfp:bf_brs_sd_vs_bfp)

v.min <- -3
v.max <- 15

plotTheme(ggplot(bfs_long_comp_bfp, aes(bf.bfp, BF_log_diff)) + 
            geom_segment(aes(x = -5, y = 0, xend = 22, yend = 0, colour = "red"), size = 1.5)+
            geom_segment(aes(x = 0, y = v.min, xend = 0, yend = v.max), size = .5) +
            geom_point(aes(shape = Comparison), size = 3) +
            facet_wrap(~Comparison) + theme_bw() +
            scale_color_discrete(guide = FALSE) +
            xlab("Bayes Factor BayesFactor Package") + ylab("Log Diff BFs"))

bfs_long_comp_brms_bridge <- bfs.df %>% 
  select(bf.brms.bridge, bf_bfp_vs_brms_bridge, bf_brms_sd_vs_brms_bridge, bf_brs_sd_vs_brms_bridge) %>%
  gather(Comparison, BF_log_diff, c(bf_bfp_vs_brms_bridge, bf_brms_sd_vs_brms_bridge, bf_brs_sd_vs_brms_bridge))

plotTheme(ggplot(bfs_long_comp_brms_bridge, aes(bf.brms.bridge, BF_log_diff)) + 
            geom_segment(aes(x = -5, y = 0, xend = 22, yend = 0, colour = "red"), size = 1.5)+
            geom_segment(aes(x = 0, y = v.min, xend = 0, yend = v.max), size = .5) +
            geom_point(aes(shape = Comparison), size = 3) +
            facet_wrap(~Comparison) + theme_bw() +
            scale_color_discrete(guide = FALSE) +
            xlab("Bayes Factor Bridge Sampling brms Package") + ylab("Log Diff BFs"))



# manually inspect outliers
filter <- bfs.df$mcmcsteps == 50000 & bfs.df$bf.cont.true<(-1.1) & bfs.df$bf.cont.true>-999
outliers <- bfs.df[filter,]
outliers$diff <- outliers$bf.cont-outliers$bf.cont.true
true.max <- exp(outliers$bf.cont.true[which(outliers$diff == max(outliers$diff))])
diff.max <- exp(outliers$diff[which(outliers$diff == max(outliers$diff))])
true.min <- exp(outliers$bf.cont.true[which(outliers$diff == min(outliers$diff))])
diff.min <- exp(outliers$diff[which(outliers$diff == min(outliers$diff))])


