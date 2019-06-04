# script to test whether fixed effect of variable with random slope can be estimated correctly
# when benchmark is a bayesian t-test (bf-package) with the "true" b-values of the subjects
# test whether bf randslope function approaches t-test results with increasing n
# therefore, more values for nsubj in this simulation, but set eff.size.cat to 0

rm(list=ls(all=TRUE))
Sys.setenv(LOCAL_CPP = "-mtune = native")

# required functions/packages
source("modelrun.R")
source("utils.R")
library(BayesFactor)
library(brms)
library(tidyverse)

nsubj <- c(20, 40)
nobs <- 10
pop.sd.cont <- 1
ncont <- 1
xcont <- seq(1,5,1)
xcont.mc <- xcont-mean(xcont)
eff.size.cont <- c(0, 0.2, 0.5, 0.8)
intercept <- 0

# mcmc values BayesRS (brs for BayesRS)
nadapt_brs <- 1000
nburn_brs <- 1000
mcmcstep_brs <- c(10000)
# mcmc values brms
mcmcstep_brms <- 10000
n_chains_brms <- 1#3

# brms specifications
fixefPrior <- c(set_prior("cauchy(0,0.353)", class="b", coef="x.cont"))
ranefPrior <- set_prior("gamma(1,0.04)", class="sd")
control = list(adapt_delta = .9)

# simulation values
nreps <- 2

# results container
bfs <- matrix(nrow = nreps*length(nsubj)*length(eff.size.cont)*length(mcmcstep_brs)*length(nobs), ncol = 9)

count <- 1
dat.str <- data.frame(iv = c("x.cont"), 
                      type = c("cont"),
                      id = c(1))

for(rep in 1:nreps){
  for(n in nsubj){
    for(k in nobs){
      for(ef.cont in eff.size.cont){
        # generate data
        vals <- data.frame(id = rep(seq(1,n,1), each = k*length(xcont.mc)),
                           eff.size.cont = ef.cont,
                           x.cont = xcont.mc,
                           b.cont = rep(rnorm(n, mean = ef.cont, sd = pop.sd.cont), each = k*length(xcont.mc)),
                           error = rnorm(n*k*length(xcont.mc), mean = rep(rnorm(n, mean = intercept), 
                                                                          each = k*length(xcont.mc)), sd = 1))
        vals$y = vals$x.cont*vals$b.cont + vals$error
        dat.model <- vals[, c("id", "y", "x.cont")]
        dat.model$id <- as.factor(dat.model$id)
        
        ## BayesFactor package:
        # extract bf for fixed effect from "true" b values with ttest of bayesfactor package
        tt.bfp <- ttestBF(x = unique(vals$b.cont), mu = 0)
        bf_bfp <- tt.bfp@bayesFactor$bf
        
        ## brms package:
        # bf via savage-dickey method (savage-dickey aka sd)
        m_brms_full_prior <- brm(y ~ x.cont + (1+x.cont|id),
                                 prior=c(fixefPrior, ranefPrior), chains=n_chains_brms, iter=mcmcstep_brms,
                                 data=dat.model, save_all_pars=TRUE, sample_prior = TRUE)
        h_cont <- hypothesis(m_brms_full_prior, "x.cont = 0")
        bf_brms_sd <- 1/h_cont$hypothesis$Evid.Ratio
        
        # bf via bridge sampler
        m_brms_full_no_prior <- brm(y ~ x.cont + (1+x.cont|id),
                                    prior=c(fixefPrior, ranefPrior), chains=n_chains_brms, iter=mcmcstep_brms,
                                    data=dat.model, save_all_pars=TRUE)
        m_brms_no_cont_fixed <- brm(y ~ (1+x.cont|id),
                                    prior=ranefPrior, chains=n_chains_brms, iter=mcmcstep_brms, 
                                    data=dat.model, save_all_pars=TRUE)
        bf_brms_bridge <- bayes_factor(m_brms_full_no_prior, m_brms_no_cont_fixed)
        
        
        ## BayesRS package:
        for(steps in mcmcstep_brs){
          # run model and save bfs
          out <- modelrun(data = dat.model, dv = "y", dat.str = dat.str, 
                          nadapt = nadapt_brs, nburn = nburn_brs, nsteps = steps, 
                          checkconv = 0)
          bf_brs_sd <- as.numeric(as.character(out[[1]]$bf))
          # add all in same df
          bfs[count,] <- c(as.integer(count), n, k, ef.cont, steps, bf_bfp, 
                           log(bf_brms_bridge$bf), log(bf_brms_sd), log(bf_brs_sd))
          print(paste0(count/nrow(bfs)*100, "% of all runs"))
          save(bfs, file = "bfs.Rda")
          count <- count + 1
        }
      }
    }
  }
}

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

plotTheme(ggplot(bfs_long_comp_bfp, aes(bf.bfp, BF_log_diff)) + 
            geom_segment(aes(x = -5, y = 0, xend = 22, yend = 0, colour = "red"), size = 1.5)+
            geom_segment(aes(x = 0, y = v.min, xend = 0, yend = v.max), size = .5) +
            geom_point(aes(shape = Comparison)) +
            facet_wrap(~Comparison) + theme_bw() +
            scale_color_discrete(guide = FALSE) +
            xlab("Bayes Factor BayesFactor Package") + ylab("Log Diff BFs"))

bfs_long_comp_brms_bridge <- bfs.df %>% 
  select(bf.brms.bridge, bf_bfp_vs_brms_bridge, bf_brms_sd_vs_brms_bridge, bf_brs_sd_vs_brms_bridge) %>%
  gather(Comparison, BF_log_diff, c(bf_bfp_vs_brms_bridge, bf_brms_sd_vs_brms_bridge, bf_brs_sd_vs_brms_bridge))

plotTheme(ggplot(bfs_long_comp_brms_bridge, aes(bf.brms.bridge, BF_log_diff)) + 
            geom_segment(aes(x = -5, y = 0, xend = 22, yend = 0, colour = "red"), size = 1.5)+
            geom_segment(aes(x = 0, y = v.min, xend = 0, yend = v.max), size = .5) +
            geom_point(aes(shape = Comparison)) +
            facet_wrap(~Comparison) + theme_bw() +
            scale_color_discrete(guide = FALSE) +
            xlab("Bayes Factor BayesFactor Package") + ylab("Log Diff BFs"))

v.min <- -16
v.max <- 16
plotTheme(ggplot(bfs_long_comp_bfp, aes(bf_brms_bridge, BF_log_diff)) + 
            labs(title = "\n") +
            labs(x="\nLog(True BF)", y="Log(comp. BF)-Log(true BF)\n") + 
            geom_segment(aes(x = -5, y = 0, xend = 22, yend = 0, colour = "red"), size = 1.5)+
            geom_point(shape = 1, size = 3) + facet_wrap(~n,nrow = 2) +
            geom_segment(aes(x = -2.3026, y = v.min, xend = -2.3026, yend = v.max), size = .5) +
            geom_segment(aes(x = 2.3026, y = v.min, xend = 2.3026, yend = v.max), size = .5) +
            geom_segment(aes(x = -1.1, y = v.min, xend = -1.1, yend = v.max), size = 1) +
            geom_segment(aes(x = 1.1, y = v.min, xend = 1.1, yend = v.max), size = 1)) +
  coord_cartesian(ylim = c(v.min, v.max)) + theme(legend.position = "hide") +
  facet_wrap(n ~ Comparison)

# manually inspect outliers
filter <- bfs.df$mcmcsteps == 50000 & bfs.df$bf.cont.true<(-1.1) & bfs.df$bf.cont.true>-999
outliers <- bfs.df[filter,]
outliers$diff <- outliers$bf.cont-outliers$bf.cont.true
true.max <- exp(outliers$bf.cont.true[which(outliers$diff == max(outliers$diff))])
diff.max <- exp(outliers$diff[which(outliers$diff == max(outliers$diff))])
true.min <- exp(outliers$bf.cont.true[which(outliers$diff == min(outliers$diff))])
diff.min <- exp(outliers$diff[which(outliers$diff == min(outliers$diff))])


