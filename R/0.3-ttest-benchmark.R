# script to test whether fixed effect of variable with random slope can be estimated correctly
# when benchmark is a bayesian t-test (bf-package) with the "true" b-values of the subjects
# test whether bf randslope function approaches t-test results with increasing n
# therefore, more values for nsubj in this simulation, but set eff.size.cat to 0

rm(list=ls(all=TRUE))

# required functions/packages
source("modelrun.R")
library(BayesFactor)
library(brms)

nsubj <- 10 #c(10,25,50,100)
nobs <- 3 #seq(3,9,3)
pop.sd.cont <- 1
ncont <- 1
xcont <- seq(1,5,1)
xcont.mc <- xcont-mean(xcont)
eff.size.cont <- .2 #c(0, 0.2, 0.5, 0.8)
intercept <- 0

# mcmc values BayesRS
nadapt_rs <- 1000
nburn_rs <- 1000
mcmcstep_rs <- c(100000)
# mcmc values brms
mcmcstep_brms <- 10000
n_chains_brms <- 3

# simulation values
nreps <- 2

# results container
bfs <- matrix(nrow = nreps*length(nsubj)*length(eff.size.cont)*length(mcmcstep_rs)*length(nobs), ncol = 8)

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
        corr.bf <- tt.bfp@bayesFactor$bf
        
        
        
        ## brms package:
        fixefPrior <- c(set_prior("cauchy(0,0.353)", class="b", coef="x.cont"))
        ranefPrior <- set_prior("gamma(1,0.04)", class="sd")
        m_brms <- brm(y ~ x.cont + (1+x.cont|id),
                               prior=c(fixefPrior, ranefPrior), chains=n_chains_brms, iter=mcmcstep_brms, 
                      data=dat.model, save_all_pars=TRUE, sample_prior = TRUE)
        h_cont <- hypothesis(m_brms, "x.cont = 0")
        bf_brms <- h_cont$hypothesis$Evid.Ratio
        
        ## BayesRS package:
        for(steps in mcmcstep_rs){
          # run model and save bfs
          out <- modelrun(data = dat.model, dv = "y", dat.str = dat.str, 
                          nadapt = nadapt_rs, nburn = nburn_rs, nsteps = steps, 
                          checkconv = 0)
          bf <- as.numeric(as.character(out[[1]]$bf))
          bfs[count,] <- c(as.integer(count), n, k, ef.cont, steps, log(bf), bf_brms, corr.bf)
          print(count/nrow(bfs)*100)
          save(bfs, file = "bfs.Rda")
          count <- count + 1
        }
      }
    }
  }
}

bfs.df <- as.data.frame(bfs)
names(bfs.df) <- c("x", "n",  "nobs", "ef.cont","mcmcsteps", "bf.cont", "bf.cont.brms", "bf.cont.true")

require(ggplot2)

# plotting scheme ####
plotTheme <- function (plot){
  plot +
    theme(panel.background = element_rect(fill = "black")) +
    theme(panel.grid=element_blank(),panel.background=element_blank()) + theme_bw() + 
    theme(plot.title = element_text(size = rel(1.75), hjust = 0)) + theme(axis.line = element_line(color ='black')) +
    theme(legend.title = element_text(size=13))+
    theme(legend.text = element_text(size = 13), legend.key = element_blank(), legend.box = 'horizontal')+
    theme(axis.title.x = element_text(size=15),axis.text.x  = element_text(size=13))+  
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(size=13))+ 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.text = element_text(size = 12))
}

#### post proc ####
bfs.df$n[bfs.df$n==20] <- "N = 20"
bfs.df$n[bfs.df$n==40] <- "N = 40"
bfs.df$n <- as.factor(bfs.df$n)

# bfs.df$bf.cont.true<10&
v.min <- -16
v.max <- 16
plotTheme(ggplot(bfs.df, aes(bf.cont.true, bf.cont-bf.cont.true)) + 
            labs(title = "\n") +
            labs(x="\nLog(True BF)", y="Log(comp. BF)-Log(true BF)\n") + 
            geom_segment(aes(x = -5, y = 0, xend = 22, yend = 0, colour = "red"), size = 1.5)+
            geom_point(shape = 1, size = 3) + facet_wrap(~n,nrow = 2) +
            geom_segment(aes(x = -2.3026, y = v.min, xend = -2.3026, yend = v.max), size = .5) +
            geom_segment(aes(x = 2.3026, y = v.min, xend = 2.3026, yend = v.max), size = .5) +
            geom_segment(aes(x = -1.1, y = v.min, xend = -1.1, yend = v.max), size = 1) +
            geom_segment(aes(x = 1.1, y = v.min, xend = 1.1, yend = v.max), size = 1)) +
  coord_cartesian(ylim = c(v.min, v.max)) + theme(legend.position = "hide")

# manually inspect outliers
filter <- bfs.df$mcmcsteps == 50000 & bfs.df$bf.cont.true<(-1.1) & bfs.df$bf.cont.true>-999
outliers <- bfs.df[filter,]
outliers$diff <- outliers$bf.cont-outliers$bf.cont.true
true.max <- exp(outliers$bf.cont.true[which(outliers$diff == max(outliers$diff))])
diff.max <- exp(outliers$diff[which(outliers$diff == max(outliers$diff))])
true.min <- exp(outliers$bf.cont.true[which(outliers$diff == min(outliers$diff))])
diff.min <- exp(outliers$diff[which(outliers$diff == min(outliers$diff))])


