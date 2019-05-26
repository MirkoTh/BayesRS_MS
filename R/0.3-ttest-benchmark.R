# script to test whether fixed effect of variable with random slope can be estimated correctly
# when benchmark is a bayesian t-test (bf-package) with the "true" b-values of the subjects
# test whether bf randslope function approaches t-test results with increasing n
# therefore, more values for nsubj in this simulation, but set eff.size.cat to 0

rm(list=ls(all=TRUE))

nsubj <- c(10,25,50,100)
nobs <- seq(3,9,3)
pop.sd.cont <- 1
ncont <- 1
xcont <- seq(1,5,1)
xcont.mc <- xcont-mean(xcont)
eff.size.cont <- c(0, 0.2, 0.5, 0.8)
intercept <- 0

# mcmc values
nadapt <- 1000
nburn <- 1000
mcmcstep <- c(100000)

# simulation values
nreps <- 50

# results container
bfs <- matrix(nrow = nreps*length(nsubj)*length(eff.size.cont)*length(mcmcstep)*length(nobs), ncol = 7)

source("modelrun.R")
require(BayesFactor)
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
        # extract bf for fixed effect from "true" b values with ttest of bayesfactor package
        tt.bfp <- ttestBF(x = unique(vals$b.cont), mu = 0)
        corr.bf <- tt.bfp@bayesFactor$bf
        
        for(steps in mcmcstep){
          # run model and save bfs
          out <- modelrun(data = dat.model, dv = "y", dat.str = dat.str, nadapt = nadapt, nburn = nburn, nsteps = steps, 
                          checkconv = 0)
          bf <- out[[1]]
          bfs[count,] <- c(as.integer(count), n, k, ef.cont, steps, as.numeric(bf), corr.bf)
          print(count/nrow(bfs)*100)
          name <- paste("bfs", count, ".Rda", sep = "")
          save(bfs, file = name)
          count <- count + 1
        }
      }
    }
  }
}

bfs.df <- as.data.frame(bfs)
names(bfs.df) <- c("x", "n",  "nobs", "ef.cont","mcmcsteps", "bf.cont", "bf.cont.true")

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


