# script to compare the BF from a fixed effect model computed with
# BayesRS package to the BF computed with the BayesFactor package

rm(list=ls(all=TRUE))

nsubj <- c(20, 40)
nobs <- 10
pop.sd.cont <- 1
ncont <- 1
xcont <- seq(1,5,1)
xcont.mc <- xcont-mean(xcont)
eff.size.cont <- c(.05,.1,.3)
intercept <- 0

dat.str <- data.frame(iv = c("x.cont"), 
                      type = c("cont"),
                      id = c(0))

# mcmc values
nadapt <- 1000
nburn <- 1000
mcmcstep <- c(50000, 100000)

# simulation values
nreps <- 50

# results container
bfs <- matrix(nrow = nreps*length(nsubj)*length(eff.size.cont)*length(mcmcstep), ncol = 6)
source("modelrun.R")
require(BayesFactor)
count <- 1

for(rep in 1:nreps){
  for(n in nsubj){
      for(ef.cont in eff.size.cont){
        
        # generate data
        vals <- data.frame(id = rep(seq(1,n,1), each = nobs*length(xcont.mc)),
                           eff.size.cont = ef.cont,
                           x.cont = xcont.mc,
                           b.cont = ef.cont,
                           error = rnorm(n*nobs*length(xcont.mc),
                                         mean = rep(rnorm(n, mean = intercept), each = nobs*length(xcont.mc)),
                                         sd = 1))
        vals$y = vals$x.cont*vals$b.cont + vals$error
        dat.model <- vals[, c("id", "y", "x.cont")]
        
        dat.model$id <- as.factor(dat.model$id)
        mf <- generalTestBF(y ~ id + x.cont, data = dat.model, whichRandom = "id",
                   iterations=100000, progress=T, neverExclude = "id")
        bf.corr.cont <- -999
        # extract bf for cont eff
        # n.b. this extracs the log of the bfs
        bf.corr.cont <- mf[1]@bayesFactor$bf-mf[2]@bayesFactor$bf

        for(steps in mcmcstep){
          # run model and save bfs

          out <- modelrun(data = dat.model, dv = "y", dat.str = dat.str, nadapt = nadapt, nburn = nburn, nsteps = steps, 
                          checkconv = 0)
          bf <- out[[1]]
          bfs[count,] <- c(as.integer(count), n, ef.cont,steps, as.numeric(bf), bf.corr.cont)
          print(count/nrow(bfs)*100)
          name <- paste("bfs", count, ".Rda", sep = "")
          save(bfs, file = name)
          count <- count + 1
        }
      }
  }
}

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

bfs.df <- as.data.frame(bfs)
names(bfs.df) <- c("x", "n", "ef.cont", "mcmcsteps", "bf.cont", "bf.cont.true")
bfs.df$n <- as.factor(bfs.df$n)
bfs.df$mcmcsteps <- as.factor(bfs.df$mcmcsteps)
levels(bfs.df$n) <- c("N = 20", "N = 40")
levels(bfs.df$mcmcsteps) <- c("Steps = 50'000", "Steps = 100'000")
x.max <- 220
y.max <- 25
y.min <- -25
plotTheme(ggplot(bfs.df, aes(bf.cont.true, bf.cont-bf.cont.true)) + 
            labs(title = "\n", x="\nLog(True BF)", y="Log(comp. BF)-Log(true BF)\n") +
            geom_point(shape = 1, size = 3) + facet_wrap(~mcmcsteps+n) +
            geom_segment(aes(x = -10, y = 0, xend = x.max, yend = 0, colour = "red"), size = 1.5)) +
  coord_cartesian(ylim = c(y.min, y.max)) + theme(legend.position = "hide") 


plotTheme(ggplot(bfs.df[bfs.df$bf.cont.true>-999&bfs.df$bf.cont.true<2.5,], aes(bf.cont.true, bf.cont-bf.cont.true)) + 
            labs(title = "\n") +
            labs(x="\nLog(True BF)", y="Log(comp. BF)-Log(true BF)\n") +
            geom_point(shape = 1, size = 3) + facet_wrap(~mcmcsteps) +
            geom_segment(aes(x = -3, y = 0, xend = 2.5, yend = 0, colour = "red"), size = 1.5)) +
  coord_cartesian(ylim = c(-2, 2)) + theme(legend.position = "hide") + 
  geom_segment(aes(x = -2.3026, y = -6, xend = -2.3026, yend = 12), size = .5) +
  geom_segment(aes(x = 2.3026, y = -6, xend = 2.3026, yend = 12), size = .5) +
  geom_segment(aes(x = -1.1, y = -6, xend = -1.1, yend = 12), size = 1) +
  geom_segment(aes(x = 1.1, y = -6, xend = 1.1, yend = 12), size = 1) +
  theme(strip.text = element_text(size = 12))


idx <- which((bfs.df$bf.cont-bfs.df$bf.cont.true)==min((bfs.df$bf.cont-bfs.df$bf.cont.true)[bfs.df$ef.cont==0.05]))
min.diff <- (bfs.df$bf.cont-bfs.df$bf.cont.true)[idx]
act.bf <- exp(bfs.df$bf.cont.true[idx])
est.bf <- exp(bfs.df$bf.cont[idx])

filter<-bfs.df$bf.cont.true<-1/0.3516795
# manually inspect outliers
outliers <- bfs.df[filter,]
outliers$diff <- outliers$bf.cont-outliers$bf.cont.true
true.max <- exp(outliers$bf.cont.true[which(outliers$diff == max(outliers$diff))])
diff.max <- exp(outliers$diff[which(outliers$diff == max(outliers$diff))])
true.min <- exp(outliers$bf.cont.true[which(outliers$diff == min(outliers$diff))])
diff.min <- exp(outliers$diff[which(outliers$diff == min(outliers$diff))])

min(bfs.df$bf.cont[filter]-bfs.df$bf.cont.true[filter])
