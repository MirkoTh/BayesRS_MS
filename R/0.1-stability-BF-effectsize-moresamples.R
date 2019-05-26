# script to test reliability of BF function with data having random slopes on one continuous variable
# and having a fixed categorical effect (2 levels)

rm(list=ls(all=TRUE))

nsubj <- c(20, 40)
nobs <- 10
pop.sd.cont <- 1
ncont <- 1
xcont <- seq(1,5,1)
xcont.mc <- xcont-mean(xcont)
eff.size.cont <- c(0, 0.2, 0.5, 0.8)
intercept <- 0

# mcmc values
nadapt = 1000
nburn = 1000
mcmcstep = c(50000, 100000)


# simulation values
nreps <- 50

# results container
bfs <- matrix(nrow = nreps*length(nsubj)*length(eff.size.cont)*length(mcmcstep), ncol = 5)
source("modelrun.R")
count <- 1
dat.str <- data.frame(iv = c("x.cont"), 
                      type = c("cont"),
                      id = c(1))

require(MASS)

for(ef.cont in eff.size.cont){
  for(n in nsubj){
    # generate data
    vals <- data.frame(id = rep(seq(1,n,1), each = nobs*length(xcont.mc)),
                       eff.size.cont = ef.cont,
                       x.cont = xcont.mc,
                       b.cont = rep(mvrnorm(n, mu = ef.cont, Sigma = pop.sd.cont, empirical = TRUE),
                                    each = nobs*length(xcont.mc)), 
                       error = rnorm(n*nobs*length(xcont.mc),
                                     mean = rep(mvrnorm(n, mu = intercept, Sigma = 1, empirical = TRUE), each = nobs*length(xcont.mc)),
                                     sd = 1))
    vals$y = intercept + vals$x.cont*vals$b.cont + vals$error
    dat.model <- vals[, c("id", "y", "x.cont")]
    # run model several times on same df
    for(steps in mcmcstep){
      for(rep in 1:nreps){
        # run model and save bfs
        out <- modelrun(data = dat.model, dv = "y", dat.str = dat.str, nadapt = nadapt, nburn = nburn, nsteps = steps, 
                        checkconv = 0)
        bf <- out[[1]]
        bfs[count,] <- c(as.integer(count), n, ef.cont, steps, as.numeric(bf))
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

plotTheme(ggplot(vals, aes(x.cont, y)) + stat_summary(fun.y = mean, geom = "line"))


bfs.df <- as.data.frame(bfs)
names(bfs.df) <- c("x", "n", "ef.cont", "mcmcsteps", "bf.cont")
bfs.df$mcmcsteps <- as.factor(bfs.df$mcmcsteps)
levels(bfs.df$mcmcsteps) <- c("Steps = 50'000", "Steps = 100'000")
require(grid)
for (ef.cont in eff.size.cont){
  for(n in nsubj){
    
    dat.efcont <- bfs.df[!is.na(bfs.df$x) & bfs.df$ef.cont == ef.cont & bfs.df$n == n,]
    bf.cont.vals <- as.data.frame(as.list(aggregate(bf.cont ~ mcmcsteps, data = dat.efcont,
                                                    FUN=function(x) c(m = mean(x), s = sd(x), min = min(x), max = max(x)))))
    names(bf.cont.vals) <- c("mcmcsteps", "m", "sd", "min", "max")
    bf.cont.vals$range <- apply(data.frame((bf.cont.vals$max/bf.cont.vals$m-1)*100,
                                           (1-bf.cont.vals$min/bf.cont.vals$m)*100), 
                                1, max)
    SD <- bf.cont.vals$sd
    BFCONT <- dat.efcont$bf.cont
    plot <- plotTheme(ggplot(dat.efcont, aes(bf.cont)) +
                        geom_histogram(colour = "black", fill = "white", bins = 200) +
                        facet_wrap(~ mcmcsteps, nrow = 2) +
                        labs(title = paste("True Effect Size = ", ef.cont, "N = ", n,"\n")) +
                        labs(x=(expression("\n"*BF[10])), y= "Count\n") +
                        geom_segment(aes(x = min(-1*max(SD), min(BFCONT)-1*max(SD)),
                                         y = 0, xend = max(max(BFCONT)+2*max(SD),1.1), yend = 0)))  +
      theme(strip.text.x = element_text(size = 12), legend.position = "hide")
    
    x11(5,4.5)
    grid.draw(plot + geom_segment(aes(x = 1, y = min(ggplot_build(plot)$data[[1]]$count),
                                  xend = 1, yend = max(ggplot_build(plot)$data[[1]]$count)), linetype = "dashed") +
      coord_cartesian(xlim = c(min(-1*max(SD), min(BFCONT)-1*max(SD)),
                               max(max(BFCONT)+3*max(SD),1.1)))+
      annotate("text", x=mean(layer_scales(plot)$x$range$range),
               y=max(ggplot_build(plot)$data[[1]]$count)-.5*max(ggplot_build(plot)$data[[1]]$count), 
               label= c(paste("m =", prettyNum(bf.cont.vals[,2], digits = 2), "; sd =", prettyNum(bf.cont.vals[,3], digits = 2),
                              "\n\u00B1", prettyNum(bf.cont.vals[,6], digits = 2), "%")),
               size = 5))
  }
}
