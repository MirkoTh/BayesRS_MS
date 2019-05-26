# script to compare random-intercept-only model to model with random intercept and random slope 

rm(list=ls(all=TRUE))

nsubj <- c(20, 40, 60)
nobs <- seq(3,9,3)
ncont <- 1
xcont <- seq(1,5,1)
xcont.mc <- xcont-mean(xcont)
eff.size.cont <- seq(0, 0.75, 0.15)
intercept <- 0
psd.cont <- 1

# mcmc values
nadapt <- 1000
nburn <- 1000
mcmcstep <- c(100000)

# simulation values
nreps <- 50

# results container
bfs <- matrix(nrow = nreps*length(nsubj)*length(eff.size.cont)*
                length(mcmcstep)*length(psd.cont), ncol = 7)
source("modelrun.R")
require(BayesFactor)
count <- 1
dat.str <- data.frame(iv = c("x.cont"), 
                      type = c("cont"),
                      id = c(0))
for(rep in 1:nreps){
  for(n in nsubj){
    for(ef.cont in eff.size.cont){
      for(k in nobs){
        # generate data
        vals <- data.frame(id = rep(seq(1,n,1), each = k*length(xcont.mc)),
                           eff.size.cont = ef.cont,
                           x.cont = xcont.mc,
                           b.cont = rep(rnorm(n, mean = ef.cont, sd = psd.cont), each = k*length(xcont.mc)),
                           error = rnorm(n*k*length(xcont.mc),
                                         mean = rep(rnorm(n, mean = intercept), each = k*length(xcont.mc)),
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
        # run model and save bfs
        out <- modelrun(data = dat.model, dv = "y", dat.str = dat.str, nadapt = nadapt, nburn = nburn, nsteps = mcmcstep, 
                        checkconv = 0)
        bf <- out[[1]]
        bfs[count,] <- c(as.integer(count), n, k, ef.cont, mcmcstep, as.numeric(bf), bf.corr.cont)
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

#### post proc ####

bfs.df <- as.data.frame(bfs)
names(bfs.df) <- c("x", "n", "nobs", "ef.cont", "mcmcsteps", "bf.cont", "bf.cont.true")
df.cont <- bfs.df[, c("n", "ef.cont", "mcmcsteps", "bf.cont","bf.cont.true")]
df.cont.l <- reshape(df.cont, varying = c("bf.cont", "bf.cont.true"),
                     v.names = "BF", times = c("New R Package", "BayesFactor Package"), direction = 'long')
names(df.cont.l)[4] <- "Method"
df.cont.l$Method <- as.factor(df.cont.l$Method)
plotTheme(ggplot(df.cont.l, aes(ef.cont, BF)) + 
            labs(title = "Comparison with fixed-effect model from BFP\n") +
            labs(x="\nEffect Size", y="Log(BF)") + theme_bw() + theme(plot.title = element_text(size = rel(1.75))) +
            facet_wrap(~Method+n, scales = "free") +
            geom_segment(aes(x = -0.1, y = 0, xend = .8, yend = 0), colour = "red", size = 1)+
            geom_point(size = 3, position = position_dodge(.05), shape = 1)) +
  scale_shape(solid = FALSE) + coord_cartesian(ylim=c(-5,40))

df.cont.l$BF <- exp(df.cont.l$BF)
df.cont.l$BF.dich <- 0
# df.cont.l$BF.dich[df.cont.l$BF< 0.3125] <- "Substantial Null"
# df.cont.l$BF.dich[df.cont.l$BF>=(0.3125)&df.cont.l$BF<=(3.2)] <- "Ambiguous"
# df.cont.l$BF.dich[df.cont.l$BF> (3.2)] <- "Substantial Alternative"

df.cont.l$BF.dich[df.cont.l$BF< 0.1] <- "Strong Null"
df.cont.l$BF.dich[df.cont.l$BF>=(0.1)&df.cont.l$BF<=(10)] <- "Ambiguous"
df.cont.l$BF.dich[df.cont.l$BF> 10] <- "Strong Alternative"

table(df.cont.l$BF.dich, df.cont.l$ef.cont, df.cont.l$Method, df.cont.l$n)/100

bf.table <- aggregate(as.factor(BF.dich) ~ ef.cont + Method + n, data = df.cont.l, FUN = table)
names(bf.table) <- c("ef.cont", "Method", "n", "BF")
bf.table2 <- as.data.frame(as.matrix(bf.table))
vars <- names(bf.table2[,c(4:6)])
bf.table2[,vars] <- apply(bf.table2[,vars],2,FUN = function(x) as.numeric(as.character(x)))
bf.table2[,vars] <- bf.table2[,vars]/100
bf.dich.long <- reshape(bf.table2, varying = vars, v.names = "Proportion",
                        times = c("Ambiguous", "Strong Alternative", "Strong Null"), direction = 'long')
names(bf.dich.long)[4] <- "Evidence"

plotTheme(ggplot(bf.dich.long, aes(ef.cont, Proportion, group = Evidence)) +
            geom_line(aes(color = Evidence)) + 
            geom_point(aes(color = Evidence), size = 7, colour = "white") + 
            geom_point(aes(color = Evidence), size = 3, shape = 1) + 
            facet_wrap(~Method+n,nrow=2) + labs(y="Proportion\n", x="\nEffect Size") +
            scale_colour_manual(values = c("red", "blue", "green"))) +
  theme(legend.position = "bottom")

aggregate(Proportion ~ Method + Evidence + ef.cont, data = bf.dich.long, FUN = mean)

# same with five BF categorization categories
df.cont.l$BF.dich <- 0
df.cont.l$BF.dich[df.cont.l$BF< 0.1] <- "Strong Null"
df.cont.l$BF.dich[df.cont.l$BF>=(0.1)&df.cont.l$BF<(0.3125)] <- "Substantial Null"
df.cont.l$BF.dich[df.cont.l$BF>=(0.3125)&df.cont.l$BF<=(3.2)] <- "Ambiguous"
df.cont.l$BF.dich[df.cont.l$BF>(3.2)&df.cont.l$BF<=(10)] <- "Substantial Alternative"
df.cont.l$BF.dich[df.cont.l$BF> (10)] <- "Strong Alternative"
table(df.cont.l$BF.dich, df.cont.l$ef.cont, df.cont.l$Method, df.cont.l$n)/50

bf.table <- aggregate(as.factor(BF.dich) ~ ef.cont + Method + n, data = df.cont.l, FUN = table)
names(bf.table) <- c("ef.cont", "Method", "n", "BF")
bf.table2 <- as.data.frame(as.matrix(bf.table))
vars <- names(bf.table2[,c(4:8)])
bf.table2[,vars] <- apply(bf.table2[,vars],2,FUN = function(x) as.numeric(as.character(x)))
bf.table2[,vars] <- bf.table2[,vars]/100
bf.dich.long <- reshape(bf.table2, varying = vars, v.names = "Proportion",
                        times = c("Ambiguous", "Strong Alternative", "Strong Null", 
                                  "Substantial Alternative", "Substantial Null"), direction = 'long')

names(bf.dich.long)[4] <- "Evidence"
bf.dich.long$Evidence <- as.factor(bf.dich.long$Evidence)
bf.dich.long$Evidence <- factor(bf.dich.long$Evidence, levels = levels(bf.dich.long$Evidence)[c(3,5,1,4,2)])

plotTheme(ggplot(bf.dich.long, aes(ef.cont, Proportion, group = Evidence)) +
            geom_line(aes(color = Evidence, linetype = Evidence)) + 
            geom_point(aes(color = Evidence), size = 7, colour = "white") + 
            geom_point(aes(color = Evidence), size = 3, shape = 1) + 
            facet_wrap(~Method+n,nrow=2) + labs(title="\n",y="Proportion\n", x="\nEffect Size") +
            scale_color_manual(values = c("red", "orange","grey","light blue", "dark blue")) +
            scale_linetype_manual(values = c("solid", "solid","dashed","solid", "solid"))) +
  theme(legend.position = "bottom")+guides(color=guide_legend(nrow=2,byrow=TRUE))

aggregate(Proportion ~ Method + Evidence + ef.cont, data = bf.dich.long, FUN = mean)

es.0<-df.cont.l[df.cont.l$ef.cont==0,]
plotTheme(ggplot(es.0,aes(log(BF)))+geom_histogram(color="black",fill="white")+facet_grid(~Method, scales = "free"))
plotTheme(ggplot(es.0[log(es.0$BF)<75,],aes(log(BF)))+
            geom_histogram(color="black",fill="white") + facet_grid(~Method, scales = "free")+
            geom_segment(aes(x = -2.3026, y = -2, xend = -2.3026, yend = 90), size = .5) +
            geom_segment(aes(x = 2.3026, y = -2, xend = 2.3026, yend = 90), size = .5) +
            geom_segment(aes(x = -1.1, y = -2, xend = -1.1, yend = 90), size = 1) +
            geom_segment(aes(x = 1.1, y = -2, xend = 1.1, yend = 90), size = 1) + 
            labs(title="\n",y="Count\n", x="\nlog(BF"))


