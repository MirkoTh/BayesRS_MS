# test a model that accounts for correlations between random slopes to a model that does not

rm(list=ls(all=TRUE))
#### simulation and mcmc parameters ####
# mcmc values
nadapt <- 1000
nburn <- 1000
nsteps <- 100000

# simulation values
nreps <- 100

#### data info for model ####
# first enter continuous variables, then categorical ones
dat.str <- data.frame(iv = c("x.cont1", "x.cat1"), 
                      type = c("cont", "cat"),
                      subject = c(1,1))
ias.subject <- matrix(0, nrow=nrow(dat.str), ncol = nrow(dat.str))
ias.subject[c(2)] <- 1
randvar.ia <- list(ias.subject)

cor.subject <- matrix(0, nrow=nrow(dat.str)+1, ncol = nrow(dat.str)+1)
cor.subject[c(2,3,6)] <- 1
corstr <- list(cor.subject)

##### simulate data for experiment ####
# expt. info
N <- c(20,100)
nobs <- 8
ntrials <- nobs * 5 * 2

# b means
# main effect cont1, cat1, cont1 x cat1 interaction
grand.mean <- 0
m.cont1 <- 0
m.cat1 <- c(.2, 1)
m.ia1 <- 0
corr <- c(.2,.5,.8)

# results container
bfs.with <- matrix(nrow = nreps*length(m.cat1)*length(corr)*length(N), ncol = 7)
bfs.without <- matrix(nrow = nreps*length(m.cat1)*length(corr)*length(N), ncol = 7)

require(MASS)
require(ggplot2)
require(grid)

count <- 1
# run simulation
for (rep in 1:nreps){
  for (corrsize in corr){
    for (catsize in m.cat1){
      for (s.size in N){
        # correlation between by-subject continuous1 slope, and by-subject categorical1 slope
        subj.cont1.cat1.corr <- mvrnorm(n = s.size, mu = c(m.cont1,catsize), 
                                        Sigma = matrix(data = c(1,corrsize,corrsize,1), nrow = 2, ncol = 2, 
                                                       byrow = TRUE), empirical = FALSE)
        subj.cont1.cat1.corr
        cor(subj.cont1.cat1.corr)
        colMeans(subj.cont1.cat1.corr)
        
        b.cont1.subj <- data.frame(subject = 1:s.size, vals = subj.cont1.cat1.corr[,1])
        b.cat1.subj <- data.frame(subject = 1:s.size, vals = subj.cont1.cat1.corr[,2])
        b.subj.rand <- data.frame(subject = 1:s.size, vals = rnorm(n = s.size, mean = 0, sd = 1))
        b.ia1.subj <- data.frame(subject = 1:s.size, vals = rnorm(n = s.size, mean = 0, sd = 1))
        
        # generate according to lin reg formula
        data.pre2 <- data.frame(subject = rep(1:s.size, each = ntrials), x.cont1 = rep(seq(-2,2), each = ntrials/5), 
                                x.cat1= rep(c(-0.5,0.5), each = ntrials/10))
        data.pre2$y <- 0
        for (i in 1:nrow(data.pre2)){
          data.pre2$y[i] <- b.subj.rand$vals[data.pre2$subject[i]==b.subj.rand$subject] +
            data.pre2$x.cont1[i] * (m.cont1+b.cont1.subj$vals[data.pre2$subject[i]==b.cont1.subj$subject]) +
            data.pre2$x.cat1[i] * (catsize+b.cat1.subj$vals[data.pre2$subject[i]==b.cat1.subj$subject]) +
            data.pre2$x.cont1[i] * data.pre2$x.cat1[i] * (m.ia1+b.ia1.subj$vals[data.pre2$subj[i]==b.ia1.subj$subject])
        }
        # add measurement error
        data.pre2$y <- data.pre2$y + rnorm(n = nrow(data.pre2), mean = 0, sd = 1)
        
        recvars <- which(names(data.pre2) %in% c("subject", "item", "x.cat1"))
        data.pre2[,recvars] <- lapply(data.pre2[,recvars], as.factor)
        
        # data ok?
        
        # x11()
        # pl <- ggplot(data.pre2, aes(x.cont1, y, group = x.cat1)) +
        #   stat_summary(fun.y = mean, geom = "line", aes(colour = x.cat1)) #+
        # facet_wrap(~ subject, nrow = 5)
        # grid.draw(pl)
        
        
        #### run model with correlations ####
        
        source("modelrun.R")
        out <- modelrun(data = data.pre2, dv = "y", dat.str = dat.str, randvar.ia = randvar.ia, nadapt = nadapt, nburn = nburn, 
                        nsteps = nsteps, checkconv = 0, mcmc.save.indiv = 1, corstr = corstr)
        bfs <- out[[1]]
        bfs.with[count,] <- c(as.integer(count), corrsize, catsize, s.size, as.numeric(bfs))
        print(count/nrow(bfs.with)*100)
        name <- paste("bfs.with", count, ".Rda", sep = "")
        save(bfs.with, file = name)
        
        
        #### run model without correlations ####

        out <- modelrun(data = data.pre2, dv = "y", dat.str = dat.str, randvar.ia = randvar.ia, nadapt = nadapt, nburn = nburn, 
                        nsteps = nsteps, checkconv = 0, mcmc.save.indiv = 1)
        bfs <- out[[1]]
        bfs.without[count,] <- c(as.integer(count), corrsize, catsize, s.size, as.numeric(bfs))
        print(count/nrow(bfs.with)*100)
        name <- paste("bfs.without", count, ".Rda", sep = "")
        save(bfs.without, file = name)
        count <- count + 1
      }
    }
  }
}

# post proc and plotting
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
bfs.with<-as.data.frame(bfs.with)
bfs.with$correlation<-"With Correlation"
bfs.without<-as.data.frame(bfs.without)
bfs.without$correlation<-"Without Correlation"
# plotting
bfs.df <- rbind(bfs.with,bfs.without)
names(bfs.df) <- c("x", "Correlation", "Effect.Cat.", "N", "BF.Cont.", "BF.Cat.", "BF.IA", "Method")
df.cont <- bfs.df[, c("N", "Effect.Cat.", "Correlation", "BF.Cont.", "BF.Cat.", "BF.IA", "Method")]
df.cont.l <- reshape(df.cont, varying = c("BF.Cat.", "BF.Cont.", "BF.IA"),
                     v.names = "BF", times = c("Cat.", "Cont.", "IA"), direction = 'long')
names(df.cont.l)[5] <- "Type.of.Effect"
df.cont.l$Method <- as.factor(df.cont.l$Method)
# same with five BF categorization categories
df.cont.l$BF.dich <- 0
df.cont.l$BF.dich[df.cont.l$BF< 0.1] <- "Strong Null"
df.cont.l$BF.dich[df.cont.l$BF>=(0.1)&df.cont.l$BF<(0.3125)] <- "Substantial Null"
df.cont.l$BF.dich[df.cont.l$BF>=(0.3125)&df.cont.l$BF<=(3.2)] <- "Ambiguous"
df.cont.l$BF.dich[df.cont.l$BF>(3.2)&df.cont.l$BF<=(10)] <- "Substantial Alternative"
df.cont.l$BF.dich[df.cont.l$BF> (10)] <- "Strong Alternative"
df.cont.l$N <- as.factor(df.cont.l$N)
levels(df.cont.l$N)<-c("N = 20","N = 100")
df.cont.l$Effect.Cat. <- as.factor(df.cont.l$Effect.Cat.)
levels(df.cont.l$Effect.Cat.)<-c("Effect Cat. = 0","Effect Cat. = 0.2", "Effect Cat. = 1")

## cont
# histogram
plotTheme(ggplot(df.cont.l[df.cont.l$Type.of.Effect=="Cont.",], aes(log(BF))) + 
            geom_histogram(fill="white", color="black") + 
            facet_wrap(~N+Method,nrow=2) +
            labs(title="\n",y="Count\n", x="\nlog(BF)")) + 
  geom_segment(aes(x = -2.3026, y = 0, xend = -2.3026, yend = 320), size = .5) +
  geom_segment(aes(x = 2.3026, y = 0, xend = 2.3026, yend = 320), size = .5) +
  geom_segment(aes(x = -1.1, y = 0, xend = -1.1, yend = 320), size = 1) +
  geom_segment(aes(x = 1.1, y = 0, xend = 1.1, yend = 320), size = 1)

## cat
# density plot
d.pl<-plotTheme(ggplot(df.cont.l[df.cont.l$Type.of.Effect=="Cat.",], aes(log(BF),group=Method,..count..)) + 
            geom_density(color="black", aes(linetype=Method)) + 
            facet_wrap(~N+Effect.Cat.,nrow=2, scales = "free") +
            labs(title="\n",y="Density\n", x="\nlog(BF)"))
tmp<-ggplot_build(d.pl)$data[[1]]
max.panel<-data.frame(N=c("N = 20", "N = 20","N = 20","N = 100","N = 100","N = 100"),
                      Effect.Cat.=rep(c("Effect Cat. = 0","Effect Cat. = 0.2", "Effect Cat. = 1"),2))
max.panel<-cbind(max.panel,aggregate(count ~ PANEL, data = ggplot_build(d.pl)$data[[1]], FUN = max))

 d.pl <- d.pl + 
  geom_segment(data=max.panel, aes(x = -2.3026, y = 0, xend = -2.3026, 
                   yend = max.panel$count), size = .5,inherit.aes = FALSE) +
  geom_segment(data=max.panel, aes(x = 2.3026, y = 0, xend = 2.3026, 
                   yend = max.panel$count), size = .5,inherit.aes = FALSE) +
  geom_segment(data=max.panel, aes(x = -1.1, y = 0, xend = -1.1, 
                   yend = max.panel$count), size = 1,inherit.aes = FALSE) +
  geom_segment(data=max.panel, aes(x = 1.1, y = 0, xend = 1.1, 
                   yend = max.panel$count), size = 1,inherit.aes = FALSE) +
   theme(legend.position = "bottom")+guides(color=guide_legend(nrow=2,byrow=TRUE))
 d.pl
 
 pdf("Corrplot.pdf",  width = 12, height = 9)
 d.pl
 dev.off()

 # aggregated plots
bf.table <- aggregate(as.factor(BF.dich) ~ N + Effect.Cat. + Method + Type.of.Effect, data = df.cont.l, FUN = table)
names(bf.table) <- c("N","Effect.Cat.","Method","Type.of.Effect", "BF")
bf.table2 <- as.data.frame(as.matrix(bf.table))
vars <- names(bf.table2[,c(5:9)])
bf.table2[,vars] <- apply(bf.table2[,vars],2,FUN = function(x) as.numeric(as.character(x)))
bf.table2[,vars] <- bf.table2[,vars]/300
bf.dich.long <- reshape(bf.table2, varying = vars, v.names = "Proportion",
                        times = c("Ambiguous", "Strong Alternative", "Strong Null", 
                                  "Substantial Alternative", "Substantial Null"), direction = 'long')

names(bf.dich.long)[5] <- "Evidence"
bf.dich.long$Evidence <- as.factor(bf.dich.long$Evidence)
bf.dich.long$Evidence <- factor(bf.dich.long$Evidence, levels = levels(bf.dich.long$Evidence)[c(3,5,1,4,2)])

bf.dich.cont<-bf.dich.long[bf.dich.long$Type.of.Effect=="Cont.",]
plotTheme(ggplot(bf.dich.cont, aes(Method, Proportion, group = Evidence)) +
            geom_line(aes(color = Evidence, linetype = Evidence)) + 
            geom_point(aes(color = Evidence), size = 7, colour = "white") + 
            geom_point(aes(color = Evidence), size = 3, shape = 1) + 
            facet_wrap(~N+Effect.Cat.,nrow=2) + labs(title="\n",y="Proportion\n", x="\nEffect Size") +
            scale_color_manual(values = c("red", "orange","grey","light blue", "dark blue")) +
            scale_linetype_manual(values = c("solid", "solid","dashed","solid", "solid"))) +
  theme(legend.position = "bottom")+guides(color=guide_legend(nrow=2,byrow=TRUE))
bf.dich.cat<-bf.dich.long[bf.dich.long$Type.of.Effect=="Cat.",]
plotTheme(ggplot(bf.dich.cat, aes(Method, Proportion, group = Evidence)) +
            geom_line(aes(color = Evidence, linetype = Evidence)) + 
            geom_point(aes(color = Evidence), size = 7, colour = "white") + 
            geom_point(aes(color = Evidence), size = 3, shape = 1) + 
            facet_wrap(~N+Effect.Cat.,nrow=2) + labs(title="\n",y="Proportion\n", x="\nEffect Size") +
            scale_color_manual(values = c("red", "orange","grey","light blue", "dark blue")) +
            scale_linetype_manual(values = c("solid", "solid","dashed","solid", "solid"))) +
  theme(legend.position = "bottom")+guides(color=guide_legend(nrow=2,byrow=TRUE))

