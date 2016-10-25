rm(list = ls())
library(doParallel)
library(FLCore)
library(rbenchmark)
# 
# 1) figure out a way to import and export FLStock objects
# 2) figure out the best way to have a function and input for optimization

source("~/fmle_bb.R")

data(ple4)
ple4SR<-as.FLSR(ple4)
#### Specifying the stock recruitment relationship and error model
model(ple4SR)<-bevholt()

setSize <- c(2, 3, 4, 5)
resultsDF <- data.frame()
# size = 5
# Loop through the set sizes and perform sub-experiment
for(s in setSize) {

  size <- 10 ^ s
  
  ple4SR_size <- propagate(ple4SR, iter = size)

  # Run performance test
  results <- benchmark(baseFLR     <- FLCore::fmle(ple4SR_size),
                       fmle_bb_par <- fmle_bb(ple4SR_size, inParallel = TRUE),
                       fmle_bb_seq <- fmle_bb(ple4SR_size, inParallel = FALSE),
                       replications = 2,
                       columns = c("test", "elapsed", "replications"))

  # Save results
  # agSum <- summary(results)
  agSumPlus <- data.frame(results, count = size) 
  
  resultsDF <- rbind(resultsDF, agSumPlus)
}
# warnings()

library(ggplot2)
g2 <- ggplot(data=resultsDF, aes(x=(count)))
g2 <- g2 + geom_line(aes(y = mean / 1000, group=expr, colour=expr), size=1) 
g2 <- g2 + scale_fill_hue(l=40)
g2 <- g2 + theme(axis.text.x = element_text(angle=30, vjust=1))
g2 <- g2 + guides(colour = guide_legend("Legend"))
g2 <- g2 + labs(title="R Performance - fmle optimization", x="Iterations", y="Mean Time (seconds)")
g2


