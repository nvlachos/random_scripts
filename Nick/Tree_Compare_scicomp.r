#!/usr/bin/env Rscript
# Load required packages
.libPaths("/scicomp/groups/OID/NCEZID/DHQP/CEMB/Rich_DIR/R_Packages/R/x86_64-pc-linux-gnu-library/3.3")
.libPaths("/scicomp/groups/OID/NCEZID/DHQP/CEMB/Rich_DIR/R_Packages/R/3.4")
library("docopt", lib.loc="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Rich_DIR/R_Packages/R/3.4")
library("treespace", lib.loc="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Rich_DIR/R_Packages/R/3.4")
library("phangorn", lib.loc="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Rich_DIR/R_Packages/R/3.4")
'Usage:
   Tree_Compare_scicomp.r [-q <query> -s <subject> -n <ncomparisons> -o <output>]

Options:
   -q Query tree [default: 0]
   -s Subject tree [default: 1]
   -n Number of comparisons [default: 10000]
   -o Output file [default: Tree_Compare.txt]

 ]' -> doc

opts <- docopt(doc)

QT <- read.tree(opts$q)
ST <- read.tree(opts$s)

if(!(is.rooted(QT))){
  QT <- midpoint(QT)
}

if(!(is.rooted(ST))){
  ST <- midpoint(ST)
}

score <- treeDist(QT, ST)

Nodes <- QT$Nnode
Leaves <- Nodes + 1
Labels <- QT$tip.label

RT <- rmtree(opts$n, Leaves, tip.label = Labels)

QDist = c(treeDist(QT, RT[[1]]))
value <- 0
for (i in 2:opts$n){
  value[i] <- treeDist(QT, RT[[i]])
  QDist = c(value[i], QDist)
}

QSD <- sd(QDist)
QMean <- mean(QDist)
QDiff <- abs(score - QMean)
QZ <- QDiff / QSD
QP <- pnorm(-abs(QZ))

SDist = c(treeDist(ST, RT[[1]]))
value <- 0
for (i in 2:opts$n){
  value[i] <- treeDist(ST, RT[[i]])
  SDist = c(value[i], SDist)
}

SSD <- sd(SDist)
SMean <- mean(SDist)
SDiff <- abs(score - SMean)
SZ <- SDiff / SSD
SP <- pnorm(-abs(SZ))

Out_Data = c(opts$q, opts$s, score, QMean, QP, SMean, SP)
write(Out_Data, file = opts$o, ncolumns = 7, sep = "\t")



