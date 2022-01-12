# compare_to_velocyto.R
#     Copyright (C) 2021  Rahul Dhodapkar
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

library(velocyto.R)
library(pagoda2)

library(Seurat)
library(arules)
library(crules)
library(ggplot2)
library(umap)
library(mstknnclust)
library(igraph)
library(cowplot)

ldat <-read.loom.matrices('./data/10X43_1.loom')
emat <- ldat$spliced
emat <- emat[,colSums(emat)>=1e3]

r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
r$adjustVariance(plot=T,do.par=T,gam.k=10)
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine');


r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)

par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Xist"],main='Xist')

emat <- ldat$spliced; nmat <- ldat$unspliced
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]; # restrict to cells that passed p2 filter
# take cluster labels
cluster.label <- r$clusters$PCA[[1]]
cell.colors <- sccore::fac2col(cluster.label)
# take embedding
emb <- r$embeddings$PCA$tSNE

cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)

show.velocity.on.embedding.cor(
    emb,
    rvel.cd,
    n=300,
    scale='sqrt',
    cell.colors=ac(cell.colors,alpha=0.5),
    cex=0.8,
    arrow.scale=5,
    show.grid.flow=TRUE,
    min.grid.cell.mass=0.5,
    grid.n=40,
    arrow.lwd=1,
    do.par=F,
    cell.border.alpha = 0.1)

# embed on UMAP coordinate space.
prcomp.out <- prcomp(emat)
orig.umap <- umap(as.matrix(prcomp.out$rotation[,1:50]))

show.velocity.on.embedding.cor(
    orig.umap$layout,
    rvel.cd,
    n=300,
    scale='sqrt',
    cell.colors=ac(cell.colors,alpha=0.5),
    cex=0.8,
    arrow.scale=5,
    show.grid.flow=TRUE,
    min.grid.cell.mass=0.5,
    grid.n=40,
    arrow.lwd=1,
    do.par=F,
    cell.border.alpha = 0.1)


# now perform our own forecasting using arules.
x <- t(emat)
y <- GenerateIncidenceMatrix(x, method='percentile', percentile.cutoff=0.99)
y1 <- as(y, 'ngCMatrix')
txs <- as(t(y1), 'transactions')

rules <- apriori(txs, 
    parameter=list(
        supp=0.01, 
        conf=0.1,
        maxlen=6),
    control=list(
        memopt=TRUE,
        load=FALSE
    ))

IDENTITY.DECAY.FACTOR <- 1

L <- as(rules@lhs@data, 'dgCMatrix')
R <- as(rules@rhs@data, 'dgCMatrix')

# weight RHS data by confidence per rule.
conf <- rules@quality$confidence

ConfRepFromPVector <- function(p, conf) {
    reps <- c()
    for (i in 2:length(p)) {
        reps <- c(reps, 
            rep(conf[[i-1]], p[[i]] - p[[i-1]])
        )
    }
    return(reps)
}

R@x <- ConfRepFromPVector(R@p, conf)

identity.rules <- as(Diagonal(nrow(L)), 'dgCMatrix')
L <- cbind(L, identity.rules)
R <- cbind(R, identity.rules * IDENTITY.DECAY.FACTOR)


W <- L %*% t(R)

x2 <- x %*% W
x3 <- x2 %*% W
x4 <- x3 %*% W
x5 <- x4 %*% W

x2.scaled <- Seurat::FastRowScale(t(as.matrix(x2)))
x2.scaled <- t(x2.scaled)
colnames(x2.scaled) <- rownames(emat)
x2.pca.coords <- x2.scaled %*% prcomp.out$x[,1:50]
pred.umap <- predict(orig.umap, x2.pca.coords)

# now plot
df <-data.frame(
    cluster=cluster.label,
    x1 = orig.umap$layout[,1],
    y1 = orig.umap$layout[,2],
    x2 = pred.umap[,1],
    y2 = pred.umap[,2]
)
df.sub <- subset(df, cluster.label==10)

ggplot(data=df, aes(
        x=x1, y=y1,
        xend=x2, yend=y2,
        color=cluster)) +
    geom_point(alpha = 0.5) +
    geom_segment(data=df.sub,
        arrow = arrow(length = unit(0.1,"cm")),
        size = 0.6,
        alpha = 0.4)

# use local community label propagation to build a
# transition table for clusters of cells.
dist.mat <- as.matrix(dist(rbind(prcomp.out$rotation[,1:50], x2.pca.coords)))
dist.mat.triplet <- as(dist.mat, 'TsparseMatrix')

knn.G <- generate.knn(
    data.frame(
        object_i = dist.mat.triplet@i,
        object_j = dist.mat.triplet@j,
        d_ij = dist.mat.triplet@x 
    ),
    suggested.k = 5
)
