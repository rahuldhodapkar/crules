# cell_state_forecast.R
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
# Data for this vignette obtained from https://pubmed.ncbi.nlm.nih.gov/31629685/
#

library(Seurat)
library(arules)
library(crules)
library(ggplot2)
library(umap)
library(cowplot)
library(FNN)

#
# Load and clean data
#

sa07890.data <- Read10X(data.dir='./data/10x/GSM3754380/SA07890')
sa07890 <- CreateSeuratObject(counts = sa07890.data, project = "sa07890", min.cells = 3, min.features = 200)
sa07890 <- NormalizeData(sa07890, normalization.method = "LogNormalize", scale.factor = 10000)
sa07890 <- FindVariableFeatures(sa07890, selection.method = "vst", nfeatures = 2000)
sa07890 <- ScaleData(sa07890, features = rownames(sa07890))
sa07890 <- RunPCA(sa07890, features = VariableFeatures(object = sa07890))
sa07890 <- RunUMAP(sa07890, dims=1:50)
sa07890 <- FindNeighbors(sa07890, dims = 1:50)
sa07890 <- FindClusters(sa07890, resolution = 1)
DimPlot(sa07890, reduction = "pca")
DimPlot(sa07890, reduction = "umap", label=T)
sa07890 <- RenameIdents(object=sa07890, 
    `0` = 'ETP',
    `1` = 'DN2',
    `2` = 'ETP',
    `3` = 'DN2',
    `4` = 'ETP',
    `5` = 'ETP',
    `6` = 'DN2',
    `7` = 'ETP',
    `8` = 'DN3',
    `9` = 'DN2',
    `10`= 'ETP',
    `11`= 'ETP',
    `12`= 'ETP')



sa09216.data <- Read10X(data.dir='./data/10x/GSM3754380/FT-SA09216')
sa09216 <- CreateSeuratObject(counts = sa09216.data, project = "sa09216", min.cells = 3, min.features = 200)
sa09216 <- NormalizeData(sa09216, normalization.method = "LogNormalize", scale.factor = 10000)
sa09216 <- FindVariableFeatures(sa09216, selection.method = "vst", nfeatures = 2000)
sa09216 <- ScaleData(sa09216, features = rownames(sa09216))
sa09216 <- RunPCA(sa09216, features = VariableFeatures(object = sa09216))
DimPlot(sa09216, reduction = "pca")

# perform forecasting with arules.
FeaturePlot(sa07890, 'Cd4')

x <- Seurat::GetAssayData(object = sa07890, slot='counts')
x <- t(x)

y <- GenerateIncidenceMatrix(x, method='percentile', percentile.cutoff=0.99)

# each ROW should be a CELL, corresponding to a TRANSACTION or ITEMSET.

y1 <- as(y, 'ngCMatrix')
txs <- as(t(y1), 'transactions')

rules <- apriori(txs, 
    parameter=list(
        supp=0.1, 
        conf=0.4,
        maxlen=3),
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

#' Internal utility function to quickly scale and center a sparse matrix
#'
#' @param M a dgCMatrix to be 
#'
#' @return row-scaled and centered matrix M
#'
ColumnScaleSparseMatrix <- function(M) {
    if(!'dgCMatrix' %in% class(M)) {
        message('WARN: RowScaleSparseMatrix - non-dgCMatrix passed for scaling')
    }
    reps <- c()
    for (i in 2:length(M@p)) {
        reps <- c(reps,
            as.vector(
                scale(
                    M@x[length(reps):(length(reps) + M@p[[i]] - M@p[[i-1]])],
                    center=TRUE,
                    scale=TRUE
                )
            )
        )
    }
    M@x <- reps
    return(M)
}

x.pca.coords <- t(as.matrix(sa07890[["RNA"]]@scale.data)[VariableFeatures(sa07890),]) %*% Loadings(sa07890[["pca"]])

x2.scaled <- Seurat::FastRowScale(t(as.matrix(x2)))
x2.scaled <- t(x2.scaled)
colnames(x2.scaled) <- rownames(sa07890[["RNA"]])
x2.pca.coords <- x2.scaled[,VariableFeatures(sa07890)] %*% Loadings(sa07890[["pca"]])

df <- data.frame(
    x1 = x.pca.coords[,1],
    y1 = x.pca.coords[,2],
    x2 = x2.pca.coords[,1],
    y2 = x2.pca.coords[,2]
)

x5.scaled <- Seurat::FastRowScale(t(as.matrix(x5)))
x5.scaled <- t(x5.scaled)
colnames(x5.scaled) <- rownames(sa07890[["RNA"]])
x5.pca.coords <- x5.scaled[,VariableFeatures(sa07890)] %*% Loadings(sa07890[["pca"]])

orig.umap <- umap(x.pca.coords)
pred.umap <- predict(orig.umap, x2.pca.coords)

df <- data.frame(
    x1 = orig.umap$layout[,1],
    y1 = orig.umap$layout[,2],
    x2 = pred.umap[,1],
    y2 = pred.umap[,2],
    Flt3 = as.vector(sa07890[['RNA']]['Flt3',]),
    Cd34 = as.vector(sa07890[['RNA']]['Cd34',]),
    Kit = as.vector(sa07890[['RNA']]['Kit',]),
    Cd44 = as.vector(sa07890[['RNA']]['Cd44',]),
    Cd25 = as.vector(sa07890[['RNA']]['Il2ra',]),
    Bcl11b = as.vector(sa07890[['RNA']]['Bcl11b',]),
    celltype = Idents(sa07890)
)
df.sub <- subset(df, celltype == 'DN2')
df.sub <- subset(df, Flt3 > 0)
df.sub2 <- subset(df, Bcl11b > 0)

df <- data.frame(
    x1 = x.pca.coords[,1],
    y1 = x.pca.coords[,2],
    x2 = x5.pca.coords[,1],
    y2 = x5.pca.coords[,2],
)

ggplot(df, aes(
        x=x1, y=y1,
        xend=x2, yend=y2,
        color=celltype)) +
    geom_point(alpha = 0.5) +
    geom_segment(data=df.sub,
        arrow = arrow(length = unit(0.1,"cm")),
        size = 0.6,
        alpha = 0.4)

ggplot(df, aes(x=x1,y=y1,color=celltype)) +
    geom_point()


p1 <- ggplot(df, aes(x=x1,y=y1,color=Bcl11b)) + 
    geom_point()
p2 <- ggplot(df, aes(x=x2,y=y2,color=Bcl11b)) + 
    geom_point()
plot_grid(p1,p2)


p1 <- ggplot(df, aes(x=x1,y=y1,color=Cd34)) + 
    geom_point()
p2 <- ggplot(df, aes(x=x2,y=y2,color=Cd34)) + 
    geom_point()
plot_grid(p1,p2)

p1 <- ggplot(df, aes(x=x1,y=y1,color=Flt3)) + 
    geom_point()
p2 <- ggplot(df, aes(x=x2,y=y2,color=Flt3)) + 
    geom_point()
plot_grid(p1,p2)

FeaturePlot(sa07890 , 'Flt3')

# now generate transition probabilities
knn.G <- get.knnx(x.pca.coords, x2.pca.coords, algorithm='kd_tree')


ix2ident <- hashmap(1:length(Idents(sa07890)), as.character(Idents(sa07890)))

new.idents <- matrix()
for (i in 1:nrow(knn.G$nn.index)) {

}
vals <- matrix(
        ix2ident[[knn.G$nn.index[,2:10]]],
        nrow=nrow(knn.G$nn.index)
    )

rownames(vals) <- Idents(sa07890)
# now compute transitions

trans.mat <- matrix(
    0, 
    nrow=length(levels(Idents(sa07890))),
    ncol=length(levels(Idents(sa07890))),
    dimnames=rep(list(levels(Idents(sa07890))), 2))
trans.mat

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

for (i in 1:nrow(vals)) {
    t1 <- rownames(vals)[[i]]
    t2 <- getmode(vals[i,])
    trans.mat[[t1,t2]] <- trans.mat[[t1,t2]] + 1
}
