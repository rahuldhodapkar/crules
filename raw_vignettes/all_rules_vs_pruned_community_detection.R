# monocytes_poc.R
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

library(Seurat)
library(arules)
library(crules)

load('./data/monocytes.RData')

########################################

# from reads
communities.raw <- readRDS('./data/exp_raw_monocyte_communities.rds')

communities.raw.penalties <- CalcCommPenaltyFromGMT(
    TrimSmallComms(communities.raw, min.community.size=10),
    gmt.filename='./data/c2.all.v7.4.symbols.gmt'
)

# from association rules

x <- Seurat::GetAssayData(object = monocytes, slot='counts')
x <- t(x)

y <- GenerateIncidenceMatrix(x, method='percentile', percentile.cutoff=0.99)

# each ROW should be a CELL, corresponding to a TRANSACTION or ITEMSET.

y1 <- as(y, 'ngCMatrix')
txs <- as(t(y1), 'transactions')

rules <- apriori(txs, 
    parameter=list(
        supp=0.1, 
        conf=0.8,
        maxlen=5),
    control=list(
        memopt=TRUE,
        load=FALSE
    ))
hypergraph.rank <- RankRulesHypergraph(rules)

# add hypergraph rank as interestingness measure
quality(rules)$hgrank <- hypergraph.rank

items <- items(rules)
items.mat <- as(items, 'ngCMatrix')
items.mat <- as(items.mat, 'dgCMatrix')

weighted.adjacency.matrix <- items.mat %*% t(items.mat)
communities.all <- DetectOverlappingCommunitiesSLPAw(
    weighted.adjacency.matrix, 
    num.iters=20, 
    rand.seed=42,
    community.threshold=0)
communities.all.penalties <- CalcCommPenaltyFromGMT(
    TrimSmallComms(communities.all, min.community.size=2),
    gmt.filename='./data/c2.all.v7.4.symbols.gmt'
)

top.hg.rules <- head(
    sort(rules, by='hgrank', decreasing=FALSE),
    n=floor(0.80 * length(rules))
)

items <- items(top.hg.rules)
items.mat <- as(items, 'ngCMatrix')
items.mat <- as(items.mat, 'dgCMatrix')

weighted.adjacency.matrix <- items.mat %*% t(items.mat)
communities.top.hg <- DetectOverlappingCommunitiesSLPAw(
    weighted.adjacency.matrix, 
    num.iters=20, 
    rand.seed=42,
    community.threshold=0)
communities.top.hg.penalties <- CalcCommPenaltyFromGMT(
    TrimSmallComms(communities.top.hg, min.community.size=2),
    gmt.filename='./data/c2.all.v7.4.symbols.gmt'
)

nrow(communities.raw.penalties)
nrow(communities.all.penalties)
nrow(communities.top.hg.penalties)

mean(communities.raw.penalties$best.penalties)
mean(communities.all.penalties$best.penalties)
mean(communities.top.hg.penalties$best.penalties)

median(communities.raw.penalties$best.penalties)
median(communities.all.penalties$best.penalties)
median(communities.top.hg.penalties$best.penalties)

sd(communities.raw.penalties$best.penalties)
sd(communities.all.penalties$best.penalties)
sd(communities.top.hg.penalties$best.penalties)

t.test(communities.all.penalties$best.penalties, communities.top.hg.penalties$best.penalties)

wilcox.test(
        communities.all.penalties$best.penalties[
            communities.all.penalties$best.penalties < (
                    mean(communities.all.penalties$best.penalties)
                    * 3 * sd(communities.all.penalties$best.penalties)
                )
        ],
        communities.raw.penalties$best.penalties[
            communities.raw.penalties$best.penalties < (
                    mean(communities.raw.penalties$best.penalties)
                    * 3 * sd(communities.raw.penalties$best.penalties)
                )
        ],
        paired=FALSE
    )


t.test(
        sqrt(communities.all.penalties$best.penalties),
        sqrt(communities.top.hg.penalties$best.penalties),
        paired=FALSE
    )

t.test(
        communities.all.penalties$best.penalties,
        communities.top.hg.penalties$best.penalties,
        paired=FALSE
    )

hist(log(communities.all.penalties$best.penalties))
hist(log(communities.top.hg.penalties$best.penalties))