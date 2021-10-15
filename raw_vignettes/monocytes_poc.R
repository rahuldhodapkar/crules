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

x <- Seurat::GetAssayData(object = monocytes, slot='counts')
x <- t(x)

y <- GenerateIncidenceMatrix(x, method='nonzero')

# each ROW should be a CELL, corresponding to a TRANSACTION or ITEMSET.

y1 <- as(y, 'ngCMatrix')
txs <- as(t(y1), 'transactions')

rules <- apriori(txs, 
    parameter=list(
        supp=0.5, 
        conf=0.8,
        maxlen=3),
    control=list(
        memopt=TRUE,
        load=FALSE
    ))
hypergraph.rank <- RankRulesHypergraph(rules)

# add hypergraph rank as interestingness measure
quality(rules)$hgrank <- hypergraph.rank

inspect(head(sort(rules, by='hgrank', decreasing=FALSE), n=100))

inspect(head(sort(rules, by='lift', decreasing=TRUE), n=10))

print("mined association rules from all monocytes")
