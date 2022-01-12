# splat_simulated_differentiation_path.R
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

library(splatter) # BiocManager::install("splatter")
library(scater)

set.seed(1)

sim.paths <- splatSimulate(
    de.prob = 0.2,
    nGenes = 1000,
    batchCells = 2000,
    method = 'paths',
    verbose = FALSE)
sim.paths <- logNormCounts(sim.paths)
sim.paths <- runPCA(sim.paths)
plotPCA(sim.paths, colour_by = "Step")

sparse.count.data <- as(t(sim.paths@assays@data$counts), 'dgCMatrix')
sparse.normed.data <- RowScaleSparseMatrix(sparse.count.data)

rules <- GenerateCellularRules(
        sparse.count.data,
        min.support = 0, min.conf = 0.25, max.rule.len = 2)

trans.mat.list <- list()
for (i in 1:25) {
    print(i)
    xprime <- ForecastStates(
        sparse.count.data, rules,
        tau=0.2, n.sim.steps = i
    )

    umap.df <- ProjectSimUMAP(sparse.normed.data,xprime)
    umap.df$Step <- sim.paths$Step

    ggplot(umap.df, aes(
            x=x1, y=y1,
            xend=x2, yend=y2,
            color=Step)) +
        geom_point(alpha = 0.5) +
        geom_segment(
            arrow = arrow(length = unit(0.1,"cm")),
            size = 0.6,
            alpha = 0.4)

    # generate markov chain
    step2id <- hashmap(c(1:25,
              26:50,
              51:75,
              76:100),
            c(rep('g1', 25),
              rep('g2', 25),
              rep('g3', 25),
              rep('g4', 25)))
    ids <- step2id[[sim.paths$Step]]

    trans.mat <- GenerateMarkovChain(ids,sparse.normed.data,xprime)

    trans.mat.list[[i]] <- trans.mat
}

print("All done!")

for(i in 1:25) {
    print(i)
    print(trans.mat.list[[i]])
}

