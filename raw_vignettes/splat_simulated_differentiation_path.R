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
library(hashmap)
library(FNN)
library(crules)
library(arules)
library(markovchain)

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
        sparse.count.data, arm.algorithm='apriori',
        supp=0,conf=25, zmax=2)

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

markov.trace <- GenerateMarkovTrace(
        ids, sparse.count.data, rules,
        tau= 0.2, n.sim.steps = 25)

markov.chain.fit <- markovchainFit(markov.trace)

print("All done!")
