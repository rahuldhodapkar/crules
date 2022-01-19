# crules.R
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


#' Create incidence matrix from dataset, wrapper function to allow for simpler
#' user interaction to various means of generating incidence matrix.
#'
#' @param x a dgCMatrix of numbers (e.g. RNA counts)
#' @param method a key to indicate which generation method is to be used.
#' @param ... all additional parameters will be passed to downstream functions
#'
#' @return a dgCMatrix incidence matrix
#'
#' @rdname GenerateIncidenceMatrix
#' @export GenerateIncidenceMatrix
#'
GenerateIncidenceMatrix <- function(x, 
                                    method='percentile',
                                    ...) {

    if (method == 'percentile') {
        return(GenerateIncidenceMatrixPercentile(x, ...))
    } else if (method == 'nonzero') {
        return(GenerateIncidenceMatrixNonzero(x, ...))
    } else {
        warning("WARN: invalid option passed to GenerateIncidenceMatrix, returning NULL.")
        return(NULL)
    }
}

#' Sets all nonzero entries to 1, creating incidence matrix.
#'
#' @param x a dgCMatrix of numbers (e.g. RNA counts)
#'
#' @return a dgCMatrix incidence matrix
#'
#' @rdname GenerateIncidenceMatrixNonzero
#' @export GenerateIncidenceMatrixNonzero
#'
GenerateIncidenceMatrixNonzero <- function(x) {
    x@x <- rep(1, length(x@x))
    return(x)
}

#' Create incidence matrix from dataset, using column-wise percentile values to 
#' generate inclusion cutoffs.
#'
#' @param x a dgCMatrix of numbers (e.g. RNA counts)
#' @param percentile.cutoff the percentile at which 
#'        to consider "prescence" in the incidence matrix.
#'        Default = 0.8.
#'
#' @return a dgCMatrix incidence matrix
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom Matrix drop0
#'
#' @rdname GenerateIncidenceMatrixPercentile
#' @export GenerateIncidenceMatrixPercentile
#'
GenerateIncidenceMatrixPercentile <- function(
    x, 
    percentile.cutoff=0.8) {
    ########################

    ordered.val.num <- floor(percentile.cutoff * ncol(x))
    updates <- unlist(pbmclapply(1:ncol(x), function(i) {
        if((x@p[i+1]) == x@p[i]) {
            return() #because column has no positive values
        }

        num.zero.entries <- ncol(x) - (x@p[i+1] - x@p[i])
        if (num.zero.entries >= ordered.val.num) {
            delimiter <- 1
        } else {
            col.i.vals <- sort(x@x[(x@p[i] + 1):(x@p[i+1])],
                                    decreasing=FALSE)
            delimiter <- col.i.vals[[ordered.val.num - num.zero.entries]]
        }

        col.vals.to.update <- c()

        for (val.i in 1:(x@p[i+1]-x@p[i])) {
            if ( x@x[x@p[i] + val.i] < delimiter ) {
                col.vals.to.update <- append(col.vals.to.update, 0)
            } else {
                col.vals.to.update <- append(col.vals.to.update, 1)
            }
        }

        return(col.vals.to.update)
    }))

    x@x <- updates
    return(as(drop0(x), 'dgCMatrix'))
}

#' Generate representation of rule hypergraph from which to derive
#' interestingness measure for rules in the style of CONFIGV
#' (Santolucito et. al, OOPSLA 2017)
#'
#' @param rules a rules output object from arules apriori or eclat
#'
#' @return a vector of length equal to the number of rules specifying 
#'         a rank score by hypergraph degree as defined in
#'         (Santolucito et. al, OOPSLA 2017)
#'
#' @importFrom arules items
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom methods as
#'
#' @rdname RankRulesHypergraph
#' @export RankRulesHypergraph
#'
RankRulesHypergraph <- function(rules) {
    items <- items(rules)
    items.mat <- as(items, 'ngCMatrix')
    items.mat <- as(items.mat, 'dgCMatrix')

    deg <- rowSums(items.mat)
    num.elems.per.rule <- colSums(items.mat)

    deg.weighted.rules <- as.vector(deg %*% items.mat)

    return(deg.weighted.rules / num.elems.per.rule)
}

#' Learn functional regulons from selected association rules on gene set
#'
#' @param rules a rules output object from arules apriori or eclat
#'
#' @return a communities object containing 
#'
#' @importFrom arules items
#' @importFrom methods as
#' @importFrom Matrix t
#'
#' @rdname RankRulesHypergraph
#' @export RankRulesHypergraph
#'
LearnFunctionalRegulons <- function(rules) {
    items <- items(rules)
    items.mat <- as(items, 'ngCMatrix')
    items.mat <- as(items.mat, 'dgCMatrix')

    weighted.adjacency.matrix <- items.mat %*% Matrix::t(items.mat)

    # ***TODO***
    return ()
}

#' Pure R implementation of SLPAw (Xie et. al 2011)
#'
#' @param W a weighted walk matrix
#' @param num.iters a predefined number of iterations after which to terminate
#'        the SLPA algorithm. Default 20, reported as empirically observed 
#'        general stability threshold in the original paper.
#' @param rand.seed an integer seed for the RNG, default=42
#' @param community.threshold a number in [0,1] to determine the thresholding
#'        value for communities to report. The function will report all
#'        communities with a total 
#'        number of labels >= floor(community.threshold * num.iters)
#'        in memory for each node.
#'
#' @return a list of lists of node IDs containing communities identified
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @rdname DetectOverlappingCommunitiesSLPAw
#' @export DetectOverlappingCommunitiesSLPAw
#'
DetectOverlappingCommunitiesSLPAw <- function(
    W, 
    num.iters=20, 
    rand.seed=42,
    community.threshold=0.25) {

    set.seed(rand.seed)

    n <- dim(W)[[1]]
    order <- 1:n

    message('Executing speaker-listener protocol')

    message('Intialization')
    mem <- list()
    for (i in 1:n) {
        mem[[i]] <- list(i)
    }
    diag(W) <- 0
    W <- drop0(W)

    message('Evolution')
    pb <- txtProgressBar(min = 0,
        label = "Running iterative label propagation",
        max = num.iters,
        style = 3)

    for (t in 1:num.iters) {
        rand.order <- sample(x=order, size=n, replace=FALSE)
        for (i in rand.order) {
            neighbors <- which(W[i,] != 0)
            neighbor.weights <- W[i,neighbors]
            if (length(neighbors) == 0) {
                next
            }

            spoken.values <- unlist(lapply(1:length(neighbors), function(j) {
                neighbor.ix <- neighbors[[j]]
                return(mem[[neighbor.ix]][[
                    sample.int(n=length(mem[[neighbor.ix]]), size=1)
                ]])
            }))

            mem[[i]] <- append(mem[[i]], 
                ListenerProbWeighted(spoken.values,weights=neighbor.weights))
        }
        setTxtProgressBar(pb, t)
    }

    message('Post-processing')
    communities <- list()
    for (i in 1:n) {
        communities[[i]] <- list()
    }

    community.min.size <- floor(num.iters * community.threshold)
    for (i in 1:n) {
        ui <- unique(mem[[i]])
        comms.i <- ui[which(
            tabulate(match(mem[[i]], ui)) > community.min.size
        )]

        if (length(comms.i) == 0) {
            next
        }

        for(c in comms.i) {
            communities[[c]] <- append(communities[[c]], i)
        }
    }

    names(communities) <- rownames(W)

    # map back to names
    communities.symb <- list()
    for (i in 1:length(communities)) {
        communities.symb[[i]] <- rownames(W)[unlist(communities[[i]])]
    }
    names(communities.symb) <- names(communities)

    return (communities.symb)
}


#' Compare communities list of lists to *.gmt symbol table from MSigDB and
#' generate a concordance score based on the Jaccard dissimilarity index.
#'
#' @param communities a list of lists containing symbols for gene communities
#' @param gmt.filename a filename for gmt against which to test communities
#'
#' @return a dataframe containing the penalty weights for each community,
#'         as well as the name of the closest gene set in the provided GMT.
#'
#' @importFrom GSA GSA.read.gmt
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @rdname CalcCommPenaltyFromGMT
#' @export CalcCommPenaltyFromGMT
#'
CalcCommPenaltyFromGMT <- function(
    communities, 
    gmt.filename) {
    ###########################
    gmt.data <- GSA.read.gmt(gmt.filename)

    best.penalties <- c()
    best.matched.genesets <- c()

    message('Calculating Community Penalties')
    pb <- txtProgressBar(min = 0,
        label = "Running iterative label propagation",
        max = length(communities),
        style = 3)

    for (c.ix in 1:length(communities)) {
        c <- communities[[c.ix]]

        min.penalty <- Inf
        min.penalty.gs.name <- NULL
        for (gs.ix in 1:length(gmt.data$genesets)) {
            gs <- gmt.data$genesets[[gs.ix]]
            gs.name <- gmt.data$geneset.names[[gs.ix]]

            penalty <-  1 - (
                length(intersect(c, gs)) / length(union(c, gs))
            )

            if (penalty < min.penalty) {
                min.penalty <- penalty
                min.penalty.gs.name <- gs.name
            }
        }

        best.penalties <- c(best.penalties, min.penalty)
        best.matched.genesets <- c(best.matched.genesets, min.penalty.gs.name)

        setTxtProgressBar(pb, c.ix)
    }

    return (data.frame(
        community.names=names(communities),
        best.penalties=best.penalties,
        best.matched.genesets=best.matched.genesets
    ))
}


#' Internal LISTENER function to be used in SLPAw algorithm.
#'
#' @param x a vector of values
#' @param weights an optional set of weights to accompany the values
#' @param ... dummy additional arguments that may be passed to other listeners
#'
#' @return the mode or weightd mode of provided values, or a random number 
#'         amongst tied modes. 
#'
#' @importFrom stats setNames
#' @importFrom stats aggregate
#'
ListenerWeightedModeOrRandom <- function(x, weights=NULL, ...) {
    if(length(x) == 0) {
        stop("ERROR: called 'ListenerWeightedModeOrRandom' with empty 'x' parameter")
    }
    if(is.null(weights)) {
        weights <- rep(1, length(x))
    }

    agg <- setNames(
        aggregate(weights, by=list(x), FUN=sum),
        c('x', 'sum.weight')
    )

    max.ixs <- which(agg$sum.weight == max(agg$sum.weight))
    return(agg$x[[ max.ixs[sample.int(n=length(max.ixs), size=1)] ]])
}

#' Internal LISTENER function to be used in SLPAw algorithm.
#'
#' @param x a vector of values
#' @param ... dummy additional arguments that may be passed to other listeners
#'
#' @return the mode or weightd mode of provided values, or a random number 
#'         amongst tied modes. 
#'
#' @importFrom stats setNames
#' @importFrom stats aggregate
#'
ListenerModeOrRandom <- function(x, ...) {
    if(length(x) == 0) {
        stop("ERROR: called 'ListenerModeOrRandom' with empty 'x' parameter")
    }

    ux <- unique(x)
    modes <- ux[which.max(tabulate(match(x, ux)))]
    return(modes[[sample.int(n=length(modes), size=1)]])
}

#' Internal LISTENER function to be used in SLPAw algorithm.
#'
#' @param x a vector of values
#' @param weights an optional set of weights to accompany the values
#' @param ... dummy additional arguments that may be passed to other listeners
#'
#' @return the mode or weightd mode of provided values, or a random number 
#'         amongst tied modes. 
#'
#' @importFrom stats setNames
#' @importFrom stats aggregate
#'
ListenerProbWeighted <- function(x, weights, ...) {
    if(length(x) == 0) {
        stop("ERROR: called 'ListenerProbWeighted' with empty 'x' parameter")
    }
    if(is.null(weights)) {
        weights <- rep(1, length(x))
    }

    return(x[[sample.int(n=length(x), size=1, prob=weights)]])
}

#' Small helper function to trim small communities from list of lists
#'
#' @param communities a vector of values
#' @param min.community.size an optional set of weights to accompany the values
#'
#' @return a pruned list of lists with only those of size greater than
#'         `min.community.size`.
#'
#' @rdname TrimSmallComms
#' @export TrimSmallComms
#'
TrimSmallComms <- function(communities, min.community.size=2) {
    ixs.to.remove <- c()
    for (i in 1:length(communities)) {
        if (length(communities[[i]]) < min.community.size) {
            ixs.to.remove <- c(ixs.to.remove, i)
        }
    }
    if(length(ixs.to.remove) > 0) {
        communities <- communities[- ixs.to.remove]
    }

    return(communities)
}

#' ngCmatric-specific manipulation, internal helper function
#' to weight values of matrix by a "confidence" level obtained
#' from the association rule learning stage. Confidence may be
#' replaced with any other metric on the rule space
#'
#' @param p the p vector from an ngCMatrix
#' @param vals the vector of values to replace per column:
#'              e.g. confidence
#'
#' @return a vector which may be assigned to the x slot of 
#'          the ngCMatrix used for parameter "p".
#'
NewValsFromPVector <- function(p, vals) {
    reps <- c()
    for (i in 2:length(p)) {
        reps <- c(reps, 
            rep(vals[[i-1]], p[[i]] - p[[i-1]])
        )
    }
    return(reps)
}

#' Fast implementation to column scale a dgCMatrix to sum to 1
#'
#' @param M a dgCMatrix to be column scaled.
#'
#' @return column-scaled matrix M.
#'
#' @rdname ColumnScaleSparseMatrix
#' @export ColumnScaleSparseMatrix
#'
ColumnScaleSparseMatrix <- function(M) {
    if(!'dgCMatrix' %in% class(M)) {
        message('WARN: ColumnScaleSparseMatrix - non-dgCMatrix passed for scaling')
    }
    M@x <- M@x / rep.int(colSums(M), diff(M@p))
    return(M)
}

#' Fast implementation to row scale a dgCMatrix to sum to 1
#'
#' @param M a dgCMatrix to be row scaled.
#'
#' @importFrom Matrix t
#'
#' @return column-scaled matrix M.
#'
#' @rdname RowScaleSparseMatrix
#' @export RowScaleSparseMatrix
#'
RowScaleSparseMatrix <- function(M) {
    if(!'dgCMatrix' %in% class(M)) {
        message('WARN: ColumnScaleSparseMatrix - non-dgCMatrix passed for scaling')
    }
    return(Matrix::t(ColumnScaleSparseMatrix(Matrix::t(M))))
}

#' Simulate the cell states using a walk matrix generated
#' from the association rules, re-normalizing the matrix at
#' each step.
#'
#' @param x a cell (row) by gene (column) ngCMatrix of expression
#'           data to perform a walk along
#' @param W the row-normalized walk matrix
#' @param nsteps the total number of time steps to run the simulation
#'
#' @return a new x' of the original matrix after nsteps iterations.
#'
WalkCellState <- function(x, W, nsteps=5) {
    xprime <- x
    for (i in 1:nsteps) {
        xprime <- xprime %*% W
    }
    return(xprime)
}

#' Generate association rules from original data E^0
#'
#' @param exp.mat a cell (row) by gene (column) expression matrix
#'                 of class 'ngCMatrix'
#' @param arm.algorithm algorithm used for frequent itemset mining,
#'                       supports 'apriori' (default) and 'eclat'.
#' @param supp Default 0, minimum support to pass to arm.algorithm in [0,100]
#' @param conf Default 40, minimum confidence to pass to arm.algorithm in [0,100]
#' @param ... additional parameters to pass for rule mining
#'
#' @importFrom arules fim4r
#' @importFrom Matrix t
#'
#' @return a rules object generated from the single cell data
#'
#' @rdname GenerateCellularRules
#' @export GenerateCellularRules
#'
GenerateCellularRules <- function(
    exp.mat,
    arm.algorithm = 'apriori',
    supp = 0,
    conf = 40,
    ...) {
    # input validation
    if (!arm.algorithm %in% c('apriori', 'eclat', 'fpgrowth')) {
        stop(paste0(
            'ERROR: GenerateCellularRules - unsupported `arm.algorithm`="',
            arm.algorithm, '" passed.'))
    }

    args <- list(...)

    # non-tunable params
    GEN.INCIDENCE.MATRIX.PERCENTILE.CUTOFF <- 0.99

    # function body
    y <- GenerateIncidenceMatrix(
        exp.mat, method='percentile', 
        percentile.cutoff=GEN.INCIDENCE.MATRIX.PERCENTILE.CUTOFF)
    y1 <- as(Matrix::t(y), 'ngCMatrix')
    txs <- as(y1, 'transactions')

    rules <- fim4r(
        transactions = txs,
        method = arm.algorithm,
        target = 'rules',
        supp = supp,
        conf = conf,
        ...
    )

    return(rules)
}

#' Creates a walk matrix from a set of rules.
#'
#' @param rules a rules object
#' @param tau defualt 0.2. Time constant for the walk matrix, in [0,1].
#'
#' @importFrom Matrix t
#'
#' @return a gene x gene walk matrix W.
GenerateWalkMatrix <- function(
        rules,
        tau = 0.2
    ) {
    L <- as(rules@lhs@data, 'dgCMatrix')
    R <- as(rules@rhs@data, 'dgCMatrix')
    conf <- rules@quality$confidence
    R@x <- NewValsFromPVector(R@p, conf)

    T <- RowScaleSparseMatrix(L %*% Matrix::t(R))
    I <- as(Diagonal(nrow(L)), 'dgCMatrix')
    W <- RowScaleSparseMatrix(tau*T + (1-tau)*I)
    return(W)
}

#' Predict cell state transition graph using data simulated using
#' association rules, weighted by rule confidence.
#'
#' @param exp.mat a cell (row) by gene (column) expression matrix
#'                 of class 'ngCMatrix'
#' @param rules a rules object
#' @param tau time constant for simulation
#' @param n.sim.steps total number of simulation steps to run.
#'
#' @return a new expression matrix E^t after `n.sim.steps` simulation
#'          steps with time constant tau
#'
#' @importFrom Matrix Diagonal
#' @importFrom Matrix t
#'
#' @rdname ForecastStates
#' @export ForecastStates
#'
ForecastStates <- function(
        exp.mat,
        rules,
        tau = 0.2,
        n.sim.steps = 5
    ) {

    W <- GenerateWalkMatrix(rules, tau)
    x <- RowScaleSparseMatrix(exp.mat)
    xprime <- WalkCellState(x, W, nsteps=n.sim.steps)

    return(xprime)
}

#' Project a normalized simulated vector of cell expression
#' onto UMAP space given an normalized expression matrix.
#' If these vectors are large (large number of genes), PCA
#' may be performed on both `exp.mat` and `x.prime` to reduce
#' computational load.
#'
#' @param exp.mat original expression data (E^0)
#' @param x.prime predicted expression data (E^t)
#'
#' @importFrom umap umap
#' @importFrom stats predict
#'
#' @return a data frame with umap coordinates for the cells
#'          and simulated states (x1,y1) -> (x2,y2) where each
#'          row corresponds to a cell in the original dataset.
#'
#' @rdname ProjectSimUMAP
#' @export ProjectSimUMAP
#'
ProjectSimUMAP <- function(
        exp.mat,
        x.prime
    ) {
    orig.umap <- umap(as.matrix(exp.mat))
    pred.umap <- predict(orig.umap, as.matrix(x.prime))

    df <- data.frame(
        x1 = orig.umap$layout[,1],
        y1 = orig.umap$layout[,2],
        x2 = pred.umap[,1],
        y2 = pred.umap[,2]
    )

    return(df)
}

#' Calculate new id vector (sigma)
#'
#' @param ids a vector of cell cluster ids to generate 
#'             transition probabilities
#' @param exp.mat original expression data (E^0)
#' @param x.prime predicted expression data (E^t)
#'
#' @importFrom FNN get.knnx
#' @importFrom hashmap hashmap
#' @importFrom nnet which.is.max
#'
#' @return vector of new ids
#'
#' @rdname InferSimulatedCellStates
#' @export InferSimulatedCellStates
#'
InferSimulatedCellStates <- function(
        ids,
        exp.mat,
        x.prime
    ) {
    knn.G <- get.knnx(
        exp.mat,
        x.prime,
        k=10,
        algorithm='kd_tree')
    ix2ident <- hashmap(
        1:length(ids), 
        ids)
    
    vals <- matrix(
        ix2ident[[knn.G$nn.index[,2:10]]],
        nrow=nrow(knn.G$nn.index)
    )
    rownames(vals) <- ids

    getmode <- function(v) {
        uniqv <- unique(v)
        uniqv[nnet::which.is.max(tabulate(match(v, uniqv)))]
    }

    sigma <- sapply(1:nrow(vals), function(i){getmode(vals[i,])})
    return(sigma)
}

#' Generate Markov trace for simulation of `exp.mat` on the provided set of `rules`
#' with time constant tau and for n.sim.steps steps.
#'
#' @param ids a vector of cell identities (e.g. cluster ids, cell types)
#' @param exp.mat a cell (row) by gene (column) expression matrix
#'                 of class 'ngCMatrix'
#' @param rules a rules object
#' @param tau time constant for simulation
#' @param n.sim.steps total number of simulation steps to run.
#'
#' @return a (cell x sim step) matrix capturing cell states at each step of the
#'          simulation for fitting to a markov model.
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @rdname GenerateMarkovTrace
#' @export GenerateMarkovTrace
#'
GenerateMarkovTrace <- function(
        ids,
        exp.mat,
        rules,
        tau = 0.2,
        n.sim.steps = 5
    ) {

    x0 <- RowScaleSparseMatrix(exp.mat)
    x.prime <- x0
    W <- GenerateWalkMatrix(rules, tau)

    pb <- txtProgressBar(min = 0,
    label = "Running cell state simulation",
    max = n.sim.steps,
    style = 3)


    markov.trace <- matrix('', ncol=n.sim.steps + 1, nrow=length(ids))
    markov.trace[,1] <- ids
    for (i in 1:n.sim.steps) {
        x.prime <- x.prime %*% W
        markov.trace[,i+1] <- InferSimulatedCellStates(ids, x0, x.prime)

        setTxtProgressBar(pb, i)
    }

    return(markov.trace)
}
