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
    return(drop0(x))
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
#'
#' @rdname RankRulesHypergraph
#' @export RankRulesHypergraph
#'
LearnFunctionalRegulons <- function(rules) {
    items <- items(rules)
    items.mat <- as(items, 'ngCMatrix')
    items.mat <- as(items.mat, 'dgCMatrix')

    weighted.adjacency.matrix <- items.mat %*% t(items.mat)

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
