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
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
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
#'
#' @rdname GenerateIncidenceMatrixPercentile
#' @export GenerateIncidenceMatrixPercentile
#'
GenerateIncidenceMatrixPercentile <- function(
    x, 
    percentile.cutoff=0.8) {
    ########################

    ordered.val.num <- floor(percentile.cutoff * ncol(x))
    pbmclapply(1:ncol(x), function(i) {
        if((x@p[i+1]) == x@p[i]) {
            return() #because column has no positive values
        }

        num.zero.entries <- ncol(x) - (x@p[i+1] - x@p[i])
        if (num.zero.entries > ordered.val.num) {
            delimiter <- 1
        } else {
            col.i.vals <- sort(x@x[(x@p[i] + 1):(x@p[i+1])],
                                    decreasing=FALSE)
            delimiter <- col.i.vals[[ordered.val.num - num.zero.entries]]
        }

        for (val.i in 1:(x@p[i+1]-x@p[i])) {
            if ( x@x[x@p[i] + val.i] < delimiter ) {
                x[x@i[x@p[i] + val.i]+1,i] <- 0
            } else {
                x[x@i[x@p[i] + val.i]+1,i] <- 1
            }
        }
    })

    return(x)
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
