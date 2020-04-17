#' A data set of counts of emails between email addresses
#' This is a non-symmetric network.
#' Nodes represent email address.
#' The edge-state ij between two email addresses i and j is the number of emails sent from i to j
#' The Groups vector is the node label from the igraph attribute "notes"
#'
#' @format A list containing
#' \describe{
#'   \item{Edges}{an edges object with each edge-state representing the number of emails between two email addresses}
#'   \item{Groups}{A vector giving a group name to which the email address belong. The order matches the edges such that Edges[i,j] is the edge-state between the nodes i and nodes j who are members of Groups[i] and Groups[j] respectively}
#' }
#' @source \url{https://cran.r-project.org/package=igraphdata}
"Enron"
