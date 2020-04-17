#' The Enron data set as extracted from igraph using the script in data-raw
#' @format An \code{edges} object of counts of emails between nodes
#' @seealso igraph
"Enron"

#' The Macaque data set as extracted from igraph using the script in data-raw
#' @format An \code{edges} object of activation counts between brain regions in a Macaque
#' @seealso igraph
"Macaque"

#' The StackOverflow data set as extracted from igraph using the script in data-raw
#' Extracted on 27/8/2019 from kaggle (login required) using:
#' library(rvest)
#' read_html("https://www.kaggle.com/stackoverflow/stack-overflow-tag-network/downloads/stack_network_links.csv/1")
#' @format An \code{edges} object of activation counts between brain regions in a Macaque
#' @source \url{https://www.kaggle.com/stackoverflow/stack-overflow-tag-network/}
#' @seealso igraph
"StackOverflow"
