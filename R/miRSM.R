## Internal function cluster from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
cluster <- function(graph, 
                    method = "MCL", 
                    expansion = 2, 
                    inflation = 2,
                    hcmethod = "average", 
                    directed = FALSE, 
                    outfile = NULL, ...) {

    #method <- match.arg(method)
    if (method == "FN") {
        graph <- simplify(graph)
        fc <- fastgreedy.community(graph, merges = TRUE, modularity = TRUE)
        membership <- membership(fc)
        if (!is.null(V(graph)$name)) {
            names(membership) <- V(graph)$name
        }
        if (!is.null(outfile)) {
            cluster.save(cbind(names(membership), membership), 
                         outfile = outfile)
        } else {
            return(membership)
        }
    } else if (method == "LINKCOMM") {
        edgelist <- get.edgelist(graph)
        if (!is.null(E(graph)$weight)) {
            edgelist <- cbind(edgelist, E(graph)$weight)
        }
        lc <- getLinkCommunities(edgelist, plot = FALSE, directed = directed,
            hcmethod = hcmethod)
        if (!is.null(outfile)) {
            cluster.save(lc$nodeclusters, outfile = outfile)
        } else {
            return(lc$nodeclusters)
        }
    } else if (method == "MCL") {
        adj <- matrix(rep(0, length(V(graph))^2), nrow = length(V(graph)),
            ncol = length(V(graph)))
        for (i in seq_along(V(graph))) {
            neighbors <- neighbors(graph, v = V(graph)$name[i], mode = "all")
            j <- match(neighbors$name, V(graph)$name, nomatch = 0)
            adj[i, j] = 1
        }
        lc <- mcl(adj, addLoops = TRUE, expansion = expansion, 
            inflation = inflation, allow1 = TRUE, 
            max.iter = 100, ESM = FALSE)
        lc$name <- V(graph)$name
        lc$Cluster <- lc$Cluster

        if (!is.null(outfile)) {
            cluster.save(cbind(lc$name, lc$Cluster), outfile = outfile)
        } else {
            result <- lc$Cluster
            names(result) <- V(graph)$name
            return(result)
        }
    } else if (method == "MCODE") {
        compx <- mcode(graph, vwp = 0.9, haircut = TRUE, fluff = TRUE,
            fdt = 0.1)
        index <- which(!is.na(compx$score))
        membership <- rep(0, vcount(graph))
        for (i in seq_along(index)) {
            membership[compx$COMPLEX[[index[i]]]] <- i
        }
        if (!is.null(V(graph)$name))
            names(membership) <- V(graph)$name
        if (!is.null(outfile)) {
            cluster.save(cbind(names(membership), membership), 
                         outfile = outfile)
            invisible(NULL)
        } else {
            return(membership)
        }
    }
}

## Internal function cluster.save from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
#' @importFrom utils write.table
cluster.save <- function(membership, 
                         outfile) {
  
    wd <- dirname(outfile)
    wd <- ifelse(wd == ".", paste(wd, "/", sep = ""), wd)
    filename <- basename(outfile)
    if ((filename == "") || (grepl(":", filename))) {
        filename <- "membership.txt"
    } else if (grepl("\\.", filename)) {
        filename <- sub("\\.(?:.*)", ".txt", filename)
    }
    write.table(membership, file = paste(wd, filename, sep = "/"), 
                row.names = FALSE, col.names = c("node", "cluster"),
                quote = FALSE)
}

## Internal function mcode.vertex.weighting from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.vertex.weighting <- function(graph, 
                                   neighbors) {
  
    stopifnot(is.igraph(graph))
    weight <- lapply(seq_len(vcount(graph)), function(i) {
        subg <- induced.subgraph(graph, neighbors[[i]])
        core <- graph.coreness(subg)
        k <- max(core)
        ### k-coreness
        kcore <- induced.subgraph(subg, which(core == k))
        if (vcount(kcore) > 1) {
            if (any(is.loop(kcore))) {
                k * ecount(kcore)/choose(vcount(kcore) + 1, 2)
            } else {
                k * ecount(kcore)/choose(vcount(kcore), 2)
            }
        } else {
            0
        }
    })

    return(unlist(weight))
}

## Internal function mcode.find.complex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.find.complex <- function(neighbors, 
                               neighbors.indx, 
                               vertex.weight,
                               vwp, 
                               seed.vertex, 
                               seen) {

    res <- .C("complex", as.integer(neighbors), as.integer(neighbors.indx),
        as.single(vertex.weight), as.single(vwp), as.integer(seed.vertex),
        seen = as.integer(seen), COMPLEX = as.integer(rep(0, length(seen))),
        PACKAGE = "miRSM")

    return(list(seen = res$seen, COMPLEX = which(res$COMPLEX != 0)))
}

## Internal function mcode.find.complexex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.find.complexex <- function(graph, 
                                 neighbors, 
                                 vertex.weight, 
                                 vwp) {

    seen <- rep(0, vcount(graph))

    neighbors <- lapply(neighbors, function(item) {
        item[-1]
    })
    neighbors.indx <- cumsum(unlist(lapply(neighbors, length)))

    neighbors.indx <- c(0, neighbors.indx)
    neighbors <- unlist(neighbors) - 1

    COMPLEX <- list()
    n <- 1
    w.order <- order(vertex.weight, decreasing = TRUE)
    for (i in w.order) {
        if (!(seen[i])) {
            res <- mcode.find.complex(neighbors, neighbors.indx, vertex.weight,
                vwp, i - 1, seen)
            if (length(res$COMPLEX) > 1) {
                COMPLEX[[n]] <- res$COMPLEX
                seen <- res$seen
                n <- n + 1
            }
        }
    }
    rm(neighbors)
    return(list(COMPLEX = COMPLEX, seen = seen))
}

## Internal function mcode.fluff.complex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.fluff.complex <- function(graph, 
                                vertex.weight, 
                                fdt = 0.8, 
                                complex.g,
                                seen) {

    seq_complex.g <- seq_along(complex.g)
    for (i in seq_complex.g) {
        node.neighbor <- unlist(neighborhood(graph, 1, complex.g[i]))
        if (length(node.neighbor) > 1) {
            subg <- induced.subgraph(graph, node.neighbor)
            if (graph.density(subg, loops = FALSE) > fdt) {
                complex.g <- c(complex.g, node.neighbor)
            }
        }
    }

    return(unique(complex.g))
}

## Internal function mcode.post.process from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.post.process <- function(graph, 
                               vertex.weight, 
                               haircut, 
                               fluff, 
                               fdt = 0.8,
                               set.complex.g, 
                               seen) {

    indx <- unlist(lapply(set.complex.g, function(complex.g) {
        if (length(complex.g) <= 2)
            0 else 1
    }))
    set.complex.g <- set.complex.g[indx != 0]
    set.complex.g <- lapply(set.complex.g, function(complex.g) {
        coreness <- graph.coreness(induced.subgraph(graph, complex.g))
        if (fluff) {
            complex.g <- mcode.fluff.complex(graph, vertex.weight, fdt,
                complex.g, seen)
            if (haircut) {
                ## coreness needs to be recalculated
                coreness <- graph.coreness(induced.subgraph(graph, complex.g))
                complex.g <- complex.g[coreness > 1]
            }
        } else if (haircut) {
            complex.g <- complex.g[coreness > 1]
        }
        return(complex.g)
    })
    set.complex.g <- set.complex.g[lapply(set.complex.g, length) > 2]
    return(set.complex.g)
}

## Internal function mcode from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode <- function(graph, 
                  vwp = 0.5, 
                  haircut = FALSE, 
                  fluff = FALSE, 
                  fdt = 0.8,
                  loops = TRUE) {

    stopifnot(is.igraph(graph))
    if (vwp > 1 | vwp < 0) {
        stop("vwp must be between 0 and 1")
    }
    if (!loops) {
        graph <- simplify(graph, remove.multiple = FALSE, remove.loops = TRUE)
    }
    neighbors <- neighborhood(graph, 1)
    W <- mcode.vertex.weighting(graph, neighbors)
    res <- mcode.find.complexex(graph, neighbors = neighbors, vertex.weight = W,
        vwp = vwp)
    COMPLEX <- mcode.post.process(graph, vertex.weight = W, haircut = haircut,
        fluff = fluff, fdt = fdt, res$COMPLEX, res$seen)
    score <- unlist(lapply(COMPLEX, function(complex.g) {
        complex.g <- induced.subgraph(graph, complex.g)
        if (any(is.loop(complex.g)))
            score <- ecount(complex.g)/choose(vcount(complex.g) + 1, 2) *
                vcount(complex.g) 
        else score <- ecount(complex.g)/choose(vcount(complex.g), 2) *
            vcount(complex.g)
        return(score)
    }))
    order_score <- order(score, decreasing = TRUE)
    return(list(COMPLEX = COMPLEX[order_score], score = score[order_score]))
}

## Internal function integer.edgelist from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
integer.edgelist <- function(network)
  # Returns an edge list with integer numbers replacing elements of the network.
{
  if(!is.character(network)){
    cn <- cbind(as.character(network[,1]),as.character(network[,2]))
  }else{
    cn <- network
  }
  nodes <- unique(as.character(t(cn)))
  ids <- seq(nodes)
  names(ids) <- nodes
  g <- matrix(ids[t(cn)],nrow(cn),ncol(cn),byrow=TRUE)
  ret <- list()
  ret$edges <- g
  ret$nodes <- ids
  return(ret)
}

## Internal function edge.duplicates from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
edge.duplicates <- function(network, verbose = TRUE)
  # Finds and removes loops, duplicate edges, and bi-directional edges.
{
  xx <- cbind(as.character(network[,1]),as.character(network[,2]))
  edges <- integer.edgelist(network)$edges
  ne <- nrow(edges)
  loops <- rep(0,ne)
  dups <- rep(0,ne)
  
  out <- .C("edgeDuplicates",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(ne), loops = as.integer(loops), dups = as.integer(dups), as.logical(verbose),
            PACKAGE = "miRSM")
  
  if(verbose){cat("\n")}
  
  loops <- which(out$loops == 1)
  dups <- which(out$dups == 1)
  inds <- unique(c(loops,dups))
  ret <- list()
  ret$inds <- inds
  if(length(inds)>0){
    ret$edges <- xx[-inds,]
  }else{
    ret$edges <- xx
  }
  if(verbose){
    if(length(loops)>0){
      cat("   Found and removed ",length(loops)," loop(s)\n",sep="")
    }
    if(length(dups)>0){
      cat("   Found and removed ",length(dups)," duplicate edge(s)\n",sep="")
    }
  }
  return(ret)
}

## Internal function getLinkCommunities from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
getLinkCommunities <- function(network, hcmethod = "average", use.all.edges = FALSE, edglim = 10^4, directed = FALSE, dirweight = 0.5, bipartite = FALSE, dist = NULL, plot = TRUE, check.duplicates = TRUE, removetrivial = TRUE, verbose = TRUE) 
  # network is an edge list. Nodes can be ASCII names or integers, but are always treated as character names in R.
  # If plot is true (default), a dendrogram and partition density score as a function of dendrogram height are plotted side-by-side.
  # When there are more than "edglim" edges, hierarchical clustering is carried out via temporary files written to disk using compiled C++ code.
{
  
  if(is.character(network) && !is.matrix(network)){
    if(file.access(network) == -1){
      stop(cat("\nfile not found: \"",network,"\"\n",sep=""))
    }else{
      network <- read.table(file = network, header = FALSE)
    }
  }
  x <- network
  rm(network)
  
  if(ncol(x)==3){
    wt <- as.numeric(as.character(x[,3]))
    if(length(which(is.na(wt)==TRUE))>0){
      stop("\nedge weights must be numerical values\n")
    }
    x <- cbind(as.character(x[,1]),as.character(x[,2]))
  }else if(ncol(x)==2){
    x <- cbind(as.character(x[,1]),as.character(x[,2]))
    wt <- NULL
  }else{
    stop("\ninput data must be an edge list with 2 or 3 columns\n")
  }
  
  if(check.duplicates){
    dret <- edge.duplicates(x, verbose = verbose)
    x <- dret$edges
    if(!is.null(wt)){
      if(length(dret$inds) > 0){
        wt <- wt[-dret$inds]
      }
    }
    rm(dret)
  }
  
  el <- x # Modified edge list returned to user.
  rm(x)
  len <- nrow(el) # Number of edges.
  nnodes <- length(unique(c(as.character(el[,1]),as.character(el[,2])))) # Number of nodes.
  
  intel <- integer.edgelist(el) # Edges with numerical node IDs.
  edges <- intel$edges
  node.names <- names(intel$nodes)
  numnodes <- length(node.names)
  
  if(bipartite){
    # Check that network is bipartite.
    big <- graph_from_edgelist(as.matrix(el), directed = directed)
    bip.test <- bipartite.mapping(big)
    if(!bip.test$res){
      stop("\nnetwork is not bi-partite; change bipartite argument to FALSE\n")
    }
    bip <- rep(1,length(bip.test$type))
    bip[which(bip.test$type==FALSE)] <- 0
    names(bip) <- V(big)$name
    bip <- bip[match(node.names, names(bip))]
    rm(big, bip.test)
  }else{
    bip <- 0
  }
  
  rm(intel)
  
  # Switch depending on size of network.
  if(len <= edglim){
    disk <- FALSE
    if(is.null(dist)){
      emptyvec <- rep(1,(len*(len-1))/2)
      if(!is.null(wt)){ weighted <- TRUE}else{ wt <- 0; weighted <- FALSE}
      if(!use.all.edges){
        dissvec <- .C("getEdgeSimilarities",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),rowlen=integer(1),weights=as.double(wt),as.logical(directed),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = as.double(emptyvec), as.logical(bipartite), as.logical(verbose),
                      PACKAGE = "miRSM")$dissvec
      }else{
        dissvec <- .C("getEdgeSimilarities_all",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),as.integer(numnodes),rowlen=integer(1),weights=as.double(wt),as.logical(FALSE),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = as.double(emptyvec), as.logical(bipartite), as.logical(verbose),
                      PACKAGE = "miRSM")$dissvec
      }
      distmatrix <- matrix(1,len,len)
      distmatrix[lower.tri(distmatrix)] <- dissvec
      colnames(distmatrix) <- 1:len
      rownames(distmatrix) <- 1:len
      distobj <- as.dist(distmatrix) # Convert into 'dist' object for hclust.
      rm(distmatrix)
    }else{	
      # Did the user provide an adequate distance matrix?
      if(!inherits(dist,"dist")){
        stop("\ndistance matrix must be of class \"dist\" (see ?as.dist)\n")
      }else if(attr(dist,which="Size") != len){
        stop("\ndistance matrix size must equal the number of edges in the input network\n")
      }else if(length(dist) != (len*(len-1))/2){
        stop("\ndistance matrix must be the lower triangular matrix of a square matrix\n")
      }
      distobj <- dist
    }
    if(verbose){
      cat("\n   Hierarchical clustering of edges...")
    }
    #if(hcmethod=="energy"){
    #	hcedges <- energy.hclust(distobj)
    #}else{
    #	hcedges <- hclust(distobj, method = hcmethod)
    #	}
    hcedges <- hclust(distobj, method = hcmethod)
    hcedges$order <- rev(hcedges$order)
    rm(distobj)
    if(verbose){cat("\n")}
  }else{
    disk <- TRUE
    if(!is.null(wt)){ weighted <- TRUE}else{ wt <- 0; weighted <- FALSE}
    if(!use.all.edges){
      rowlen <- .C("getEdgeSimilarities",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),rowlen=integer(len-1),weights=as.double(wt),as.logical(directed),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = double(1), as.logical(bipartite), as.logical(verbose),
                   PACKAGE = "miRSM")$rowlen
    }else{
      rowlen <- .C("getEdgeSimilarities_all",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),as.integer(numnodes),rowlen=integer(len-1),weights=as.double(wt),as.logical(FALSE),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = double(1), as.logical(bipartite), as.logical(verbose),
                   PACKAGE = "miRSM")$rowlen
    }
    if(verbose){cat("\n")}
    hcobj <- .C("hclustLinkComm",as.integer(len),as.integer(rowlen),heights = single(len-1),hca = integer(len-1),hcb = integer(len-1), as.logical(verbose),
                PACKAGE = "miRSM")
    if(verbose){cat("\n")}
    hcedges<-list()
    hcedges$merge <- cbind(hcobj$hca, hcobj$hcb)
    hcedges$height <- hcobj$heights
    
    hcedges$order <- .C("hclustPlotOrder",as.integer(len),as.integer(hcobj$hca),as.integer(hcobj$hcb),order=integer(len),
                        PACKAGE = "miRSM")$order
    hcedges$order <- rev(hcedges$order)
    hcedges$method <- "single"
    class(hcedges) <- "hclust"
    
  }
  
  # Calculate link densities, cut the tree, and extract optimal clusters.
  
  hh <- unique(round(hcedges$height, digits = 5)) # Round to 5 digits to prevent numerical instability affecting community formation.
  countClusters <- function(x,ht){return(length(which(ht==x)))}
  clusnums <- sapply(hh, countClusters, ht = round(hcedges$height, digits = 5)) # Number of clusters at each height.
  numcl <- length(clusnums)
  
  ldlist <- .C("getLinkDensities",as.integer(hcedges$merge[,1]), as.integer(hcedges$merge[,2]), as.integer(edges[,1]), as.integer(edges[,2]), as.integer(len), as.integer(clusnums), as.integer(numcl), pdens = double(length(hh)), heights = as.double(hh), pdmax = double(1), csize = integer(1), as.logical(removetrivial), as.logical(bipartite), as.integer(bip), as.logical(verbose),
               PACKAGE = "miRSM")
  
  pdens <- c(0,ldlist$pdens)
  heights <- c(0,hh)
  pdmax <- ldlist$pdmax
  csize <- ldlist$csize
  
  if(csize == 0){
    stop("\nno clusters were found in this network; maybe try a larger network\n")
  }
  
  if(verbose){
    cat("\n   Maximum partition density = ",max(pdens),"\n")
  }
  
  # Read in optimal clusters from a file.
  clus <- list()
  for(i in 1:csize){
    if(verbose){
      mes<-paste(c("   Finishing up...1/4... ",floor((i/csize)*100),"%"),collapse="")
      cat(mes,"\r")
      flush.console()
    }
    clus[[i]] <- scan(file = "linkcomm_clusters.txt", nlines = 1, skip = i-1, quiet = TRUE)
  }
  
  file.remove("linkcomm_clusters.txt")
  
  # Extract nodes for each edge cluster.
  ecn <- data.frame()
  ee <- data.frame()
  lclus <- length(clus)
  for(i in 1:lclus){
    if(verbose){
      mes<-paste(c("   Finishing up...2/4... ",floor((i/lclus)*100),"%"),collapse="")
      cat(mes,"\r")
      flush.console()
    }
    ee <- rbind(ee,cbind(el[clus[[i]],],i))
    nodes <- node.names[unique(c(edges[clus[[i]],]))]
    both <- cbind(nodes,rep(i,length(nodes)))
    ecn <- rbind(ecn,both)
  }
  colnames(ecn) <- c("node","cluster")
  colnames(ee) <- c("node1","node2","cluster")
  
  # Extract the node-size of each edge cluster and order largest to smallest.
  ss <- NULL
  unn <- unique(ecn[,2])
  lun <- length(unn)
  for(i in 1:length(unn)){
    if(verbose){
      mes<-paste(c("   Finishing up...3/4... ",floor((i/lun)*100),"%"),collapse="")
      cat(mes,"\r")
      flush.console()
    }
    ss[i] <- length(which(ecn[,2]==unn[i]))
  }
  names(ss) <- unn
  ss <- sort(ss,decreasing=T)
  
  # Extract the number of edge clusters that each node belongs to.
  unn <- unique(ecn[,1])
  
  iecn <- as.integer(as.factor(ecn[,1]))
  iunn <- unique(iecn)
  lunn <- length(iunn)
  nrows <- nrow(ecn)
  
  oo <- rep(0,lunn)
  
  oo <- .C("getNumClusters", as.integer(iunn), as.integer(iecn), counts = as.integer(oo), as.integer(lunn), as.integer(nrows), as.logical(verbose),
           PACKAGE = "miRSM")$counts
  
  names(oo) <- unn
  
  if(verbose){cat("\n")}
  
  pdplot <- cbind(heights,pdens)
  
  # Add nodeclusters of size 0.
  missnames <- setdiff(node.names,names(oo))
  m <- rep(0,length(missnames))
  names(m) <- missnames
  oo <- append(oo,m)
  
  all <- list()
  
  all$numbers <- c(len,nnodes,length(clus)) # Number of edges, nodes, and clusters.
  all$hclust <- hcedges # Return the 'hclust' object. To plot the dendrogram: 'plot(lcobj$hclust,hang=-1)'
  all$pdmax <- pdmax # Partition density maximum height.
  all$pdens <- pdplot # Add data for plotting Partition Density as a function of dendrogram height.
  all$nodeclusters <- ecn # n*2 character matrix of node names and the cluster ID they belong to.
  all$clusters <- clus # Clusters of edge IDs arranged as a list of lists.
  all$edges <- ee # Edges and the clusters they belong to, arranged so we can easily put them into an edge attribute file for Cytoscape.
  all$numclusters <- sort(oo,decreasing=TRUE) # The number of clusters that each node belongs to (named vector where the names are node names).
  all$clustsizes <- ss # Cluster sizes sorted largest to smallest (named vector where names are cluster IDs).
  all$igraph <- graph_from_edgelist(el, directed = directed) # igraph graph.
  all$edgelist <- el # Edge list.
  all$directed <- directed # Logical indicating if graph is directed or not.
  all$bipartite <- bipartite # Logical indicating if graph is bipartite or not.
  
  class(all) <- "linkcomm"
  
  if(plot){
    if(verbose){
      cat("   Plotting...\n")
    }
    if(len < 1500){ # Will be slow to plot dendrograms for large networks.
      if(len < 500){
        all <- plot(all, type="summary", droptrivial = removetrivial, verbose = verbose)
      }else{ # Slow to reverse order of large dendrograms.
        all <- plot(all, type="summary", right = FALSE, droptrivial = removetrivial, verbose = verbose)
      }
    }else if(len <= edglim){
      oldpar <- par(no.readonly = TRUE)
      par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,2.1))
      plot(hcedges,hang=-1,labels=FALSE)
      abline(pdmax,0,col='red',lty=2)
      plot(pdens,heights,type='n',xlab='Partition Density',ylab='Height')
      lines(pdens,heights,col='blue',lwd=2)
      abline(pdmax,0,col='red',lty=2)
      par(oldpar)
    }else{
      plot(heights,pdens,type='n',xlab='Height',ylab='Partition Density')
      lines(heights,pdens,col='blue',lwd=2)
      abline(v = pdmax,col='red',lwd=2)
    }
  }
  
  return(all)
  
}

## Internal function plot.linkcomm from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plot.linkcomm <- function(x, type = "", ...)
  # S3 method for "plot" generic function.
  # x is a "linkcomm" object.
{
  switch(type,
         summary = plotLinkCommSumm(x, ...),
         members = plotLinkCommMembers(x, ...),
         dend = plotLinkCommDend(x, ...),
         graph = plotLinkCommGraph(x, ...),
         commsumm = plotLinkCommSummComm(x, ...)
  )
}

## Internal function plotLinkCommSumm from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plotLinkCommSumm <- function(x, col = TRUE, pal = brewer.pal(9,"Set1"), right = TRUE, droptrivial = TRUE, verbose = TRUE, ...)
  # x is a "linkcomm" object.
{
  oldpar <- par(no.readonly = TRUE)
  # Set up and colour clusters in dendrogram.
  if(is.null(x$dendr)){
    dd <- as.dendrogram(x$hclust)
  }else{
    dd <- x$dendr
  }
  if(col && is.null(x$dendr)){
    cl <- unlist(x$clusters)
    crf <- colorRampPalette(pal,bias=1)
    cols <- crf(length(x$clusters))
    cols <- sample(cols,length(x$clusters),replace=FALSE)
    numnodes <- nrow(x$hclust$merge) + length(which(x$hclust$merge[,1]<0)) + length(which(x$hclust$merge[,2]<0))
    dd <- dendrapply(dd,.COL,height=x$pdmax,clusters=cl,cols=cols,labels=FALSE,numnodes=numnodes,droptrivial = droptrivial,verbose=verbose)
    if(verbose){cat("\n")}
    assign("i",0,environment(.COL))
    assign("memb",0,environment(.COL))
    assign("first",0,environment(.COL))
    assign("left",0,environment(.COL))
  }
  if(is.null(x$dendr) && right){ dd <- rev(dd)}
  grid.newpage()
  plot.new()
  # Set margin.
  margin<-unit(0.045,"npc")
  pushViewport(viewport(x=margin,y=margin,width=unit(1,"npc")-2*margin,height=unit(1,"npc")-2*margin,just=c("left","bottom")))
  pushViewport(viewport(layout=grid.layout(nrow=4,ncol=4,widths=unit(c(1,0.01,0.2,0.09),units=c("null",rep("native",3))),heights=unit(c(0.04,0.79,0.025,0.05),units=rep("npc",4)),respect=TRUE)))
  # Plot dendrogram using base plot function.
  pushViewport(viewport(layout.pos.row=1:3,layout.pos.col=1))
  #return(gridPLT())
  gpl <- c(0.009,0.71,0.13,0.902)
  par(oma=rep(0,4),mar=rep(0,4),ann=FALSE,omd=c(0,1,0,1),pty="m",mgp=rep(0,3),fig = gpl,xpd=NA,new=TRUE)
  plot(dd,axes=FALSE,leaflab="none")
  popViewport(1)
  # Title.
  pushViewport(viewport(layout.pos.row=1,layout.pos.col=1:3))
  title <- grid.text("Link Communities Dendrogram",x = unit(0.5,"npc"),y = unit(2,"npc"),draw=FALSE,name="title")
  title <- editGrob(title,gp = gpar(fontsize=14))
  grid.draw(title)
  popViewport(1)
  # Plot link partition densities.
  numzeros <- -1*log10(max(x$pdens[,2]))
  if(numzeros <= 1){ # Prevent partition density axis from being rounded to 0.
    rr <- 1
  }else{
    rr <- trunc(numzeros)+1
  }
  if(round(max(x$pdens[,2]), digits = rr) > max(x$pdens[,2])){
    xscale_add <- round(max(x$pdens[,2]), digits = rr) + 0.05*max(x$pdens[,2]) # Add 5% to part density x-axis.
    xaxs_max <- round(max(x$pdens[,2]), digits = rr)
  }else{
    xscale_add <- max(x$pdens[,2]) + 0.075*max(x$pdens[,2]) # 7.5% of max partition density added to x-axis.
    xaxs_max <- round(max(x$pdens[,2]), digits = (rr+1))
  }
  pushViewport(viewport(layout.pos.row=2, layout.pos.col=3, xscale=c(0,xscale_add),yscale=c(0,1)))
  ph <- x$pdens[,1]/max(x$pdens[,1])
  max <- x$pdmax/max(x$pdens[,1])
  
  grid.lines(x$pdens[,2],ph,gp=gpar(col='blue',lwd=2),default.units="native")
  xticks <- seq(0, xaxs_max, length.out=3)
  xa <- grid.xaxis(at = xticks,draw=FALSE,name="xa")
  xa <- editGrob(xa,gPath="ticks",y1 = unit(-0.02,"npc"))
  xa <- editGrob(xa,gPath="labels",gp = gpar(fontsize=10),y = unit(-0.04,"npc"))
  grid.draw(xa)
  xl <- grid.text("Partition Density",x=unit(0.5, "npc"), y = unit(-0.08, "npc"),draw=FALSE,name="xl")
  xl <- editGrob(xl,gp = gpar(fontsize=10))
  grid.draw(xl)
  popViewport(1)
  pushViewport(viewport(layout.pos.row=2,layout.pos.col=1:3))
  grid.lines(x=c(0,1),y=c(max,max),gp = gpar(col="red",lty=2,lwd=2))
  popViewport(1)
  pushViewport(viewport(layout.pos.row=2,layout.pos.col=4))
  yticks <- c(0,0.2,0.4,0.6,0.8,1)
  ya <- grid.yaxis(at = yticks,name="ya",draw=FALSE)
  ya <- editGrob(ya,gPath="ticks",x1 = unit(0.1,"npc"))
  ya <- editGrob(ya,gPath="labels",gp = gpar(fontsize=10),x = unit(0.2,"npc"))
  ya <- editGrob(ya,gPath="labels",just = c("left","centre"))
  if(max(x$pdens[,1]<1)){
    yu <- seq(0,round(max(x$pdens[,1]),2),length.out=6)
    roundS <- function(x){return(round(x,2))}
    yu <- sapply(yu,roundS)
    ya <- editGrob(ya,gPath="labels",label=as.character(yu))
  }
  grid.draw(ya)
  yl <- grid.text("Height",x=unit(0.8, "npc"), y = unit(0.5, "npc"),rot=90,draw=FALSE,name="yl")
  yl <- editGrob(yl,gp = gpar(fontsize=10))
  grid.draw(yl)
  popViewport(1)
  # Summary statistics.
  pushViewport(viewport(layout.pos.row=4,layout.pos.col=1))
  summ <- paste("# edges = ",x$numbers[1],",   ","# nodes = ",x$numbers[2],"\n# clusters = ",x$numbers[3],",   Largest cluster = ",x$clustsizes[1]," nodes\nHclust method: ",x$hclust$method)
  ne <- grid.text(summ,x=unit(0.5,"npc"),y=unit(0.1,"npc"),draw=FALSE,name="ne")
  ne <- editGrob(ne,gp = gpar(fontsize=11))
  grid.draw(ne)
  popViewport(1)
  # Return linkcomm object with dendrogram so we don't have to generate it again in the future.
  if(is.null(x$dendr)){
    x$dendr <- dd
    return(x)
  }
  popViewport(0)
  par(oldpar)
}

## Internal function plotLinkCommMembers from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plotLinkCommMembers <- function(x, nodes = head(names(x$numclusters),10), pal = brewer.pal(11,"Spectral"), shape = "rect", total=TRUE, fontsize=11, nspace = 3.5, maxclusters = 20)
  # Plots community membership matrix using a community-specific colour scheme.
  # x is a "linkcomm" object.
{
  # Construct community matrix.
  comms <- unique(x$nodeclusters[as.character(x$nodeclusters[,1])%in%nodes,2]) # Community (cluster) IDs.
  if(length(comms) > maxclusters){
    comms <- comms[1:maxclusters]
  }
  commatrix <- getCommunityMatrix(x,nodes=nodes)
  crf <- colorRampPalette(pal,bias=1)
  cols <- crf(length(comms))
  grid.newpage()
  # Set margin.
  if(total){
    C <- 2; R <- 3
    nodesums <- apply(commatrix,1,sum)
    commsums <- apply(commatrix,2,sum)
  }else{
    C <- 1; R <- 2
  }
  margin<-unit(0.1,"lines")
  pushViewport(viewport(x=1,y=1,width=unit(1,"npc")-2*margin,height=unit(1,"npc")-2*margin,just=c("right","top")))
  pushViewport(viewport(layout=grid.layout(nrow=length(nodes)+R,ncol=length(comms)+C,widths=unit(c(nspace,rep(1,length(comms)+C-1)),rep("null",length(comms)+C)),heights=unit(rep(1,length(nodes)+R),rep("null",length(nodes)+R)),respect=TRUE)))
  # Titles.
  pushViewport(viewport(layout.pos.row=1,layout.pos.col=2:length(comms)+1))
  ctitle <- grid.text("Community Membership",x = unit(0.5,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="ctitle")
  ctitle <- editGrob(ctitle,gp = gpar(fontsize=14))
  grid.draw(ctitle)
  popViewport(1)
  # Draw membership coloured squares/circles/polygons.
  for(i in 1:(length(nodes)+R-2)){
    if(i != length(nodes)+1){
      pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=1))
      nname <- grid.text(as.character(nodes[i]),x = unit(0.9,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="nname")
      nname <- editGrob(nname,gp = gpar(fontsize=fontsize),just="right")
      grid.draw(nname)
      popViewport(1)
    }
    for(j in 1:(length(comms)+C-1)){
      if(total && j == length(comms)+1 && i != length(nodes)+1){
        pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=j+1))
        ntot <- grid.text(nodesums[i],x = unit(0.5,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="ntot")
        ntot <- editGrob(ntot,gp = gpar(fontsize=12))
        grid.draw(ntot)
        popViewport(1)
        if(i == 1){
          pushViewport(viewport(layout.pos.row=2,layout.pos.col=j+1))
          rt <- grid.text(expression(Sigma),x = unit(0.5,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="rt")
          rt <- editGrob(rt,gp = gpar(fontsize=12))
          grid.draw(rt)
          popViewport(1)
        }
      }else{
        if(i == 1 && j != length(comms)+1){
          pushViewport(viewport(layout.pos.row=2,layout.pos.col=j+1))
          rtitle <- grid.text(comms[j],x = unit(0.5,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="rtitle")
          rtitle <- editGrob(rtitle,gp = gpar(fontsize=12))
          grid.draw(rtitle)
          popViewport(1)
        }
        if(total && i == length(nodes)+1 && j != length(comms)+1){
          if(j==1){
            pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=1))
            ct <- grid.text(expression(Sigma),x = unit(0.9,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="ct")
            ct <- editGrob(ct,gp = gpar(fontsize=12))
            grid.draw(ct)
            popViewport(1)
          }
          pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=j+1))
          ctot <- grid.text(commsums[j],x = unit(0.5,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="ctot")
          ctot <- editGrob(ctot,gp = gpar(fontsize=12))
          grid.draw(ctot)
          popViewport(1)
        }else if(i != length(nodes)+1 && j != length(comms)+1){
          if(commatrix[i,j] == 1){
            fill <- cols[j]
          }else{
            fill <- "white"
          }
          if(shape=="rect"){
            pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=j+1))
            grid.rect(gp=gpar(fill=fill,col="grey"),width = unit(0.9,"npc"), height = unit(0.9,"npc"),draw=TRUE)
            popViewport(1)
          }else if(shape=="circle"){
            pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=j+1))
            grid.circle(x=0.5,y=0.5,r=0.45,gp=gpar(fill=fill,col="grey"),draw=TRUE)
            popViewport(1)
          }
        }
      }
    }
  }
}

## Internal function plotLinkCommDend from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plotLinkCommDend <- function(x, col=TRUE, pal = brewer.pal(9,"Set1"), height=x$pdmax, right = FALSE, labels=FALSE, plotcut=TRUE, droptrivial = TRUE, leaflab = "none", verbose = TRUE, ...)
  # x is a "linkcomm" object.
{
  dd <- as.dendrogram(x$hclust)
  if(col){
    cl <- unlist(x$clusters)
    crf <- colorRampPalette(pal,bias=1)
    cols <- crf(length(x$clusters))
    cols <- sample(cols,length(x$clusters),replace=FALSE)
    numnodes <- nrow(x$hclust$merge) + length(which(x$hclust$merge[,1]<0)) + length(which(x$hclust$merge[,2]<0))
    dd <- dendrapply(dd, .COL, height=height, clusters=cl, cols=cols, labels=labels, numnodes = numnodes, droptrivial = droptrivial, verbose = verbose)
    cat("\n")
    assign("i",0,environment(.COL))
    assign("memb",0,environment(.COL))
    assign("first",0,environment(.COL))
    assign("left",0,environment(.COL))
  }
  if(right){
    dd <- rev(dd)
  }
  plot(dd,ylab="Height", leaflab = leaflab, ...)
  if(plotcut){
    abline(h=height,col='red',lty=2,lwd=2)
  }
  #ll <- sapply(x$clusters,length)
  #maxnodes <- length(unique(x$nodeclusters[x$nodeclusters[,2]%in%which(ll==max(ll)),1]))
  summ <- paste("# clusters = ",length(x$clusters),"\nLargest cluster = ",x$clustsizes[1]," nodes")
  mtext(summ, line = -28)
}

## Internal function from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
.COL<-local({
  
  memb <- 0
  first <- 0
  i <- 0
  left <- 0
  
  colorHclusters <<- function(x, height, clusters, cols, labels, numnodes, droptrivial, verbose)
    # Adds colours to edges that belong to clusters below "height" in the dendrogram.
    # Clusters gives leaf IDs for clusters that should be coloured.
    # x is a node in the tree.
  {
    left <<- left + 1
    if(verbose){
      out <- paste(c("   Colouring dendrogram... ",floor((left/numnodes)*100),"%"),collapse="")
      cat(out,"\r")
      flush.console()
    }
    
    if(round(attributes(x)$height,digits=5) > height){
      return(x)
    }else{
      if(is.leaf(x)){
        if(is.na(match(as.numeric(attributes(x)$label),clusters))){
          if(!labels){
            attributes(x)$label <- NULL
          }
          return(x)
        }
      }
      a <- attributes(x)
      if(memb == 0){
        memb <<- attributes(x)$members
        if(droptrivial == TRUE && memb == 2){
          memb <<- 0
          first <<- 1
        }else{
          i <<- i+1
          first <<- 1 # Because we don't colour the edge leading to the first node in a cluster.
        }
      }
      if(first == 0){
        attr(x,"edgePar") <- c(a$edgePar,list(col = cols[i], lwd = 2))
      }
      if(is.leaf(x)){
        if(!labels){
          attributes(x)$label <- NULL
        }
        memb <<- memb-1
      }
      first <<- 0
    }
    
    return(x)
  }
  
})

## Internal function plotLinkCommGraph from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plotLinkCommGraph <- function(x, clusterids = 1:length(x$clusters), nodes = NULL, layout = layout.fruchterman.reingold, pal = brewer.pal(7,"Set2"), random = TRUE, node.pies = TRUE, pie.local = TRUE, vertex.radius = 0.03, scale.vertices = 0.05, edge.color = NULL, vshape = "none", vsize = 15, ewidth = 3, margin = 0, vlabel.cex = 0.8, vlabel.color = "black", vlabel.family = "Helvetica", vertex.color = "palegoldenrod", vlabel = TRUE, col.nonclusters = "black", jitter = 0.2, circle = TRUE, printids = TRUE, cid.cex = 1, shownodesin = 0, showall = FALSE, verbose = TRUE, ...)
  # x is a "linkcomm" object.
{
  if(length(nodes) > 0){
    clusterids <- which.communities(x, nodes = nodes)
  }
  clusters <- x$clusters[clusterids]
  miss <- setdiff(x$hclust$order,unlist(clusters))
  crf <- colorRampPalette(pal,bias=1)
  cols <- crf(length(clusters))
  if(random){
    cols <- sample(cols,length(clusters),replace=FALSE)
  }
  if(showall){
    # Add single edge "clusters".
    single <- setdiff(1:x$numbers[1],unlist(clusters))
    ll <- length(clusters)
    for(i in 1:length(single)){
      clusters[[(i+ll)]] <- single[i]
    }
    cols <- append(cols, rep(col.nonclusters, length(single)))
  }
  drawcircle <- FALSE
  if(inherits(layout,"character")){
    if(layout == "spencer.circle"){
      if(length(clusters) > length(x$clusters[1:x$numbers[3]])){
        clusterids <- 1:x$numbers[3]
      }
      ord <- orderCommunities(x, clusterids = clusterids, verbose = FALSE)
      clusters <- ord$ordered
      clusterids <- ord$clusids
      layout <- layout.spencer.circle(x, clusterids = clusterids, jitter = jitter, verbose = verbose)$nodes
      drawcircle <- TRUE
    }
  }
  names(cols) <- clusterids
  if(length(unlist(clusters)) < nrow(x$edgelist) || length(miss) == 0){
    # Convert old clus ids into new ones.
    edges <- x$edgelist[unlist(clusters),]
    ig <- graph.edgelist(edges, directed=x$directed)
    clen <- sapply(clusters,length)
    j<-1
    # Colour edges according to community membership.
    for(i in 1:length(clusters)){
      newcids <- j:sum(clen[1:i])
      E(ig)[newcids]$color <- cols[i]
      j <- tail(newcids,1)+1
    }
  }else{
    ig <- x$igraph
    for(i in 1:length(clusters)){
      E(ig)[clusters[[i]]]$color <- cols[i]
    }
  }
  
  if(shownodesin == 0){
    vnames <- V(ig)$name
  }else{ # Show nodes that belong to more than x number of communities.
    vnames <- V(ig)$name
    inds <- NULL
    for(i in 1:length(vnames)){
      if(x$numclusters[which(names(x$numclusters)==vnames[i])] < shownodesin){
        inds <- append(inds,i)
      }
    }
    vnames[inds] <- ""
  }
  if(vlabel==FALSE){
    vnames = NA
  }
  
  dev.hold(); on.exit(dev.flush())
  oldpar <- par(no.readonly = TRUE)
  par(mar = c(4,4,2,2))
  
  if(!node.pies){
    plot(ig, layout=layout, vertex.shape=vshape, edge.width=ewidth, vertex.label=vnames, vertex.label.family=vlabel.family, vertex.label.color=vlabel.color, vertex.size=vsize, vertex.color=vertex.color, margin=margin, vertex.label.cex = vlabel.cex, ...)
  }else{
    nodes <- V(ig)$name
    # Get node community membership by edges.
    if(pie.local){
      edge.memb <- numberEdgesIn(x, clusterids = clusterids, nodes = nodes)
    }else{
      edge.memb <- numberEdgesIn(x, nodes = nodes)
    }
    
    cat("   Getting node layout...")
    if(inherits(layout,"function")){
      lay <- layout(ig)
    }else{
      lay <- layout
    }
    lay <- layout.norm(lay, xmin=-1, xmax=1, ymin=-1, ymax=1)
    rownames(lay) <- V(ig)$name
    cat("\n")
    node.pies <- .nodePie(edge.memb=edge.memb, layout=lay, nodes=nodes, edges=100, radius=vertex.radius, scale=scale.vertices)
    cat("\n")
    # Plot graph.
    if(is.null(edge.color)){
      plot(ig, layout=lay, vertex.shape="none", vertex.label=NA, vertex.label.dist=1, edge.width=ewidth, vertex.label.color=vlabel.color, ...)
    }else{
      plot(ig, layout=lay, vertex.shape="none", vertex.label=NA, vertex.label.dist=1, edge.width=ewidth, vertex.label.color=vlabel.color, edge.color=edge.color, ...)
    }
    labels <- list()
    # Plot node pies and node names.
    for(i in 1:length(node.pies)){
      yp <- NULL
      for(j in 1:length(node.pies[[i]])){
        seg.col <- cols[which(names(cols)==names(edge.memb[[i]])[j])]
        polygon(node.pies[[i]][[j]][,1], node.pies[[i]][[j]][,2], col = seg.col)
        yp <- append(yp, node.pies[[i]][[j]][,2])
      }
      lx <- lay[which(rownames(lay)==names(node.pies[i])),1] + 0.1
      ly <- max(yp) + 0.02 # Highest point of node pie.
      labels[[i]] <- c(lx, ly)
    }
    # Plot node names after nodes so they overlay them.
    for(i in 1:length(labels)){
      text(labels[[i]][1], labels[[i]][2], labels = vnames[which(nodes==names(node.pies[i]))], cex = vlabel.cex, col = vlabel.color)
    }
  }
  
  if(circle && drawcircle){
    # Add circle for Spencer layout.
    cx<-NULL; for(i in 1:100){cx[i]<-1.25*cos(i*(2*pi)/100)}
    cy<-NULL; for(i in 1:100){cy[i]<-1.25*sin(i*(2*pi)/100)}
    polygon(cx-0.08,cy-0.08, border="grey",lwd=2)
    # Add community anchor points and cluster IDs.
    for(i in 1:length(clusters)){
      px <- 1.1*cos(i*(2*pi)/length(clusters))
      py <- 1.1*sin(i*(2*pi)/length(clusters))
      points(px-0.08,py-0.08, pch = 20, col = cols[i])
      if(printids){
        tx <- 1.3*cos(i*(2*pi)/length(clusters))
        ty <- 1.3*sin(i*(2*pi)/length(clusters))
        text(tx-0.08,ty-0.08, labels = clusterids[i], col = cols[i], cex = cid.cex, font=2)
      }
    }
  }
  par(oldpar)
  
}


## Internal function layout.spencer.circle from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
layout.spencer.circle <- function(x, clusterids = 1:x$numbers[3], verbose = TRUE, jitter = 0.2)
  # Returns x-y node coordinates for Rob Spencer's circular layout of link communities together with x-y coordinates for the community anchors.
  # x is a "linkcomm" object.
{
  clusters <- x$clusters[clusterids]
  edges <- x$edgelist[unlist(clusters),]
  ig <- graph.edgelist(edges, directed=FALSE)
  # Put communities in dendrogram order.
  clusters <- orderCommunities(x, clusterids = clusterids, verbose = verbose)$ordered
  # Set up community anchor points in Cartesian coordinates around unit circle (communities evenly spaced).
  xy_anchors <- matrix(0,length(clusters),2)
  for(i in 1:length(clusters)){
    xy_anchors[i,] <- c(cos(i*(2*pi)/length(clusters)), sin(i*(2*pi)/length(clusters)))
  }
  # Calculate community membership percentages per node.
  nodes <- c(x$edgelist[unlist(clusters),1],x$edgelist[unlist(clusters),2])
  node_names <- unique(nodes)
  xy_nodes <- matrix(0,length(node_names),2)
  for(i in 1:length(node_names)){
    if(verbose){
      mes <- paste(c("   Calculating node co-ordinates for Spencer circle...",floor(i/(length(node_names))*100),"%"),collapse="")
      cat(mes,'\r')
      flush.console()
    }
    freqs <- NULL
    total <- length(which(nodes==node_names[i]))
    for(j in 1:length(clusters)){
      freqs <- length(which(c(x$edgelist[clusters[[j]],1],x$edgelist[clusters[[j]],2])==node_names[i]))/total
      # Update x-y coordinates for this node.
      xy_nodes[i,1] <- sum(xy_nodes[i,1], freqs*xy_anchors[j,1])
      xy_nodes[i,2] <- sum(xy_nodes[i,2], freqs*xy_anchors[j,2])
    }
    # Add random jitter if this node has identical x-y coordinates to an earlier node.
    if(duplicated(xy_nodes)[i]){
      xy_nodes[i,1] <- sum(xy_nodes[i,1], runif(1, min = -jitter, max = jitter))
      xy_nodes[i,2] <- sum(xy_nodes[i,2], runif(1, min = -jitter, max = jitter))
    }
  }
  
  rownames(xy_nodes) <- node_names
  
  xy_nodes <- xy_nodes[match(V(ig)$name,rownames(xy_nodes)),]
  xy_nodes <- xy_nodes[!is.na(xy_nodes[,1]),]
  
  if(verbose){cat("\n")}
  
  xy <- list()
  xy$nodes <- xy_nodes
  xy$anchors <- xy_anchors
  
  return(xy)
  
}

## Internal function plotLinkCommSummComm from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plotLinkCommSummComm <- function(x, clusterids = 1:x$numbers[3], summary = "conn", pie = FALSE, col = TRUE, pal = brewer.pal(11,"Spectral"), random = FALSE, verbose = TRUE, ...) 
  # Plots pie or bar chart summarising sizes of communities in terms of nodes, link density, community connectedness, or community modularity.
  # x is a "linkcomm" object.
{
  if(col){
    crf <- colorRampPalette(pal,bias=1)
    cols <- crf(length(clusterids))
    if(random){
      cols <- sample(cols,length(clusterids),replace=FALSE)
    }
  }else{
    cols <- "lightblue"
  }
  # Extract number of nodes per community.
  nums <- NULL
  if(summary == "nodes"){
    for(i in 1:length(clusterids)){
      nums[i] <- length(unique(c(x$edgelist[x$clusters[[clusterids[i]]],1],x$edgelist[x$clusters[[clusterids[i]]],2])))
    }
    main <- "Node density per community"
  }else if(summary == "ld"){
    nums <- LinkDensities(x, clusterids = clusterids)
    main <- "Link density per community"
  }else{
    nums <- getCommunityConnectedness(x, clusterids = clusterids, conn = summary, verbose = verbose)
    if(summary == "conn"){
      main <- "Community Connectedness"
    }else{
      main <- "Community Modularity"
    }
  }
  names(nums) <- clusterids
  
  if(pie){
    pie(nums, col = cols, main = main, ...)
  }else{
    barplot(nums, xlab = "Community", ylab = main, col = cols)
    abline(h=0)
  }
  
}

## Internal function CLIQUE from subspace package
## (https://github.com/cran/subspace).
# The CLIQUE Algorithm for Subspace Clustering
CLIQUE  <- function(data,xi=10,tau=0.2) {
  arr <- java_object_from_data(data)
  #Now that the data is in the correct format, we can call into our Java Code that will then call into the
  #actual implementation of the CLIQUE Algorithm
  res <- rJava::.jcall("ClusteringApplier",returnSig="[Li9/subspace/base/Cluster;",method="clique",arr,as.integer(xi),tau,evalArray=F)
  #We can then turn the Java Clustering Object that was returned into an R-Friendly S3-Object
  res <- r_clusters_from_java_clusters(res)
  return(res)
}

## Internal function r_clusters_from_java_clusters from subspace package
## (https://github.com/cran/subspace).
#This function turns Cluster[] from Java into a list of more suitable S3 Objects.
r_clusters_from_java_clusters <- function(clus) {
  
  if(rJava::is.jnull(clus)) {
    warning("An error occured in the clustering function. Therefore, NULL is returned.")
    return(NULL)
  }
  
  subspace_matrix <- rJava::.jcall("ClusteringApplier",returnSig="[[Z",method="extract_subspace",clus,simplify=T)
  objects_matrix <- rJava::.jcall("ClusteringApplier",returnSig="[[I",method="extract_objects",clus,simplify=T)
  
  if(nrow(subspace_matrix)==0){
    warning("No subspace Clusters were generated. NULL is being returned. This is probably due to the parameters given to the clustering algorithm. Try a set of parameters that is more likely to produce many clusters")
    return(NULL)
  }
  
  res <- lapply(1:nrow(subspace_matrix),function(index){
    objects <- as.vector(objects_matrix[index,])
    #Add 1 because Java uses 0 as first index but R uses 1
    objects <- objects+1
    #Filter out those indices that were added in "extract_objects" to make the objects matrix rectangular
    objects <- objects[objects>0]
    return(subspace_cluster(subspace=as.vector(subspace_matrix[index,]),
                            objects=objects))
  })
  class(res) <- append(class(res),"subspace_clustering")
  return(res)
}

## Internal function java_object_from_data from subspace package
## (https://github.com/cran/subspace).
#Turns the input data that is in the format of a matrix or data frame into a Reference to a 
#java double[][] Object.
java_object_from_data <- function(data) {
  #To achieve this, the input data matrix (or data frame) is first turned into a vector column by column. This
  #Vector is then passed into a java function that also needs to know the number of columns of the original
  #matrix to reconstruct it as a double[][]. This is much faster than producing a double[][] with .jarray.
  res <- rJava::.jcall("JavaObjectFromDataConverter",returnSig="[[D",method="matrix_from_array",
                       as.vector(as.matrix(data)),ncol(data),evalArray=F)
  return(res)
}

## Internal function subspace_cluster from subspace package
## (https://github.com/cran/subspace).
subspace_cluster <- function(subspace,objects) {
  res <- list(subspace=subspace,objects=objects)
  class(res) <- append(class(res),"subspace_cluster")
  return(res)
}

## Internal function CandModgenes for extracting candidate module genes
CandModgenes <- function(ceRExp, 
                         mRExp = NULL, 
                         Modulegenes, 
                         num.ModuleceRs = 2, 
                         num.ModulemRs = 2){
    
    if(is.null(mRExp)){
      ceR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(unique(Modulegenes[[i]]) %in%
                                                                           colnames(ceRExp))))
      index <- which(ceR_Num >= num.ModuleceRs)
      CandidateModulegenes <- lapply(index, function(i) unique(Modulegenes[[i]]))
    } else {
      ceR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
                                                                           colnames(ceRExp))))
      mR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
                                                                          colnames(mRExp))))
      index <- which(ceR_Num >= num.ModuleceRs & mR_Num >= num.ModulemRs)
      CandidateModulegenes <- lapply(index, function(i) Modulegenes[[i]]) 
    }
    
    CandidateModulegenes <- lapply(seq_along(index), function(i) GeneSet(CandidateModulegenes[[i]], 
        setName = paste("Module", i, sep=" ")))
    CandidateModulegenes <- GeneSetCollection(CandidateModulegenes)
    
    return(CandidateModulegenes)
}

## Internal function module_group_sim_matrix for Calculating similarity matrix between two list of module groups
module_group_sim_matrix <- function(Module.group1, 
                                    Module.group2, 
                                    sim.method = "Simpson"){
  
  m <- length(Module.group1)
  n <- length(Module.group2)
  Sim <- matrix(NA, m, n)
  
  if (sim.method == "Simpson") {
    for (i in seq(m)){
      for (j in seq(n)){	    
        overlap_vertex <- length(intersect(Module.group1[[i]], Module.group2[[j]]))
        min_vertex <- min(length(Module.group1[[i]]), length(Module.group2[[j]]))
        vertex_Sim <- overlap_vertex/min_vertex           
        Sim[i, j] <- vertex_Sim
      }
    }
  } else if (sim.method == "Jaccard") {
    for (i in seq(m)){
      for (j in seq(n)){	    
        overlap_vertex <- length(intersect(Module.group1[[i]], Module.group2[[j]]))
        union_vertex <- length(union(Module.group1[[i]], Module.group2[[j]]))
        vertex_Sim <- overlap_vertex/union_vertex
        Sim[i, j] <- vertex_Sim
      }
    }
  } else if (sim.method == "Lin") {
    for (i in seq(m)){
      for (j in seq(n)){
        overlap_vertex <- length(intersect(Module.group1[[i]], Module.group2[[j]]))
        sum_vertex <- length(Module.group1[[i]]) + length(Module.group2[[j]])
        vertex_Sim <- 2 * overlap_vertex/sum_vertex           
        Sim[i, j] <- vertex_Sim
      }
    }
  }    
  
  return(Sim)
}

## Internal function cluster from miRspongeR package
## Disease enrichment analysis of modules
moduleDEA <- function(Modulelist, OrgDb = "org.Hs.eg.db", padjustvaluecutoff = 0.05,
                      padjustedmethod = "BH") {
  
  entrezIDs <- lapply(seq_along(Modulelist), function(i) bitr(Modulelist[[i]], fromType = "SYMBOL",
                                                              toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID)
  
  entrezIDs <- lapply(seq_along(Modulelist), function(i) as.character(entrezIDs[[i]]))
  
  enrichDOs <- lapply(seq_along(Modulelist), function(i) enrichDO(entrezIDs[[i]], pvalueCutoff = padjustvaluecutoff,
                                                                  pAdjustMethod = padjustedmethod))
  
  enrichDGNs <- lapply(seq_along(Modulelist), function(i) enrichDGN(entrezIDs[[i]], pvalueCutoff = padjustvaluecutoff,
                                                                    pAdjustMethod = padjustedmethod))
  
  enrichNCGs <- lapply(seq_along(Modulelist), function(i) enrichNCG(entrezIDs[[i]], pvalueCutoff = padjustvaluecutoff,
                                                                    pAdjustMethod = padjustedmethod))
  
  return(list(enrichDOs, enrichDGNs, enrichNCGs))
}

## Internal function cluster from miRspongeR package
## Functional GO, KEGG and Reactome enrichment analysis of modules
moduleFEA <- function(Modulelist, ont = "BP", KEGGorganism = "hsa", Reactomeorganism = "human",
                      OrgDb = "org.Hs.eg.db", padjustvaluecutoff = 0.05, padjustedmethod = "BH") {
  
  entrezIDs <- lapply(seq_along(Modulelist), function(i) bitr(Modulelist[[i]], fromType = "SYMBOL",
                                                              toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID)
  
  entrezIDs <- lapply(seq_along(Modulelist), function(i) as.character(entrezIDs[[i]]))
  
  enrichGOs <- lapply(seq_along(Modulelist), function(i) enrichGO(entrezIDs[[i]], OrgDb = OrgDb,
                                                                  ont = ont, pvalueCutoff = padjustvaluecutoff, pAdjustMethod = padjustedmethod))
  
  enrichKEGGs <- lapply(seq_along(Modulelist), function(i) enrichKEGG(entrezIDs[[i]], organism = KEGGorganism,
                                                                      pvalueCutoff = padjustvaluecutoff, pAdjustMethod = padjustedmethod))
  
  enrichReactomes <- lapply(seq_along(Modulelist), function(i) enrichPathway(entrezIDs[[i]], organism = Reactomeorganism,
                                                                             pvalueCutoff = padjustvaluecutoff, pAdjustMethod = padjustedmethod))
  
  return(list(enrichGOs, enrichKEGGs, enrichReactomes))
  
}

#' Identification of co-expressed gene modules from matched ceRNA and mRNA
#' expression data or single gene expression data using WGCNA package
#'
#' @title module_WGCNA
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param RsquaredCut Desired minimum scale free topology fitting index 
#' R^2 with interval [0 1].
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom WGCNA pickSoftThreshold
#' @importFrom WGCNA adjacency
#' @importFrom WGCNA TOMdist
#' @importFrom WGCNA standardColors
#' @importFrom flashClust flashClust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp[, seq_len(80)], 
#'     mRExp[, seq_len(80)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Langfelder P, Horvath S. WGCNA: an R package for weighted 
#' correlation network analysis. BMC Bioinformatics. 2008, 9:559.#' 
module_WGCNA <- function(ceRExp, 
                         mRExp = NULL, 
                         RsquaredCut = 0.9, 
                         num.ModuleceRs = 2,
                         num.ModulemRs = 2) {

    if(is.null(mRExp)){
      ExpData <- assay(ceRExp)  
    } else {
      ExpData <- cbind(assay(ceRExp), assay(mRExp))
    }
  
    Optimalpower <- pickSoftThreshold(ExpData, RsquaredCut = RsquaredCut)$powerEstimate
    adjacencymatrix <- adjacency(ExpData, power = Optimalpower)
    dissTOM <- TOMdist(adjacencymatrix)
    hierTOM <- flashClust(as.dist(dissTOM), method = "average")

    # The function cutreeDynamic colors each gene by the branches that
    # result from choosing a particular height cutoff.
    colorh <- cutreeDynamic(hierTOM, method = "tree") + 1
    StandColor <- c("grey", standardColors(n = NULL))
    colorh <- unlist(lapply(seq_len(length(colorh)), function(i) StandColor[colorh[i]]))
    colorlevels <- unique(colorh)
    colorlevels <- colorlevels[-which(colorlevels == "grey")]

    Modulegenes <- lapply(seq_len(length(colorlevels)), function(i) colnames(ExpData)[which(colorh ==
        colorlevels[i])])
    
    CandidateModulegenes <- CandModgenes(ceRExp, mRExp = mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}


#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data or single gene expression data using GFA package
#'
#' @title module_GFA
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param StrengthCut Desired minimum strength (absolute value of 
#' association with interval [0 1]) for each bicluster.
#' @param iter.max The total number of Gibbs sampling steps 
#' (default 1000).
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom GFA normalizeData
#' @importFrom GFA getDefaultOpts
#' @importFrom GFA gfa
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_GFA <- module_GFA(ceRExp[seq_len(20), seq_len(15)],
#'     mRExp[seq_len(20), seq_len(15)], iter.max = 3000)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Bunte K, Lepp\'{a}aho E, Saarinen I, Kaski S. 
#' Sparse group factor analysis for biclustering of multiple data sources. Bioinformatics. 2016, 32(16):2457-63.
#' @references Lepp\'{a}aho E, Ammad-ud-din M, Kaski S. GFA: 
#' exploratory analysis of multiple data sources with group factor 
#' analysis. J Mach Learn Res. 2017, 18(39):1-5.
module_GFA <- function(ceRExp, 
                       mRExp = NULL, 
                       StrengthCut = 0.9, 
                       iter.max = 5000,
                       num.ModuleceRs = 2, 
                       num.ModulemRs = 2) {
    
    if(is.null(mRExp)){
    ExpData <- list(assay(ceRExp))
    names(ExpData) = c("ceRNA expression")
    } else {
    ExpData <- list(assay(ceRExp), assay(mRExp))
    names(ExpData) = c("ceRNA expression", "mRNA expression")
    }
  
    # Normalize the data - here we assume that every feature is equally
    # important
    norm <- normalizeData(ExpData, type = "scaleFeatures")

    # Get the model options to detect bicluster structure
    opts <- getDefaultOpts(bicluster = TRUE)

    # Check for sampling chain convergence
    opts$convergenceCheck <- TRUE
    opts$iter.max <- iter.max

    # Infer the model
    res <- gfa(norm$train, opts = opts)

    # Extract gene index of each bicluster, using stength cutoff (absolute
    # value of association)
    BCresnum <- lapply(seq_len(dim(res$W)[2]), function(i) which(abs(res$W[,
        i]) >= StrengthCut))

    # Extract genes of each bicluster
    if(is.null(mRExp)){
    Modulegenes <- lapply(seq_along(BCresnum), function(i) colnames(assay(ceRExp))[BCresnum[[i]]])  
    } else {
    Modulegenes <- lapply(seq_along(BCresnum), function(i) colnames(cbind(assay(ceRExp),
        assay(mRExp)))[BCresnum[[i]]])
    }
    
    CandidateModulegenes <- CandModgenes(ceRExp, mRExp = mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}


#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data or single gene expression data using igraph package
#'
#' @title module_igraph
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param cor.method The method of calculating correlation selected, 
#' including 'pearson' (default), 'kendall', 'spearman'.
#' @param pos.p.value.cutoff The significant p-value cutoff of 
#' positive correlation.
#' @param cluster.method The clustering method selected in 
#' \pkg{igraph} package, including 'betweenness', 'greedy' (default), 
#' 'infomap', 'prop', 'eigen', 'louvain', 'walktrap'.
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom igraph graph_from_biadjacency_matrix
#' @importFrom igraph cluster_edge_betweenness
#' @importFrom igraph cluster_fast_greedy
#' @importFrom igraph cluster_infomap
#' @importFrom igraph cluster_label_prop
#' @importFrom igraph cluster_leading_eigen
#' @importFrom igraph cluster_louvain
#' @importFrom igraph cluster_walktrap
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_igraph <- module_igraph(ceRExp[, seq_len(10)],
#'     mRExp[, seq_len(10)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Csardi G, Nepusz T. The igraph software package for 
#' complex network research, InterJournal, Complex Systems. 2006:1695.
module_igraph <- function(ceRExp, 
                          mRExp = NULL, 
                          cor.method = "pearson", 
                          pos.p.value.cutoff = 0.01,
                          cluster.method = "greedy", 
                          num.ModuleceRs = 2, 
                          num.ModulemRs = 2) {

    cor.binary <- cor_binary(ceRExp, mRExp = mRExp, cor.method = cor.method, pos.p.value.cutoff = pos.p.value.cutoff)
    cor.binary.graph <- graph_from_biadjacency_matrix(cor.binary)

    if (cluster.method == "betweenness") {
        Modulegenes <- cluster_edge_betweenness(cor.binary.graph)
    } else if (cluster.method == "greedy") {
        Modulegenes <- cluster_fast_greedy(cor.binary.graph)
    } else if (cluster.method == "infomap") {
        Modulegenes <- cluster_infomap(cor.binary.graph)
    } else if (cluster.method == "prop") {
        Modulegenes <- cluster_label_prop(cor.binary.graph)
    } else if (cluster.method == "eigen") {
        Modulegenes <- cluster_leading_eigen(cor.binary.graph)
    } else if (cluster.method == "louvain") {
        Modulegenes <- cluster_louvain(cor.binary.graph)
    } else if (cluster.method == "walktrap") {
        Modulegenes <- cluster_walktrap(cor.binary.graph)
    }

    CandidateModulegenes <- CandModgenes(ceRExp, mRExp = mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}


#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data or single gene expression data using ProNet package
#'
#' @title module_ProNet
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param cor.method The method of calculating correlation selected, 
#' including 'pearson' (default), 'kendall', 'spearman'.
#' @param pos.p.value.cutoff The significant p-value cutoff of 
#' positive correlation
#' @param cluster.method The clustering method selected in 
#' \pkg{ProNet} package, including 'FN', 'MCL' (default), 
#' 'LINKCOMM', 'MCODE'.
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @import igraph
#' @import RColorBrewer
#' @import grid
#' @importFrom Rcpp evalCpp
#' @importFrom MCL mcl
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.flush
#' @importFrom grDevices dev.hold
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_ProNet <- module_ProNet(ceRExp[, seq_len(10)],
#'     mRExp[, seq_len(10)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Clauset A, Newman ME, Moore C. Finding community 
#' structure in very large networks. Phys Rev E Stat Nonlin Soft 
#' Matter Phys., 2004, 70(6 Pt 2):066111.
#' @references Enright AJ, Van Dongen S, Ouzounis CA. An efficient 
#' algorithm for large-scale detection of protein families. 
#' Nucleic Acids Res., 2002, 30(7):1575-84.
#' @references Kalinka AT, Tomancak P. linkcomm: an R package 
#' for the generation, visualization, and analysis of link 
#' communities in networks of arbitrary size and type. 
#' Bioinformatics, 2011, 27(14):2011-2.
#' @references Bader GD, Hogue CW. An automated method for 
#' finding molecular complexes in large protein interaction 
#' networks. BMC Bioinformatics, 2003, 4:2.
module_ProNet <- function(ceRExp, 
                          mRExp = NULL, 
                          cor.method = "pearson", 
                          pos.p.value.cutoff = 0.01,
                          cluster.method = "MCL", 
                          num.ModuleceRs = 2, 
                          num.ModulemRs = 2) {

    cor.binary <- cor_binary(ceRExp, mRExp = mRExp, cor.method = cor.method, pos.p.value.cutoff = pos.p.value.cutoff)
    cor.binary.graph <- graph_from_biadjacency_matrix(cor.binary)

    if (cluster.method == "FN" | cluster.method == "MCL") {
        network_Cluster <- cluster(cor.binary.graph, method = cluster.method)

        Modulegenes <- lapply(seq_len(max(network_Cluster)), function(i) rownames(as.matrix(network_Cluster))[which(network_Cluster ==
            i)])
    } else if (cluster.method == "LINKCOMM") {
        edgelist <- get.edgelist(cor.binary.graph)
        network_Cluster <- getLinkCommunities(edgelist)$nodeclusters
        Modulegenes <- lapply(seq_len(max(c(network_Cluster$cluster))),
            function(i) as.character(network_Cluster$node[which(c(network_Cluster$cluster) ==
                i)]))
    } else if (cluster.method == "MCODE") {
        network_Cluster <- cluster(cor.binary.graph, method = cluster.method) +
            1
        Modulegenes <- lapply(seq_len(max(network_Cluster)), function(i) rownames(as.matrix(network_Cluster))[which(network_Cluster ==
            i)])
    }

    CandidateModulegenes <- CandModgenes(ceRExp, mRExp = mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}


#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data or single gene expression data using NMF package
#'
#' @title module_NMF
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param NMF.algorithm Specification of the NMF algorithm, 
#' including 'brunet' (default), 'Frobenius', 'KL', 'lee', 'nsNMF', 
#' 'offset', 'siNMF', 'snmf/l', 'snmf/r'.
#' @param num.modules The number of modules to be identified.
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom NMF nmf
#' @importFrom NMF predict
#' @importFrom NMF nneg
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' # Reimport NMF package to avoid conflicts with DelayedArray package
#' library(NMF)
#' modulegenes_NMF <- module_NMF(ceRExp[, seq_len(10)],
#'     mRExp[, seq_len(10)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references  Gaujoux R, Seoighe C. A flexible R package for 
#' nonnegative matrix factorization. BMC Bioinformatics. 2010, 11:367.
module_NMF <- function(ceRExp, 
                       mRExp = NULL, 
                       NMF.algorithm = "brunet", 
                       num.modules = 10,
                       num.ModuleceRs = 2, 
                       num.ModulemRs = 2) {

    if(is.null(mRExp)){
    ExpData <- assay(ceRExp)  
    } else {
    ExpData <- cbind(assay(ceRExp), assay(mRExp))
    }

    # Run NMF algorithm with rank num.modules, negative values are transformed
    # into 0 if exist in expression data
    res <- nmf(nneg(ExpData), rank = num.modules, method = NMF.algorithm)

    # Predict column clusters
    Cluster.membership <- predict(res)

    # Extract genes of each cluster
    Modulegenes <- lapply(seq_len(num.modules), function(i) colnames(ExpData)[which(Cluster.membership ==
        i)])

    CandidateModulegenes <- CandModgenes(ceRExp, mRExp = mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}

#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data or single gene expression data using a series of clustering packages, 
#' including stats, flashClust, dbscan, subspace, mclust, SOMbrero and ppclust packages.
#' 
#' @title module_clust 
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param cluster.method Specification of the clustering method, 
#' including 'kmeans'(default), 'hclust', 'dbscan' , 'clique', 
#' 'gmm', 'som' and 'fcm'.
#' @param num.modules Parameter of the number of modules to be identified
#' for the 'kmeans', 'hclust', 'gmm' and 'fcm' methods. Parameter of the number
#' of intervals for the 'clique' method. For the 'dbscan' and 'som' methods,
#' no need to set the parameter.
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom rJava .jcall
#' @importFrom stats kmeans
#' @importFrom stats dist
#' @importFrom stats cutree
#' @importFrom flashClust flashClust
#' @importFrom dbscan optics
#' @importFrom dbscan dbscan
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom SOMbrero trainSOM
#' @importFrom ppclust fcm
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_clust <- module_clust(ceRExp[, seq_len(30)],
#'     mRExp[, seq_len(30)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Forgy EW. Cluster analysis of multivariate 
#' data: efficiency vs interpretability of classifications. 
#' Biometrics, 1965, 21:768-769.
#' @references Hartigan JA, Wong MA. 
#' Algorithm AS 136: A K-means clustering algorithm. 
#' Applied Statistics, 1979, 28:100-108.
#' @references Lloyd SP. Least squares quantization in PCM. 
#' Technical Note, Bell Laboratories. Published in 1982 
#' in IEEE Transactions on Information Theory, 1982, 28:128-137.
#' @references MacQueen J. Some methods for classification 
#' and analysis of multivariate observations. 
#' In Proceedings of the Fifth Berkeley Symposium on 
#' Mathematical Statistics and Probability, 
#' eds L. M. Le Cam & J. Neyman, 1967, 1, pp.281-297. 
#' Berkeley, CA: University of California Press.
#' @references Langfelder P, Horvath S. Fast R Functions for 
#' Robust Correlations and Hierarchical Clustering. 
#' Journal of Statistical Software. 2012, 46(11):1-17.
#' @references Ester M, Kriegel HP, Sander J, Xu X. A density-based 
#' algorithm for discovering clusters in large spatial databases with 
#' noise, Proceedings of 2nd International Conference on Knowledge Discovery and
#' Data Mining (KDD-96), 1996, 96(34): 226-231.
#' @references Campello RJGB, Moulavi D, Sander J. 
#' Density-based clustering based on hierarchical density estimates,
#' Pacific-Asia conference on knowledge discovery and data mining. 
#' Springer, Berlin, Heidelberg, 2013: 160-172.
#' @references Agrawal R, Gehrke J, Gunopulos D, Raghavan P. 
#' Automatic subspace clustering of high dimensional data for 
#' data mining applications. In Proc. ACM SIGMOD, 1998.
#' @references Scrucca L, Fop M, Murphy TB, Raftery AE. 
#' mclust 5: clustering, classification and density estimation using 
#' Gaussian finite mixture models The R Journal 8/1, 2016, pp. 205-233.
#' @references Kohonen T. Self-Organizing Maps. 
#' Berlin/Heidelberg: Springer-Verlag, 3rd edition, 2001. 
#' @references Dunn JC. A fuzzy relative of the ISODATA process 
#' and its use in detecting compact well-separated clusters. Journal of Cybernetics, 
#' 1973, 3(3):32-57.
#' @references Bezdek JC. Cluster validity with fuzzy sets. Journal of Cybernetics, 1974, 3: 58-73.
#' @references Bezdek JC. Pattern recognition with fuzzy objective function 
#' algorithms. Plenum, NY, 1981. 
module_clust <- function(ceRExp, 
                         mRExp = NULL, 
                         cluster.method = "kmeans", 
                         num.modules = 10,
                         num.ModuleceRs = 2, 
                         num.ModulemRs = 2) {
  
  if(is.null(mRExp)){
    ExpData <- assay(ceRExp)  
  } else {
    ExpData <- cbind(assay(ceRExp), assay(mRExp))
  }
  
  if (cluster.method == "kmeans") {
    res <- kmeans(t(ExpData), centers = num.modules, iter.max = 100)
  } else if (cluster.method == "hclust") {
    diss <- dist(t(ExpData))
    hc <- flashClust(diss, method = "average")
    res <- cutree(hc, k = num.modules)
  } else if (cluster.method == "dbscan") {
    eps <- optics(t(ExpData))$eps
    res <- dbscan(t(ExpData), eps = eps)
  } else if (cluster.method == "clique") {
    res <- CLIQUE(t(ExpData), xi = num.modules)
  } else if (cluster.method == "gmm") {
    res <- Mclust(t(ExpData), G = num.modules)$classification
  } else if (cluster.method == "som") {
    res <- trainSOM(t(ExpData))$clustering
  } else if (cluster.method == "fcm") {
    res <- fcm(t(ExpData), centers = num.modules)$cluster
  } 
  
  # Extract genes of each cluster
  if (cluster.method == "kmeans") {
    Cluster.membership <- res$cluster
    Modulegenes <- lapply(seq_len(num.modules), function(i) 
      colnames(ExpData)[which(Cluster.membership == i)])
  }
  
  if (cluster.method == "hclust" | cluster.method == "gmm") {
    Cluster.membership <- res
    Modulegenes <- lapply(seq_len(num.modules), function(i) 
      names(res)[which(Cluster.membership == i)])
  }
  
  if (cluster.method == "dbscan" ) {
    Cluster.membership <- res$cluster
    Modulegenes <- lapply(seq_len(max(Cluster.membership)), function(i) 
      colnames(ExpData)[which(Cluster.membership == i)])
  }
  
  if (cluster.method == "clique") {
    Modulegenes <- lapply(seq(length(res)), function(i) 
      colnames(ExpData)[res[[i]]$objects])
  }
  
  if (cluster.method == "som" ) {
    Cluster.membership <- res
    Modulegenes <- lapply(seq_len(max(Cluster.membership)), function(i) 
      names(res)[which(Cluster.membership == i)])
  }
  
  if (cluster.method == "fcm" ) {
    Cluster.membership <- res
    Modulegenes <- lapply(seq_len(num.modules), function(i) 
      colnames(ExpData)[which(Cluster.membership == i)])
  }
                                                                                          
  
  CandidateModulegenes <- CandModgenes(ceRExp, mRExp = mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
                                       num.ModulemRs = num.ModulemRs)
  
  return(CandidateModulegenes)
}


#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data or single gene expression data using a series  
#' of biclustering packages, including fabia,  
#' BicARE and isa2
#'
#' @title module_biclust
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param BCmethod Specification of the biclustering method, 
#' including 'fabia' (default), 'fabiap', 
#' 'fabias', 'mfsc', 'nmfdiv', 'nmfeu', 'nmfsc', 'FLOC' and 'isa'.
#' @param num.modules The number of modules to be identified. For the 
#' 'isa' method, no need to set the parameter.  
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @import methods
#' @importFrom fabia fabia
#' @importFrom fabia fabiap
#' @importFrom fabia fabias
#' @importFrom fabia mfsc
#' @importFrom fabia nmfdiv
#' @importFrom fabia nmfeu
#' @importFrom fabia nmfsc
#' @importFrom fabia extractBic
#' @importFrom BicARE FLOC
#' @importFrom BicARE bicluster
#' @importFrom isa2 isa
#' @importFrom isa2 isa.biclust
#' @importFrom Biobase ExpressionSet
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_biclust <- module_biclust(ceRExp[, seq_len(30)],
#'     mRExp[, seq_len(30)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Hochreiter S, Bodenhofer U, Heusel M, Mayr A, 
#' Mitterecker A, Kasim A, Khamiakova T, Van Sanden S, Lin D, 
#' Talloen W, Bijnens L, G\'{o}hlmann HW, Shkedy Z, Clevert DA. 
#' FABIA: factor analysis for bicluster acquisition. 
#' Bioinformatics. 2010, 26(12):1520-7.
#' @references Yang J, Wang H, Wang W, Yu, PS. An improved 
#' biclustering method for analyzing gene expression. 
#' Int J Artif Intell Tools. 2005, 14(5): 771-789.
#' @references Bergmann S, Ihmels J, Barkai N. Iterative 
#' signature algorithm for the analysis of large-scale gene 
#' expression data. Phys Rev E Stat Nonlin Soft Matter Phys. 
#' 2003, 67(3 Pt 1):031902.
#' @references Sill M, Kaiser S, Benner A, Kopp-Schneider A. 
#' Robust biclustering by sparse singular value decomposition 
#' incorporating stability selection. Bioinformatics. 2011, 
#' 27(15):2089-97.
module_biclust <- function(ceRExp, 
                           mRExp = NULL, 
                           BCmethod = "fabia", 
                           num.modules = 10,
                           num.ModuleceRs = 2, 
                           num.ModulemRs = 2) {

    if(is.null(mRExp)){
    ExpData <- assay(ceRExp)  
    } else {
    ExpData <- cbind(assay(ceRExp), assay(mRExp))
    }

    if (BCmethod == "fabia") {
        BCres <- fabia(t(ExpData), p = num.modules)
    } else if (BCmethod == "fabiap") {
        BCres <- fabiap(t(ExpData), p = num.modules)
    } else if (BCmethod == "fabias") {
        BCres <- fabias(t(ExpData), p = num.modules)
    } else if (BCmethod == "mfsc") {
        BCres <- mfsc(t(ExpData), p = num.modules)
    } else if (BCmethod == "nmfdiv") {
        BCres <- nmfdiv(t(ExpData), p = num.modules)
    } else if (BCmethod == "nmfeu") {
        BCres <- nmfeu(t(ExpData), p = num.modules)
    } else if (BCmethod == "nmfsc") {
        BCres <- nmfsc(t(ExpData), p = num.modules)
    } else if (BCmethod == "FLOC") {
        ExpData <- t(ExpData)
        BCres <- FLOC(ExpData, k = num.modules)
    } else if (BCmethod == "isa") {
        BCres <- isa(ExpData)
        BCres <- isa.biclust(BCres)
    } 

    # Extract genes of each bicluster
    if (BCmethod ==  "isa") {
        BCresnum <- biclusternumber(BCres)
        Modulegenes <- lapply(seq_along(BCresnum), function(i) colnames(ExpData)
            [BCresnum[[i]]$Cols])
    }

    if (BCmethod == "fabia" | BCmethod == "fabiap" | BCmethod == "fabias" |
        BCmethod == "mfsc" | BCmethod == "nmfdiv" | BCmethod == "nmfeu" |
        BCmethod == "nmfsc") {
        Modulegenes <- lapply(seq_len(num.modules), function(i) extractBic(BCres)$bic[i,
            ]$bixn)
    }

    if (BCmethod == "FLOC") {
        Modulegenes <- lapply(seq_len(num.modules), function(i) rownames(bicluster(BCres,
            i, graph = FALSE)))
    }

    CandidateModulegenes <- CandModgenes(ceRExp, mRExp = mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}


#' Generation of positively correlated binary matrix between 
#' ceRNAs, or ceRNAs and mRNAs
#'
#' @title cor_binary
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param cor.method The method of calculating correlation selected, 
#' including 'pearson' (default), 'kendall', 'spearman'.
#' @param pos.p.value.cutoff The significant p-value cutoff of 
#' positive correlation.
#' @import SummarizedExperiment
#' @importFrom WGCNA cor
#' @importFrom WGCNA corPvalueFisher
#' @export
#' @return A binary matrix.
#'
#' @examples
#' data(BRCASampleData)
#' cor_binary_matrix <- cor_binary(ceRExp, mRExp)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Langfelder P, Horvath S. WGCNA: an R package for 
#' weighted correlation network analysis. BMC Bioinformatics. 
#' 2008, 9:559.
cor_binary <- function(ceRExp, 
                       mRExp = NULL, 
                       cor.method = "pearson", 
                       pos.p.value.cutoff = 0.01) {
    
    if(is.null(mRExp)){
    cor.r <- cor(assay(ceRExp), method = cor.method)
    } else {
    cor.r <- cor(assay(ceRExp), assay(mRExp), method = cor.method)
    }
    cor.pvalue <- corPvalueFisher(cor.r, nSamples = dim(ceRExp)[1])

    index1 <- which(cor.r > 0)
    index2 <- c(which(cor.r <= 0), which(cor.r %in% NA))
    index3 <- which(cor.pvalue < pos.p.value.cutoff)
    index4 <- c(which(cor.pvalue >= pos.p.value.cutoff), which(cor.pvalue %in%
        NA))

    cor.r[index1] <- 1
    cor.r[index2] <- 0
    cor.pvalue[index3] <- 1
    cor.pvalue[index4] <- 0

    cor.binary <- cor.r * cor.pvalue

    return(cor.binary)
}

## Constructing miRNA-target binary matrix using putative miRNA-target interactions
Bindingmatrix <- function(miRExp = NULL, 
                          ceRExp, 
                          mRExp = NULL, 
                          miRTarget){
  
  miRTarget <- assay(miRTarget)
  mir <- as.character(miRTarget[, 1])
  gene <- as.character(miRTarget[, 2])
  
  if(!is.null(miRExp)){
    miRNA <- colnames(miRExp)
  } else {
    miRNA <- unique(mir)
  }
  
  if(is.null(mRExp)){
    target <- colnames(ceRExp)
  } else {
    target <- c(colnames(ceRExp), colnames(mRExp))
  }
  
  rep <- replicate(length(miRNA), mir)
  edge <- matrix(FALSE, length(miRNA) + length(target), length(miRNA) + length(target))
  
  for (i in 1:length(mir)){    
    if (length(which(rep[i, ] == miRNA)>0)){
      
      match1 <- which(rep[i, ] == miRNA, arr.ind = TRUE)      
      rep2 <- replicate(length(target), gene[i])
      match2 <- which(rep2 == target, arr.ind = TRUE)
      edge[match1, match2 + length(miRNA)] <- TRUE
      
    }
  }
  
  edge <- edge + t(edge)
  edge <- edge!=0
  edge <- edge[seq(target) + length(miRNA), seq(miRNA)] * 1
  colnames(edge) <- miRNA
  rownames(edge) <- target
  
  return(edge)
}

## Identify miRNA sponge modules using sensitivity canonical correlation (SCC) method
miRSM_SCC <- function(miRExp = NULL, 
                      ceRExp, 
                      mRExp = NULL, 
                      miRTarget, 
                      CandidateModulegenes,  
                      typex = "standard", 
                      typez = "standard", 
                      nperms = 100, 
                      num_shared_miRNAs = 3,
                      pvalue.cutoff = 0.05, 
                      CC.cutoff = 0.8, 
                      SCC.cutoff = 0.1) {
  
  if(!is.null(miRExp)){
    miRNames <- colnames(miRExp)
  }
  ceRNames <- colnames(ceRExp)
  if(!is.null(mRExp)){
    mRNames <- colnames(mRExp)
  }
  CandidateModulegenes <- geneIds(CandidateModulegenes)  
  
  Res <- c()
  miRTarget <- assay(miRTarget)
  if(is.null(mRExp) & length(CandidateModulegenes) < 2){
    index <- NULL
    
  } else if (!is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% ceRNames)), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNAs:ceRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Canonical correlation between a group of ceRNAs and another group of ceRNAs
        perm.out_ceR1_ceR2 <- CCA.permute(assay(ceRExp)[, which(ceRNames %in%
                                                                  CandidateModulegenes[[comb_index[i, 1]]])], assay(ceRExp)[, which(ceRNames %in%
                                                                                                                                      CandidateModulegenes[[comb_index[i, 2]]])], typex = typex, typez = typez,
                                          nperms = nperms)
        out_ceR1_ceR2 <- CCA(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                             assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])],
                             typex = typex, typez = typez, penaltyx = perm.out_ceR1_ceR2$bestpenaltyx,
                             penaltyz = perm.out_ceR1_ceR2$bestpenaltyz, v = perm.out_ceR1_ceR2$v.init)
        M6 <- out_ceR1_ceR2$cor
        
        # Canonical correlation between a group of miRNAs and a group of ceRNAs
        perm.out_miR_ceR1 <- CCA.permute(assay(miRExp)[, which(miRNames %in% tmp3)],
                                         assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                                         typex = typex, typez = typez, nperms = nperms)
        out_miR_ceR1 <- CCA(assay(miRExp)[, which(miRNames %in% tmp3)],
                            assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                            typex = typex, typez = typez, penaltyx = perm.out_miR_ceR1$bestpenaltyx,
                            penaltyz = perm.out_miR_ceR1$bestpenaltyz, v = perm.out_miR_ceR1$v.init)
        M7 <- out_miR_ceR1$cor
        
        # Canonical correlation between a group of miRNAs and another group of
        # ceRNAs
        perm.out_miR_ceR2 <- CCA.permute(assay(miRExp)[, which(miRNames %in% tmp3)], 
                                         assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])],
                                         typex = typex, typez = typez, nperms = nperms)
        out_miR_ceR2 <- CCA(assay(miRExp)[, which(miRNames %in% tmp3)],
                            assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])],
                            typex = typex, typez = typez, penaltyx = perm.out_miR_ceR2$bestpenaltyx,
                            penaltyz = perm.out_miR_ceR2$bestpenaltyz, v = perm.out_miR_ceR2$v.init)
        M8 <- out_miR_ceR2$cor
        
        # Calculate partial canonical correlation between a group of ceRNAs
        # and another group of ceRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity canonical correlation between a group of
        # ceRNAs and another group of ceRNAs on condition a group of miRNAs
        M10 <- M6 - M9
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
        M9 <- NA
        M10 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Canonical correlation of ceRNA1:ceRNA2", "Canonical correlation of miRNAs:ceRNA1",
                       "Canonical correlation of miRNAs:ceRNA2", "Partial canonical correlation of ceRNA1:ceRNA2",
                       "Sensitivity canonical correlation of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Canonical correlation of ceRNA1:ceRNA2"] > CC.cutoff &
                     Res[, "Sensitivity canonical correlation of ceRNA1:ceRNA2"] > SCC.cutoff)
    
  } else if (is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% ceRNames), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNAs:ceRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Canonical correlation between a group of ceRNAs and another group of ceRNAs
        perm.out_ceR1_ceR2 <- CCA.permute(assay(ceRExp)[, which(ceRNames %in%
                                                                  CandidateModulegenes[[comb_index[i, 1]]])], assay(ceRExp)[, which(ceRNames %in%
                                                                                                                                      CandidateModulegenes[[comb_index[i, 2]]])], typex = typex, typez = typez,
                                          nperms = nperms)
        out_ceR1_ceR2 <- CCA(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                             assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])],
                             typex = typex, typez = typez, penaltyx = perm.out_ceR1_ceR2$bestpenaltyx,
                             penaltyz = perm.out_ceR1_ceR2$bestpenaltyz, v = perm.out_ceR1_ceR2$v.init)
        M6 <- out_ceR1_ceR2$cor        
        
      } else {
        M6 <- NA        
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Canonical correlation of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Canonical correlation of ceRNA1:ceRNA2"] > CC.cutoff)
    
  } else if (is.null(miRExp) & !is.null(mRExp)){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% c(ceRNames, mRNames)), ]
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Canonical correlation between a group of ceRNAs and a group of mRNAs
        perm.out_ceR_mR <- CCA.permute(assay(ceRExp)[, which(ceRNames %in%
                                                               CandidateModulegenes[[i]])], assay(mRExp)[, which(mRNames %in%
                                                                                                                   CandidateModulegenes[[i]])], typex = typex, typez = typez,
                                       nperms = nperms)
        out_ceR_mR <- CCA(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                          assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])],
                          typex = typex, typez = typez, penaltyx = perm.out_ceR_mR$bestpenaltyx,
                          penaltyz = perm.out_ceR_mR$bestpenaltyz, v = perm.out_ceR_mR$v.init)
        M6 <- out_ceR_mR$cor        
        
      } else {
        M6 <- NA        
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Canonical correlation of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Canonical correlation of ceRNAs:mRNAs"] > CC.cutoff)
    
  } else {
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% c(ceRNames, mRNames))), ]  
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Canonical correlation between a group of ceRNAs and a group of mRNAs
        perm.out_ceR_mR <- CCA.permute(assay(ceRExp)[, which(ceRNames %in%
                                                               CandidateModulegenes[[i]])], assay(mRExp)[, which(mRNames %in%
                                                                                                                   CandidateModulegenes[[i]])], typex = typex, typez = typez,
                                       nperms = nperms)
        out_ceR_mR <- CCA(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                          assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])],
                          typex = typex, typez = typez, penaltyx = perm.out_ceR_mR$bestpenaltyx,
                          penaltyz = perm.out_ceR_mR$bestpenaltyz, v = perm.out_ceR_mR$v.init)
        M6 <- out_ceR_mR$cor
        
        # Canonical correlation between a group of miRNAs and a group of mRNAs
        perm.out_miR_mR <- CCA.permute(assay(miRExp)[, which(miRNames %in% tmp3)],
                                       assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])],
                                       typex = typex, typez = typez, nperms = nperms)
        out_miR_mR <- CCA(assay(miRExp)[, which(miRNames %in% tmp3)],
                          assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])],
                          typex = typex, typez = typez, penaltyx = perm.out_miR_mR$bestpenaltyx,
                          penaltyz = perm.out_miR_mR$bestpenaltyz, v = perm.out_miR_mR$v.init)
        M7 <- out_miR_mR$cor
        
        # Canonical correlation between a group of miRNAs and a group of
        # ceRNAs
        perm.out_miR_ceR <- CCA.permute(assay(miRExp)[, which(miRNames %in% tmp3)], 
                                        assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                                        typex = typex, typez = typez, nperms = nperms)
        out_miR_ceR <- CCA(assay(miRExp)[, which(miRNames %in% tmp3)],
                           assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                           typex = typex, typez = typez, penaltyx = perm.out_miR_ceR$bestpenaltyx,
                           penaltyz = perm.out_miR_ceR$bestpenaltyz, v = perm.out_miR_ceR$v.init)
        M8 <- out_miR_ceR$cor
        
        # Calculate partial canonical correlation between a group of ceRNAs
        # and a group of mRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity canonical correlation between a group of
        # ceRNAs and a group of mRNAs on condition a group of miRNAs
        M10 <- M6 - M9
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
        M9 <- NA
        M10 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Canonical correlation of ceRNAs:mRNAs", "Canonical correlation of miRNAs:mRNAs",
                       "Canonical correlation of miRNAs:ceRNAs", "Partial canonical correlation of ceRNAs:mRNAs",
                       "Sensitivity canonical correlation of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Canonical correlation of ceRNAs:mRNAs"] > CC.cutoff &
                     Res[, "Sensitivity canonical correlation of ceRNAs:mRNAs"] > SCC.cutoff)
  }
  
  if (length(index) == 0) {
    Result <- "No miRNA sponge modules identified!"
  } else {
    if(is.null(mRExp)){
      miRSM_genes <- lapply(index, function(i) list(CandidateModulegenes[[comb_index[i, 1]]], CandidateModulegenes[[comb_index[i, 2]]]))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA1", "ceRNA2")
      }
    } else {
      miRSM_genes <- lapply(index, function(i) list(intersect(CandidateModulegenes[[i]], ceRNames), intersect(CandidateModulegenes[[i]], mRNames)))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA", "mRNA")
      }
    }
    names(miRSM_genes) <- paste("miRSM", seq_along(index), sep=" ")
    Res <- Res[index, ]
    if (length(index) > 1) {
      rownames(Res) <- paste("miRSM", seq_along(index), sep = " ")
    }
    Result <- list(Res, miRSM_genes)
    names(Result) <- c("Group competition of miRNA sponge modules", "miRNA sponge modules")
  }
  return(Result)
}

## Identify miRNA sponge modules using sensitivity distance correlation (SDC) method 
miRSM_SDC <- function(miRExp = NULL, 
                      ceRExp, 
                      mRExp = NULL, 
                      miRTarget, 
                      CandidateModulegenes,  
                      num_shared_miRNAs = 3, 
                      pvalue.cutoff = 0.05, 
                      DC.cutoff = 0.8, 
                      SDC.cutoff = 0.1) {    
  
  if(!is.null(miRExp)){
    miRNames <- colnames(miRExp)
  }
  ceRNames <- colnames(ceRExp)
  if(!is.null(mRExp)){
    mRNames <- colnames(mRExp)
  }
  CandidateModulegenes <- geneIds(CandidateModulegenes)
  
  Res <- c()
  miRTarget <- assay(miRTarget)
  if(is.null(mRExp) & length(CandidateModulegenes) < 2){
    index <- NULL
    
  } else if (!is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% ceRNames)), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNAs:ceRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Calculate distance correlation between a group of ceRNAs and another group of ceRNAs      
        M6 <- dcor(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                   assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # Calculate partial distance correlation between a group of ceRNAs
        # and another group of ceRNAs on condition a group of miRNAs
        M7 <- abs(pdcor(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                        assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])],
                        assay(miRExp)[, which(miRNames %in% tmp3)]))
        
        # Calculate sensitivity distance correlation between a group of
        # ceRNAs and another group of ceRNAs on condition a group of miRNAs
        M8 <- M6 - M7
        
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA      
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Distance correlation of ceRNA1:ceRNA2", "Partial distance correlation of ceRNA1:ceRNA2",
                       "Sensitivity distance correlation of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Distance correlation of ceRNA1:ceRNA2"] > DC.cutoff &
                     Res[, "Sensitivity distance correlation of ceRNA1:ceRNA2"] > SDC.cutoff)
    
  } else if (is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% ceRNames), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNAs:ceRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Calculate distance correlation between a group of ceRNAs and another group of ceRNAs      
        M6 <- dcor(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                   assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
      } else {
        M6 <- NA         
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Distance correlation of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Distance correlation of ceRNA1:ceRNA2"] > DC.cutoff)
    
  } else if (is.null(miRExp) & !is.null(mRExp)){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% c(ceRNames, mRNames)), ]
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {        
        
        # Calculate distance correlation between a group of ceRNAs
        # and a group of mRNAs
        M6 <- dcor(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                   assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])        
        
      } else {
        M6 <- NA        
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Distance correlation of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Distance correlation of ceRNAs:mRNAs"] > DC.cutoff)
    
  } else {
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% c(ceRNames, mRNames))), ]  
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {        
        
        # Calculate distance correlation between a group of ceRNAs
        # and a group of mRNAs
        M6 <- dcor(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                   assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # Calculate partial distance correlation between a group of ceRNAs
        # and a group of mRNAs on condition a group of miRNAs
        M7 <- abs(pdcor(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                        assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])],
                        assay(miRExp)[, which(miRNames %in% tmp3)]))
        
        # Calculate sensitivity distance correlation between a group of
        # ceRNAs and a group of mRNAs on condition a group of miRNAs
        M8 <- M6 - M7
        
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Distance correlation of ceRNAs:mRNAs", "Partial distance correlation of ceRNAs:mRNAs",
                       "Sensitivity distance correlation of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Distance correlation of ceRNAs:mRNAs"] > DC.cutoff &
                     Res[, "Sensitivity distance correlation of ceRNAs:mRNAs"] > SDC.cutoff)
  }
  
  if (length(index) == 0) {
    Result <- "No miRNA sponge modules identified!"
  } else {
    if(is.null(mRExp)){
      miRSM_genes <- lapply(index, function(i) list(CandidateModulegenes[[comb_index[i, 1]]], CandidateModulegenes[[comb_index[i, 2]]]))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA1", "ceRNA2")
      }
    } else {
      miRSM_genes <- lapply(index, function(i) list(intersect(CandidateModulegenes[[i]], ceRNames), intersect(CandidateModulegenes[[i]], mRNames)))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA", "mRNA")
      }
    }
    names(miRSM_genes) <- paste("miRSM", seq_along(index), sep=" ")
    Res <- Res[index, ]
    if (length(index) > 1) {
      rownames(Res) <- paste("miRSM", seq_along(index), sep = " ")
    }
    Result <- list(Res, miRSM_genes)
    names(Result) <- c("Group competition of miRNA sponge modules", "miRNA sponge modules")
  }
  return(Result)
}

## Identify miRNA sponge modules using sensitivity RV coefficient (SRVC) method
miRSM_SRVC <- function(miRExp = NULL, 
                       ceRExp, 
                       mRExp = NULL, 
                       miRTarget, 
                       CandidateModulegenes,  
                       num_shared_miRNAs = 3, 
                       pvalue.cutoff = 0.05, 
                       RVC.cutoff = 0.8, 
                       SRVC.cutoff = 0.1, 
                       RV_method = "RV") {    
  
  if(!is.null(miRExp)){
    miRNames <- colnames(miRExp)
  }
  ceRNames <- colnames(ceRExp)
  if(!is.null(mRExp)){
    mRNames <- colnames(mRExp)
  }
  CandidateModulegenes <- geneIds(CandidateModulegenes)
  
  Res <- c()
  miRTarget <- assay(miRTarget)
  if(is.null(mRExp) & length(CandidateModulegenes) < 2){
    index <- NULL
    
  } else if (!is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% ceRNames)), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNA1:ceRNA2
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs & RV_method == "RV") {
        
        # RV coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- RV(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                 assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNA1        
        M7 <- RV(assay(miRExp)[, which(miRNames %in% tmp3)],
                 assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNA2         
        M8 <- RV(assay(miRExp)[, which(miRNames %in% tmp3)],
                 assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # Calculate partial RV coefficient between a group of ceRNAs
        # and another group of ceRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity RV coefficient between a group of
        # ceRNAs and another group of ceRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RV2") {
        
        # RV coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- RV2(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNA1        
        M7 <- RV2(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNA2         
        M8 <- RV2(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # Calculate partial RV coefficient between a group of ceRNAs
        # and another group of ceRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity RV coefficient between a group of
        # ceRNAs and another group of ceRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RVadjMaye") {
        
        # RV coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- RVadjMaye(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                        assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNA1        
        M7 <- RVadjMaye(assay(miRExp)[, which(miRNames %in% tmp3)],
                        assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNA2         
        M8 <- RVadjMaye(assay(miRExp)[, which(miRNames %in% tmp3)],
                        assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # Calculate partial RV coefficient between a group of ceRNAs
        # and another group of ceRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity RV coefficient between a group of
        # ceRNAs and another group of ceRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RVadjGhaziri") {
        
        # RV coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- RVadjGhaziri(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                           assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNA1        
        M7 <- RVadjGhaziri(assay(miRExp)[, which(miRNames %in% tmp3)],
                           assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNA2         
        M8 <- RVadjGhaziri(assay(miRExp)[, which(miRNames %in% tmp3)],
                           assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # Calculate partial RV coefficient between a group of ceRNAs
        # and another group of ceRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity RV coefficient between a group of
        # ceRNAs and another group of ceRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
        M9 <- NA
        M10 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
      Res <- rbind(Res, tmp)
    }
    
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "RV coefficient of ceRNA1:ceRNA2", "RV coefficient of miRNAs:ceRNA1",
                       "RV coefficient of miRNAs:ceRNA2", "Partial RV coefficient of ceRNA1:ceRNA2",
                       "Sensitivity RV coefficient of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "RV coefficient of ceRNA1:ceRNA2"] > RVC.cutoff &
                     Res[, "Sensitivity RV coefficient of ceRNA1:ceRNA2"] > SRVC.cutoff)
    
  } else if (is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% ceRNames), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNA1:ceRNA2
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs & RV_method == "RV") {
        
        # RV coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- RV(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                 assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RV2") {
        
        # RV coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- RV2(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])        
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RVadjMaye") {
        
        # RV coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- RVadjMaye(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                        assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])        
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RVadjGhaziri") {
        
        # RV coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- RVadjGhaziri(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                           assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])        
        
      } else {
        M6 <- NA       
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
    }
    
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "RV coefficient of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "RV coefficient of ceRNA1:ceRNA2"] > RVC.cutoff)
    
  } else if (is.null(miRExp) & !is.null(mRExp)){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% c(ceRNames, mRNames)), ]
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs & RV_method == "RV") {
        
        # RV coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- RV(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                 assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RV2") {
        
        # RV coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- RV2(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                  assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])        
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RVadjMaye") {
        
        # RV coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- RVadjMaye(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                        assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])        
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RVadjGhaziri") {
        
        # RV coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- RVadjGhaziri(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                           assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])        
        
      } else {
        M6 <- NA       
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "RV coefficient of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "RV coefficient of ceRNAs:mRNAs"] > RVC.cutoff)
    
  } else {  
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% c(ceRNames, mRNames))), ]
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs & RV_method == "RV") {
        
        # RV coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- RV(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                 assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # RV coefficient between a group of miRNAs and a group of mRNAs        
        M7 <- RV(assay(miRExp)[, which(miRNames %in% tmp3)],
                 assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNAs         
        M8 <- RV(assay(miRExp)[, which(miRNames %in% tmp3)],
                 assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])
        
        # Calculate partial RV coefficient between a group of ceRNAs
        # and a group of mRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity RV coefficient between a group of
        # ceRNAs and a group of mRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RV2") {
        
        # RV coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- RV2(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                  assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # RV coefficient between a group of miRNAs and a group of mRNAs        
        M7 <- RV2(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNAs         
        M8 <- RV2(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])
        
        # Calculate partial RV coefficient between a group of ceRNAs
        # and a group of mRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity RV coefficient between a group of
        # ceRNAs and a group of mRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RVadjMaye") {
        
        # RV coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- RVadjMaye(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                        assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # RV coefficient between a group of miRNAs and a group of mRNAs        
        M7 <- RVadjMaye(assay(miRExp)[, which(miRNames %in% tmp3)],
                        assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNAs         
        M8 <- RVadjMaye(assay(miRExp)[, which(miRNames %in% tmp3)],
                        assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])
        
        # Calculate partial RV coefficient between a group of ceRNAs
        # and a group of mRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity RV coefficient between a group of
        # ceRNAs and a group of mRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else if (M3 >= num_shared_miRNAs & RV_method == "RVadjGhaziri") {
        
        # RV coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- RVadjGhaziri(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                           assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # RV coefficient between a group of miRNAs and a group of mRNAs        
        M7 <- RVadjGhaziri(assay(miRExp)[, which(miRNames %in% tmp3)],
                           assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # RV coefficient between a group of miRNAs and a group of ceRNAs         
        M8 <- RVadjGhaziri(assay(miRExp)[, which(miRNames %in% tmp3)],
                           assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])
        
        # Calculate partial RV coefficient between a group of ceRNAs
        # and a group of mRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity RV coefficient between a group of
        # ceRNAs and a group of mRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
        M9 <- NA
        M10 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "RV coefficient of ceRNAs:mRNAs", "RV coefficient of miRNAs:mRNAs",
                       "RV coefficient of miRNAs:ceRNAs", "Partial RV coefficient of ceRNAs:mRNAs",
                       "Sensitivity RV coefficient of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "RV coefficient of ceRNAs:mRNAs"] > RVC.cutoff &
                     Res[, "Sensitivity RV coefficient of ceRNAs:mRNAs"] > SRVC.cutoff)
  }
  
  if (length(index) == 0) {
    Result <- "No miRNA sponge modules identified!"
  } else {
    if(is.null(mRExp)){
      miRSM_genes <- lapply(index, function(i) list(CandidateModulegenes[[comb_index[i, 1]]], CandidateModulegenes[[comb_index[i, 2]]]))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA1", "ceRNA2")
      }
    } else {
      miRSM_genes <- lapply(index, function(i) list(intersect(CandidateModulegenes[[i]], ceRNames), intersect(CandidateModulegenes[[i]], mRNames)))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA", "mRNA")
      }
    }
    names(miRSM_genes) <- paste("miRSM", seq_along(index), sep=" ")
    Res <- Res[index, ]
    if (length(index) > 1) {
      rownames(Res) <- paste("miRSM", seq_along(index), sep = " ")
    }
    Result <- list(Res, miRSM_genes)
    names(Result) <- c("Group competition of miRNA sponge modules", "miRNA sponge modules")
  }
  return(Result)
}

## Identify miRNA sponge modules using sponge module identification (SM) method in 
## the reference "Zhang J, Le TD, Liu L, Li J. Identifying miRNA sponge modules using biclustering and regulatory scores. 
## BMC Bioinformatics. 2017 Mar 14;18(Suppl 3):44. doi: 10.1186/s12859-017-1467-5.".
## Note that miRNA-mRNA correlation matrix using Pearson method and miRNA-mRNA context++ score matrix using putative miRNA-target binding information
## are converted into two miRNA-mRNA binary matrices.
## a and b are the contributions of expression data and miRNA-target binding information, respectively.
miRSM_SM <- function(miRExp = NULL, 
                     ceRExp, 
                     mRExp = NULL, 
                     miRTarget, 
                     a = 0.5, 
                     b = 0.5, 
                     BCmethod = "BCPlaid", 
                     pvalue.cutoff = 0.05) {
  
  Binding_edge <- Bindingmatrix(miRExp = miRExp, ceRExp, mRExp = mRExp, miRTarget)
  
  if(is.null(miRExp)){
    colScore <- Binding_edge
    
  } else {
    if(is.null(mRExp)){
      Cor.Pvalue <- WGCNA::corAndPvalue(assay(ceRExp), assay(miRExp))$p
    } else {
      Cor.Pvalue <- WGCNA::corAndPvalue(cbind(assay(ceRExp), assay(mRExp)), assay(miRExp))$p
    }
    index1 <- which(Cor.Pvalue < pvalue.cutoff)
    index2 <- c(which(Cor.Pvalue >= pvalue.cutoff), 
                which(Cor.Pvalue %in% NA))
    Cor.Pvalue[index1] <- 1
    Cor.Pvalue[index2] <- 0
    Cor_edge <- Cor.Pvalue
    
    colScore <- a*Cor_edge + b*Binding_edge
  }
  
  if (BCmethod=="BCBimax"){
    colScore <- binarize(colScore)
    BCres <- biclust(colScore, BCBimax())
  } else if (BCmethod=="BCCC") {
    BCres <- biclust(colScore, BCCC())
  } else if (BCmethod=="BCPlaid") {
    BCres <- biclust(colScore, BCPlaid())
  } else if (BCmethod=="BCQuest") {
    BCres <- biclust(colScore, BCQuest())
  } else if (BCmethod=="BCSpectral") {
    BCres <- biclust(colScore, BCSpectral())
  } else if (BCmethod=="BCXmotifs") {
    colScore <- discretize(colScore)
    BCres <- biclust(colScore, BCXmotifs())
  } 
  
  BCresnum <- biclusternumber(BCres)
  Modules <- lapply(seq_along(BCresnum), function(i) rownames(Binding_edge)
                    [BCresnum[[i]]$Rows])
  
  return(Modules)
}

## Identify miRNA sponge modules using sensitivity similarity index (SSI) method
miRSM_SSI <- function(miRExp = NULL, 
                      ceRExp, 
                      mRExp = NULL, 
                      miRTarget, 
                      CandidateModulegenes,  
                      num_shared_miRNAs = 3, 
                      pvalue.cutoff = 0.05, 
                      SI.cutoff = 0.8, 
                      SSI.cutoff = 0.1) {    
  
  if(!is.null(miRExp)){
    miRNames <- colnames(miRExp)
  }
  ceRNames <- colnames(ceRExp)
  if(!is.null(mRExp)){
    mRNames <- colnames(mRExp)
  }
  CandidateModulegenes <- geneIds(CandidateModulegenes)
  
  Res <- c()
  miRTarget <- assay(miRTarget)
  if(is.null(mRExp) & length(CandidateModulegenes) < 2){
    index <- NULL
    
  } else if (!is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% ceRNames)), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNA1:ceRNA2
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Similarity index between a group of ceRNAs and another group of ceRNAs       
        M6 <- SMI(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])], 1, 1)
        
        # Similarity index between a group of miRNAs and a group of ceRNA1        
        M7 <- SMI(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])], 1, 1)
        
        # Similarity index between a group of miRNAs and a group of ceRNA2         
        M8 <- SMI(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])], 1, 1)
        
        # Calculate partial similarity index between a group of ceRNAs
        # and another group of ceRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity similarity index between a group of
        # ceRNAs and another group of ceRNAs on condition a group of miRNAs
        M10 <- M6 - M9      
        
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
        M9 <- NA
        M10 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
      Res <- rbind(Res, tmp)
    }
    
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Similarity index of ceRNA1:ceRNA2", "Similarity index of miRNAs:ceRNA1",
                       "Similarity index of miRNAs:ceRNA2", "Partial similarity index of ceRNA1:ceRNA2",
                       "Sensitivity similarity index of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Similarity index of ceRNA1:ceRNA2"] > SI.cutoff &
                     Res[, "Sensitivity similarity index of ceRNA1:ceRNA2"] > SSI.cutoff)
    
  } else if (is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% ceRNames), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNA1:ceRNA2
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Similarity index between a group of ceRNAs and another group of ceRNAs       
        M6 <- SMI(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])], 1, 1)        
        
      } else {
        M6 <- NA       
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
    }
    
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Similarity index of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Similarity index of ceRNA1:ceRNA2"] > SI.cutoff)
    
  } else if (is.null(miRExp) & !is.null(mRExp)){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% c(ceRNames, mRNames)), ]
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Similarity index between a group of ceRNAs and a group of mRNAs       
        M6 <- SMI(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                  assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])], 1, 1)        
        
      } else {
        M6 <- NA       
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Similarity index of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Similarity index of ceRNAs:mRNAs"] > SI.cutoff)
    
  } else {  
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% c(ceRNames, mRNames))), ]
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Similarity index between a group of ceRNAs and a group of mRNAs       
        M6 <- SMI(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                  assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])], 1, 1)
        
        # Similarity index between a group of miRNAs and a group of mRNAs        
        M7 <- SMI(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])], 1, 1)
        
        # Similarity index between a group of miRNAs and a group of ceRNAs         
        M8 <- SMI(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])], 1, 1)
        
        # Calculate partial similarity index between a group of ceRNAs
        # and a group of mRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity similarity index between a group of
        # ceRNAs and a group of mRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
        M9 <- NA
        M10 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Similarity index of ceRNAs:mRNAs", "Similarity index of miRNAs:mRNAs",
                       "Similarity index of miRNAs:ceRNAs", "Partial similarity index of ceRNAs:mRNAs",
                       "Sensitivity similarity index of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Similarity index of ceRNAs:mRNAs"] > SI.cutoff &
                     Res[, "Sensitivity similarity index of ceRNAs:mRNAs"] > SSI.cutoff)
  }
  
  if (length(index) == 0) {
    Result <- "No miRNA sponge modules identified!"
  } else {
    if(is.null(mRExp)){
      miRSM_genes <- lapply(index, function(i) list(CandidateModulegenes[[comb_index[i, 1]]], CandidateModulegenes[[comb_index[i, 2]]]))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA1", "ceRNA2")
      }
    } else {
      miRSM_genes <- lapply(index, function(i) list(intersect(CandidateModulegenes[[i]], ceRNames), intersect(CandidateModulegenes[[i]], mRNames)))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA", "mRNA")
      }
    }
    names(miRSM_genes) <- paste("miRSM", seq_along(index), sep=" ")
    Res <- Res[index, ]
    if (length(index) > 1) {
      rownames(Res) <- paste("miRSM", seq_along(index), sep = " ")
    }
    Result <- list(Res, miRSM_genes)
    names(Result) <- c("Group competition of miRNA sponge modules", "miRNA sponge modules")
  }
  return(Result)
}

## Identify miRNA sponge modules using sensitivity generalized coefficient of determination (SGCD) method
miRSM_SGCD <- function(miRExp = NULL, 
                       ceRExp, 
                       mRExp = NULL, 
                       miRTarget, 
                       CandidateModulegenes,  
                       num_shared_miRNAs = 3, 
                       pvalue.cutoff = 0.05, 
                       GCD.cutoff = 0.8, 
                       SGCD.cutoff = 0.1) {    
  
  if(!is.null(miRExp)){
    miRNames <- colnames(miRExp)
  }
  ceRNames <- colnames(ceRExp)
  if(!is.null(mRExp)){
    mRNames <- colnames(mRExp)
  }
  CandidateModulegenes <- geneIds(CandidateModulegenes)
  
  Res <- c()
  miRTarget <- assay(miRTarget)
  if(is.null(mRExp) & length(CandidateModulegenes) < 2){
    index <- NULL
    
  } else if (!is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% ceRNames)), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNA1:ceRNA2
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Generalized coefficient of determination between a group of ceRNAs and another group of ceRNAs       
        M6 <- GCD(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # Generalized coefficient of determination between a group of miRNAs and a group of ceRNA1        
        M7 <- GCD(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])])
        
        # Generalized coefficient of determination between a group of miRNAs and a group of ceRNA2         
        M8 <- GCD(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
        # Calculate partial generalized coefficient of determination between a group of ceRNAs
        # and another group of ceRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity generalized coefficient of determination between a group of
        # ceRNAs and another group of ceRNAs on condition a group of miRNAs
        M10 <- M6 - M9      
        
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
        M9 <- NA
        M10 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
      Res <- rbind(Res, tmp)
    }
    
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Generalized coefficient of determination of ceRNA1:ceRNA2", "Generalized coefficient of determination of miRNAs:ceRNA1",
                       "Generalized coefficient of determination of miRNAs:ceRNA2", "Partial generalized coefficient of determination of ceRNA1:ceRNA2",
                       "Sensitivity generalized coefficient of determination of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Generalized coefficient of determination of ceRNA1:ceRNA2"] > GCD.cutoff &
                     Res[, "Sensitivity generalized coefficient of determination of ceRNA1:ceRNA2"] > SGCD.cutoff)
    
  } else if (is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% ceRNames), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNA1:ceRNA2
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Generalized coefficient of determination between a group of ceRNAs and another group of ceRNAs       
        M6 <- GCD(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])
        
      } else {
        M6 <- NA        
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
    }
    
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Generalized coefficient of determination of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Generalized coefficient of determination of ceRNA1:ceRNA2"] > GCD.cutoff)
    
  } else if (is.null(miRExp) & !is.null(mRExp)){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% c(ceRNames, mRNames)), ]
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Generalized coefficient of determination between a group of ceRNAs and a group of mRNAs       
        M6 <- GCD(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                  assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])        
        
      } else {
        M6 <- NA        
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Generalized coefficient of determination of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Generalized coefficient of determination of ceRNAs:mRNAs"] > GCD.cutoff)
    
  } else {  
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% c(ceRNames, mRNames))), ]
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs) {
        
        # Generalized coefficient of determination between a group of ceRNAs and a group of mRNAs       
        M6 <- GCD(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                  assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # Generalized coefficient of determination between a group of miRNAs and a group of mRNAs        
        M7 <- GCD(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
        
        # Generalized coefficient of determination between a group of miRNAs and a group of ceRNAs         
        M8 <- GCD(assay(miRExp)[, which(miRNames %in% tmp3)],
                  assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])
        
        # Calculate partial generalized coefficient of determination between a group of ceRNAs
        # and a group of mRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity generalized coefficient of determination between a group of
        # ceRNAs and a group of mRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
        M9 <- NA
        M10 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "Generalized coefficient of determination of ceRNAs:mRNAs", "Generalized coefficient of determination of miRNAs:mRNAs",
                       "Generalized coefficient of determination of miRNAs:ceRNAs", "Partial generalized coefficient of determination of ceRNAs:mRNAs",
                       "Sensitivity generalized coefficient of determination of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "Generalized coefficient of determination of ceRNAs:mRNAs"] > GCD.cutoff &
                     Res[, "Sensitivity generalized coefficient of determination of ceRNAs:mRNAs"] > SGCD.cutoff)
  }
  
  if (length(index) == 0) {
    Result <- "No miRNA sponge modules identified!"
  } else {
    if(is.null(mRExp)){
      miRSM_genes <- lapply(index, function(i) list(CandidateModulegenes[[comb_index[i, 1]]], CandidateModulegenes[[comb_index[i, 2]]]))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA1", "ceRNA2")
      }
    } else {
      miRSM_genes <- lapply(index, function(i) list(intersect(CandidateModulegenes[[i]], ceRNames), intersect(CandidateModulegenes[[i]], mRNames)))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA", "mRNA")
      }
    }
    names(miRSM_genes) <- paste("miRSM", seq_along(index), sep=" ")
    Res <- Res[index, ]
    if (length(index) > 1) {
      rownames(Res) <- paste("miRSM", seq_along(index), sep = " ")
    }
    Result <- list(Res, miRSM_genes)
    names(Result) <- c("Group competition of miRNA sponge modules", "miRNA sponge modules")
  }
  return(Result)
}

## Identify miRNA sponge modules using sensitivity Coxhead's or Rozeboom's coefficient (SCRC) method
miRSM_SCRC <- function(miRExp = NULL, 
                       ceRExp, 
                       mRExp = NULL, 
                       miRTarget, 
                       CandidateModulegenes,  
                       num_shared_miRNAs = 3, 
                       pvalue.cutoff = 0.05, 
                       CRC.cutoff = 0.8, 
                       SCRC.cutoff = 0.1, 
                       CRC_method = "Coxhead") {    
  
  if(!is.null(miRExp)){
    miRNames <- colnames(miRExp)
  }
  ceRNames <- colnames(ceRExp)
  if(!is.null(mRExp)){
    mRNames <- colnames(mRExp)
  }
  CandidateModulegenes <- geneIds(CandidateModulegenes)
  
  Res <- c()
  miRTarget <- assay(miRTarget)
  if(is.null(mRExp) & length(CandidateModulegenes) < 2){
    index <- NULL
    
  } else if (!is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% ceRNames)), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNA1:ceRNA2
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs & CRC_method == "Coxhead") {
        
        # Coxhead's coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- Coxhead(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                      assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])[1]
        
        # Coxhead's coefficient between a group of miRNAs and a group of ceRNA1        
        M7 <- Coxhead(assay(miRExp)[, which(miRNames %in% tmp3)],
                      assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])])[1]
        
        # Coxhead's coefficient between a group of miRNAs and a group of ceRNA2         
        M8 <- Coxhead(assay(miRExp)[, which(miRNames %in% tmp3)],
                      assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])[1]
        
        # Calculate partial Coxhead's coefficient between a group of ceRNAs
        # and another group of ceRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity Coxhead's coefficient between a group of
        # ceRNAs and another group of ceRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else if (M3 >= num_shared_miRNAs & CRC_method == "Rozeboom") {
        
        # Rozeboom's coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- Rozeboom(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                       assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])[1]
        
        # Rozeboom's coefficient between a group of miRNAs and a group of ceRNA1        
        M7 <- Rozeboom(assay(miRExp)[, which(miRNames %in% tmp3)],
                       assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])])[1]
        
        # Rozeboom's coefficient between a group of miRNAs and a group of ceRNA2         
        M8 <- Rozeboom(assay(miRExp)[, which(miRNames %in% tmp3)],
                       assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])[1]
        
        # Calculate partial Rozeboom's coefficient between a group of ceRNAs
        # and another group of ceRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity Rozeboom's coefficient between a group of
        # ceRNAs and another group of ceRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
        M9 <- NA
        M10 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
      Res <- rbind(Res, tmp)
    }
    
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "CRC coefficient of ceRNA1:ceRNA2", "CRC coefficient of miRNAs:ceRNA1",
                       "CRC coefficient of miRNAs:ceRNA2", "Partial CRC coefficient of ceRNA1:ceRNA2",
                       "Sensitivity CRC coefficient of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "CRC coefficient of ceRNA1:ceRNA2"] > CRC.cutoff &
                     Res[, "Sensitivity CRC coefficient of ceRNA1:ceRNA2"] > SCRC.cutoff)
    
  } else if (is.null(miRExp) & is.null(mRExp) & length(CandidateModulegenes) > 1){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% ceRNames), ]
    
    comb_index <- t(combn(length(CandidateModulegenes), 2))
    
    for (i in seq(nrow(comb_index))){
      # Calculate significance of miRNAs shared by each ceRNA1:ceRNA2
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 1]]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[comb_index[i, 2]]], ceRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs & CRC_method == "Coxhead") {
        
        # Coxhead's coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- Coxhead(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                      assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])[1]        
        
      } else if (M3 >= num_shared_miRNAs & CRC_method == "Rozeboom") {
        
        # Rozeboom's coefficient between a group of ceRNAs and another group of ceRNAs       
        M6 <- Rozeboom(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 1]]])],
                       assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[comb_index[i, 2]]])])[1]        
        
      } else {
        M6 <- NA       
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
    }
    
    colnames(Res) <- c("#miRNAs regulating ceRNA1", "#miRNAs regulating ceRNA2",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "CRC coefficient of ceRNA1:ceRNA2")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "CRC coefficient of ceRNA1:ceRNA2"] > CRC.cutoff)
    
  } else if (is.null(miRExp) & !is.null(mRExp)){
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% c(ceRNames, mRNames)), ]
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(unique(miRTargetCandidate[, 1]))
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs & CRC_method == "Coxhead") {
        
        # Coxhead's coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- Coxhead(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                      assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])[1]        
        
      } else if (M3 >= num_shared_miRNAs & CRC_method == "Rozeboom") {
        
        # Rozeboom's coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- Rozeboom(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                       assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])[1]        
        
      } else {
        M6 <- NA       
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "CRC coefficient of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "CRC coefficient of ceRNAs:mRNAs"] > CRC.cutoff)
    
  } else {  
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% miRNames), 
                                              which(miRTarget[, 2] %in% c(ceRNames, mRNames))), ]
    
    for (i in seq_along(CandidateModulegenes)) {
      # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
      tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], ceRNames)), 1])
      M1 <- length(tmp1)
      tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                                intersect(CandidateModulegenes[[i]], mRNames)), 1])
      M2 <- length(tmp2)
      tmp3 <- intersect(tmp1, tmp2)
      M3 <- length(tmp3)
      M4 <- length(miRNames)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= num_shared_miRNAs & CRC_method == "Coxhead") {
        
        # Coxhead's coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- Coxhead(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                      assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])[1]
        
        # Coxhead's coefficient between a group of miRNAs and a group of mRNAs        
        M7 <- Coxhead(assay(miRExp)[, which(miRNames %in% tmp3)],
                      assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])[1]
        
        # Coxhead's coefficient between a group of miRNAs and a group of ceRNAs         
        M8 <- Coxhead(assay(miRExp)[, which(miRNames %in% tmp3)],
                      assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])[1]
        
        # Calculate partial Coxhead's coefficient between a group of ceRNAs
        # and a group of mRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity Coxhead's coefficient between a group of
        # ceRNAs and a group of mRNAs on condition a group of miRNAs
        M10 <- M6 - M9
        
      } else if (M3 >= num_shared_miRNAs & CRC_method == "Rozeboom") {
        
        # Rozeboom's coefficient between a group of ceRNAs and a group of mRNAs       
        M6 <- Rozeboom(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                       assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])[1]
        
        # Rozeboom's coefficient between a group of miRNAs and a group of mRNAs        
        M7 <- Rozeboom(assay(miRExp)[, which(miRNames %in% tmp3)],
                       assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])[1]
        
        # Rozeboom's coefficient between a group of miRNAs and a group of ceRNAs         
        M8 <- Rozeboom(assay(miRExp)[, which(miRNames %in% tmp3)],
                       assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])[1]
        
        # Calculate partial Rozeboom's coefficient between a group of ceRNAs
        # and a group of mRNAs on condition a group of miRNAs
        M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
        
        # Calculate sensitivity Rozeboom's coefficient between a group of
        # ceRNAs and a group of mRNAs on condition a group of miRNAs
        M10 <- M6 - M9        
        
      } else {
        M6 <- NA
        M7 <- NA
        M8 <- NA
        M9 <- NA
        M10 <- NA
      }
      
      tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
      Res <- rbind(Res, tmp)
      
    }
    colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                       "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                       "CRC coefficient of ceRNAs:mRNAs", "CRC coefficient of miRNAs:mRNAs",
                       "CRC coefficient of miRNAs:ceRNAs", "Partial CRC coefficient of ceRNAs:mRNAs",
                       "Sensitivity CRC coefficient of ceRNAs:mRNAs")
    index <- which(Res[, "#Shared miRNAs"] >= num_shared_miRNAs &
                     Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                     Res[, "CRC coefficient of ceRNAs:mRNAs"] > CRC.cutoff &
                     Res[, "Sensitivity CRC coefficient of ceRNAs:mRNAs"] > SCRC.cutoff)
  }
  
  if (length(index) == 0) {
    Result <- "No miRNA sponge modules identified!"
  } else {
    if(is.null(mRExp)){
      miRSM_genes <- lapply(index, function(i) list(CandidateModulegenes[[comb_index[i, 1]]], CandidateModulegenes[[comb_index[i, 2]]]))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA1", "ceRNA2")
      }
    } else {
      miRSM_genes <- lapply(index, function(i) list(intersect(CandidateModulegenes[[i]], ceRNames), intersect(CandidateModulegenes[[i]], mRNames)))
      for(i in seq_along(index)){
        names(miRSM_genes[[i]]) <- c("ceRNA", "mRNA")
      }
    }
    names(miRSM_genes) <- paste("miRSM", seq_along(index), sep=" ")
    Res <- Res[index, ]
    if (length(index) > 1) {
      rownames(Res) <- paste("miRSM", seq_along(index), sep = " ")
    }
    Result <- list(Res, miRSM_genes)
    names(Result) <- c("Group competition of miRNA sponge modules", "miRNA sponge modules")
  }
  return(Result)
}

#' Calculating similarity between two list of module groups  
#'
#' @title module_group_sim
#' @param Module.group1 List object, the first list of module group.
#' @param Module.group2 List object, the second list of module group.
#' @param sim.method Methods for calculating similatiry between two modules, select one of three methods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @export
#' @return Similarity between two list of module groups
#'
#' @examples
#' library(GSEABase)
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp) 
#' modulegenes_igraph <- module_igraph (ceRExp, mRExp) 
#' Sim <- module_group_sim(geneIds(modulegenes_WGCNA), geneIds(modulegenes_igraph))
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Simpson E H. Measurement of diversity. Nature, 1949, 163(4148): 688-688.
#' @references Jaccard P. The distribution of the flora in the alpine zone. 1. New phytologist, 1912, 11(2): 37-50.
#' @references Lin D. An information-theoretic definition of similarity. in: Icml. 1998, 98(1998): 296-304.
module_group_sim <- function(Module.group1, 
                             Module.group2, 
                             sim.method = "Simpson"){
  
  if(class(Module.group1) != "list" | class(Module.group2) != "list") {
    stop("Please check your input module group! The input module group should be list object! \n")
  } else if (class(Module.group1[[1]]) == "list" | class(Module.group2[[1]]) == "list"){
    Module.group1 <- lapply(seq(Module.group1), function(i) unique(unlist(Module.group1[[i]])))
    Module.group2 <- lapply(seq(Module.group2), function(i) unique(unlist(Module.group2[[i]])))    
  }     
  
  m <- length(Module.group1)
  n <- length(Module.group2)
  
  Sim <- module_group_sim_matrix(Module.group1, Module.group2, sim.method = sim.method)
  
  if (m < n) {        
    GS <- mean(unlist(lapply(seq(m), function(i) Sim[i, max.col(Sim)[i]])))*m/n
  } else if (m == n) {
    GS <- mean(c(unlist(lapply(seq(m), function(i) Sim[i, max.col(Sim)[i]])), 
                 unlist(lapply(seq(n), function(i) Sim[max.col(t(Sim))[i], i]))))
  } else if (m > n) {
    GS <- mean(unlist(lapply(seq(n), function(i) Sim[max.col(t(Sim))[i], i])))*n/m
  }
  
  return(GS)
}

#' Inferring differential modules between two list of module groups 
#'
#' @title diff_module
#' @param Module.group1 List object, the first list of module group.
#' @param Module.group2 List object, the second list of module group.
#' @param sim.cutoff Similarity cutoff between modules, the interval is [0 1].
#' @param sim.method Methods for calculating similatiry between two modules, select one of three methods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @export
#' @return A list of differential modules
#'
#' @examples
#' library(GSEABase)
#' data(BRCASampleData)
#' modulegenes_WGCNA_all <- module_WGCNA(ceRExp, mRExp)
#' modulegenes_WGCNA_1 <- module_WGCNA(ceRExp[-1, ], mRExp[-1, ])
#' Differential_module <- diff_module(geneIds(modulegenes_WGCNA_all), geneIds(modulegenes_WGCNA_1))
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
diff_module <- function(Module.group1, 
                        Module.group2,
                        sim.cutoff = 0.8, 
                        sim.method = "Simpson"){
  
  if(class(Module.group1) != "list" | class(Module.group2) != "list") {
    stop("Please check your input module group! The input module group should be list object! \n")
  } else if (class(Module.group1[[1]]) == "list" | class(Module.group2[[1]]) == "list"){
    Module.group1.interin <- lapply(seq(Module.group1), function(i) unique(unlist(Module.group1[[i]])))
    Module.group2.interin <- lapply(seq(Module.group2), function(i) unique(unlist(Module.group2[[i]])))
    Sim <- module_group_sim_matrix(Module.group1.interin, Module.group2.interin, sim.method = sim.method)
  } else {
    Sim <- module_group_sim_matrix(Module.group1, Module.group2, sim.method = sim.method)
  }   
  
  row.max.index <- apply(Sim, 1, function(x){which.max(x)})
  col.max.index <- apply(Sim, 2, function(x){which.max(x)})
  row.max <- unlist(lapply(seq(row.max.index), function(i) Sim[i, row.max.index[i]]))
  col.max <- unlist(lapply(seq(col.max.index), function(i) Sim[col.max.index[i], i]))
  
  Diff_Module_row <- lapply(which(row.max < sim.cutoff), function(i) Module.group1[[i]])
  Diff_Module_col <- lapply(which(col.max < sim.cutoff), function(i) Module.group2[[i]])
  Diff_Module <- append(Diff_Module_row, Diff_Module_col)
  
  if (length(Diff_Module) == 0) {
    Diff_Module <- "No differential modules identified!"
  } else {
    names(Diff_Module) <- paste("miRSM", seq(Diff_Module), sep=" ")
  }
  
  return(Diff_Module)
}


#' Identify miRNA sponge modules using sensitivity canonical correlation (SCC), 
#' sensitivity distance correlation (SDC),
#' sensitivity RV coefficient (SRVC), sensitivity similarity index (SSI), 
#' sensitivity generalized coefficient of determination (SGCD), 
#' sensitivity Coxhead's or Rozeboom's coefficient (SCRC), 
#' and sponge module (SM) methods.
#'
#' @title miRSM
#' @param miRExp NULL (default) or a SummarizedExperiment object. miRNA expression data: 
#' rows are samples and columns are miRNAs.
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param miRTarget A SummarizedExperiment object. Putative 
#' miRNA-target binding information.
#' @param CandidateModulegenes List object: a list of candidate 
#' miRNA sponge modules. Only for the SCC, SDC, SRVC, SSI, SGCD and SCRC methods.
#' @param typex The columns of x unordered (type='standard') or 
#' ordered (type='ordered'). Only for the SCC method.
#' @param typez The columns of z unordered (type='standard') or 
#' ordered (type='ordered'). Only for the SCC method.
#' @param nperms The number of permutations. Only for the SCC method.
#' @param method The method selected to identify miRNA sponge 
#' modules, including 'SCC', 'SDC', 'SRVC', 'SM', 'SSI', 'SGCD' and 'SCRC'.
#' @param num_shared_miRNAs The number of common miRNAs shared 
#' by a group of ceRNAs and mRNAs. Only for the SCC, SDC, SRVC, SSI, 
#' SGCD and SCRC methods. 
#' @param pvalue.cutoff The p-value cutoff of significant sharing 
#' of common miRNAs by a group of ceRNAs and mRNAs or significant correlation.
#' @param MC.cutoff The cutoff of matrix correlation (canonical correlation, 
#' distance correlation and RV coefficient). Only for the SCC, SDC, SRVC, 
#' SSI, SGCD and SCRC methods.
#' @param SMC.cutoff The cutoff of sensitivity matrix correlation
#' (sensitivity canonical correlation, sensitivity distance correlation 
#' and sensitivity RV coefficient). Only for the SCC, SDC, SRVC, SSI, 
#' SGCD and SCRC methods when miRExp is not NULL.
#' @param RV_method the method of calculating RV coefficients. Select
#' one of 'RV', 'RV2', 'RVadjMaye' and 'RVadjGhaziri' methods.
#' Only for the SRVC method.
#' @param BCmethod Specification of the biclustering method, 
#' including 'BCBimax', 'BCCC', 'BCPlaid' (default), 'BCQuest', 
#' 'BCSpectral', 'BCXmotifs'. Only for the SM method. 
#' @param CRC_method the method of calculating matrix correlation. Select
#' one of 'Coxhead' and 'Rozeboom' methods.
#' Only for the SCRC method.
#' @import SummarizedExperiment
#' @importFrom PMA CCA.permute
#' @importFrom PMA CCA
#' @importFrom energy dcor
#' @importFrom energy pdcor
#' @importFrom MatrixCorrelation RV
#' @importFrom MatrixCorrelation RV2
#' @importFrom MatrixCorrelation RVadjMaye
#' @importFrom MatrixCorrelation RVadjGhaziri
#' @importFrom MatrixCorrelation SMI
#' @importFrom MatrixCorrelation GCD
#' @importFrom MatrixCorrelation Coxhead
#' @importFrom MatrixCorrelation Rozeboom
#' @importFrom stats phyper
#' @importFrom GSEABase geneIds
#' @importFrom WGCNA corAndPvalue
#' @importFrom biclust biclust
#' @importFrom biclust binarize
#' @importFrom biclust discretize
#' @importFrom biclust BCBimax
#' @importFrom biclust BCCC
#' @importFrom biclust BCPlaid
#' @importFrom biclust BCQuest
#' @importFrom biclust BCSpectral
#' @importFrom biclust BCXmotifs
#' @importFrom biclust biclusternumber
#' @export
#' @return List object: Group competition of miRNA sponge modules,
#' and miRNA sponge modules.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_igraph <- module_igraph(ceRExp[, seq_len(10)], 
#'     mRExp[, seq_len(10)])
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_igraph_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_igraph, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Witten DM, Tibshirani R, Hastie T. A penalized matrix 
#' decomposition, with applications to sparse principal components 
#' and canonical correlation analysis. Biostatistics. 
#' 2009, 10(3):515-34.
#' @references Szekely GJ, Rizzo ML. Partial distance 
#' correlation with methods for dissimilarities. Annals of Statistics. 
#' 2014, 42(6):2382-2412.
#' @references Szekely GJ, Rizzo ML, Bakirov NK. 
#' Measuring and Testing Dependence by Correlation of Distances, 
#' Annals of Statistics, 2007, 35(6):2769-2794.
#' @references Robert P, Escoufier Y. A unifying tool for 
#' linear multivariate statistical methods: the RV-Coefficient. 
#' Applied Statistics, 1976, 25(3):257-265.
#' @references Smilde AK, Kiers HA, Bijlsma S, Rubingh CM, 
#' van Erk MJ. Matrix correlations for high-dimensional 
#' data: the modified RV-coefficient. Bioinformatics, 
#' 2009, 25(3):401-405.
#' @references Maye CD, Lorent J, Horgan GW. 
#' Exploratory analysis of multiple omics datasets using 
#' the adjusted RV coefficient". Stat Appl Genet Mol Biol., 
#' 2011, 10, 14.
#' @references EIGhaziri A, Qannari EM. Measures 
#' of association between two datasets; Application to sensory data, 
#' Food Quality and Preference, 2015, 40(A):116-124.
#' @references Indahl UG, Næs T, Liland KH. A similarity index for 
#' comparing coupled matrices. Journal of Chemometrics. 2018; 32:e3049.
#' @references Yanai H. Unification of various techniques of multivariate 
#' analysis by means of generalized coefficient of determination (GCD). 
#' Journal of Behaviormetrics, 1974, 1(1): 45-54.
#' @references Coxhead P. Measuring the relationship between two sets 
#' of variables. British journal of mathematical and statistical psychology, 
#' 1974, 27(2): 205-212.
#' @references Rozeboom WW. Linear correlations between sets of variables.
#' Psychometrika, 1965, 30(1): 57-71.
miRSM <- function(miRExp = NULL, 
                  ceRExp, 
                  mRExp = NULL, 
                  miRTarget, 
                  CandidateModulegenes,
                  typex = "standard", 
                  typez = "standard", 
                  nperms = 100, 
                  method = c("SCC", "SDC", "SRVC", "SM", "SSI", "SGCD", "SCRC"),
                  num_shared_miRNAs = 3, 
                  pvalue.cutoff = 0.05, 
                  MC.cutoff = 0.8,
                  SMC.cutoff = 0.1, 
                  RV_method = c("RV", "RV2", "RVadjMaye", "RVadjGhaziri"),
                  BCmethod = "BCPlaid",
                  CRC_method = c("Coxhead", "Rozeboom")) {

    if (method == "SCC") {
        Res <- miRSM_SCC(miRExp = miRExp, ceRExp, mRExp = mRExp, miRTarget, CandidateModulegenes,
            typex = typex, typez = typez, nperms = nperms, num_shared_miRNAs = num_shared_miRNAs,
            pvalue.cutoff = pvalue.cutoff, CC.cutoff = MC.cutoff, 
            SCC.cutoff = SMC.cutoff)
    } else if (method == "SDC") {
        Res <- miRSM_SDC(miRExp = miRExp, ceRExp, mRExp = mRExp, miRTarget, CandidateModulegenes,
            num_shared_miRNAs = num_shared_miRNAs,
            pvalue.cutoff = pvalue.cutoff, DC.cutoff = MC.cutoff, 
            SDC.cutoff = SMC.cutoff)
    } else if (method == "SRVC") {
        Res <- miRSM_SRVC(miRExp = miRExp, ceRExp, mRExp = mRExp, miRTarget, CandidateModulegenes,
            num_shared_miRNAs = num_shared_miRNAs,
            pvalue.cutoff = pvalue.cutoff, RVC.cutoff = MC.cutoff, 
            SRVC.cutoff = SMC.cutoff, RV_method = RV_method)
    } else if (method == "SM") {
        Res <- miRSM_SM(miRExp = miRExp, ceRExp, mRExp = mRExp, miRTarget, 
            BCmethod = BCmethod, pvalue.cutoff = pvalue.cutoff)
    } else if (method == "SSI") {
        Res <- miRSM_SSI(miRExp = miRExp, ceRExp, mRExp = mRExp, miRTarget, CandidateModulegenes,
                       num_shared_miRNAs = num_shared_miRNAs,
                       pvalue.cutoff = pvalue.cutoff, SI.cutoff = MC.cutoff, 
                       SSI.cutoff = SMC.cutoff)
    } else if (method == "SGCD") {
        Res <- miRSM_SGCD(miRExp = miRExp, ceRExp, mRExp = mRExp, miRTarget, CandidateModulegenes,
                       num_shared_miRNAs = num_shared_miRNAs,
                       pvalue.cutoff = pvalue.cutoff, GCD.cutoff = MC.cutoff, 
                       SGCD.cutoff = SMC.cutoff)
    } else if (method == "SCRC") {
        Res <- miRSM_SCRC(miRExp = miRExp, ceRExp, mRExp = mRExp, miRTarget, CandidateModulegenes,
                        num_shared_miRNAs = num_shared_miRNAs,
                        pvalue.cutoff = pvalue.cutoff, CRC.cutoff = MC.cutoff, 
                        SCRC.cutoff = SMC.cutoff, CRC_method = CRC_method)
    } 

    return(Res)
}


#' Inferring sample-specific miRNA sponge modules
#'
#' @title miRSM_SS
#' @param Modulelist.all List object, modules using all of samples.
#' @param Modulelist.exceptk List object, modules using all of samples excepting sample k.
#' @param sim.cutoff Similarity cutoff between modules, the interval is [0 1].
#' @param sim.method Methods for calculating similatiry between two modules, select one of three methods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @export
#' @return A list of sample-specific miRNA sponge modules
#'
#' @examples
#' data(BRCASampleData)
#' nsamples <- 3
#' modulegenes_all <- module_igraph(ceRExp[, 151:300], mRExp[, 151:300])
#' modulegenes_exceptk <- lapply(seq(nsamples), function(i) 
#'                               module_WGCNA(ceRExp[-i, seq(150)], 
#'                               mRExp[-i, seq(150)]))
#'  
#' miRSM_SRVC_all <- miRSM(miRExp, ceRExp[, 151:300], mRExp[, 151:300], 
#'                         miRTarget, modulegenes_all, 
#'                         method = "SRVC", SMC.cutoff = 0.01, 
#'                         RV_method = "RV")
#' miRSM_SRVC_exceptk <- lapply(seq(nsamples), function(i) miRSM(miRExp[-i, ], 
#'                              ceRExp[-i, seq(150)], mRExp[-i, seq(150)], 
#'                              miRTarget, modulegenes_exceptk[[i]],                                     
#'                              method = "SRVC",
#'                              SMC.cutoff = 0.01, RV_method = "RV"))
#' 
#' Modulegenes_all <- miRSM_SRVC_all[[2]]
#' Modulegenes_exceptk <- lapply(seq(nsamples), function(i) 
#'                               miRSM_SRVC_exceptk[[i]][[2]])
#' 
#' Modules_SS <- miRSM_SS(Modulegenes_all, Modulegenes_exceptk)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
miRSM_SS <- function(Modulelist.all,
                     Modulelist.exceptk,
                     sim.cutoff = 0.8,
                     sim.method = "Simpson") {
  
  if(class(Modulelist.all) != "list" | class(Modulelist.exceptk) != "list") {
    stop("Please check your input module group! The input module group should be list object! \n")
  } 
  
  Res <- lapply(seq(Modulelist.exceptk), function(i) diff_module(Modulelist.all, Modulelist.exceptk[[i]], sim.cutoff = sim.cutoff, sim.method = sim.method))
  
  names(Res) <- paste("Sample", seq(Res), sep=" ")
  return(Res)
}

#' Functional analysis of miRNA sponge modules, including functional 
#' enrichment and disease enrichment analysis
#'
#' @title module_FA
#' @param Modulelist List object: a list of miRNA sponge modules.
#' @param GOont One of 'MF', 'BP', and 'CC' subontologies.
#' @param KEGGorganism Organism, supported organism listed 
#' in http://www.genome.jp/kegg/catalog/org_list.html.
#' @param Reactomeorganism Organism, one of 'human', 'rat', '
#' mouse', 'celegans', 'yeast', 'zebrafish', 'fly'.
#' @param OrgDb OrgDb
#' @param padjustvaluecutoff A cutoff value of adjusted p-values.
#' @param padjustedmethod Adjusted method of p-values, can select 
#' one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 
#' 'fdr', 'none'.
#' @param Analysis.type The type of functional analysis selected, 
#' including 'FEA' (functional enrichment analysis) and 'DEA' 
#' (disease enrichment analysis).
#' @import org.Hs.eg.db
#' @importFrom clusterProfiler bitr
#' @importFrom clusterProfiler enrichGO
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom clusterProfiler compareCluster
#' @importFrom DOSE enrichDO
#' @importFrom DOSE enrichDGN
#' @importFrom DOSE enrichNCG
#' @importFrom ReactomePA enrichPathway
#' @export
#' @return List object: a list of enrichment analysis results.
#'
#' @examples
#' \dontrun{
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_SRVC_FEA <- module_FA(miRSM_WGCNA_SRVC_genes, Analysis.type = 'FEA')
#' miRSM_WGCNA_SRVC_DEA <- module_FA(miRSM_WGCNA_SRVC_genes, Analysis.type = 'DEA')
#' }
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Zhang J, Liu L, Xu T, Xie Y, Zhao C, Li J, Le TD (2019). 
#' “miRspongeR: an R/Bioconductor package for the identification and analysis of 
#' miRNA sponge interaction networks and modules.” BMC Bioinformatics, 20, 235.
#' @references Zhang J, Liu L, Zhang W, Li X, Zhao C, Li S, Li J, Le TD. 
#' miRspongeR 2.0: an enhanced R package for exploring miRNA sponge regulation. 
#' Bioinform Adv. 2022 Sep 2;2(1):vbac063.
#' @references Yu G, Wang L, Han Y, He Q (2012). 
#' “clusterProfiler: an R package for comparing biological themes among gene clusters.” 
#' OMICS: A Journal of Integrative Biology, 16(5), 284-287.
module_FA <- function(Modulelist, 
                      GOont = "BP", 
                      KEGGorganism = "hsa",
                      Reactomeorganism = "human", 
                      OrgDb = "org.Hs.eg.db", 
                      padjustvaluecutoff = 0.05,
                      padjustedmethod = "BH", 
                      Analysis.type = c("FEA", "DEA")) {
    
  if(class(Modulelist) != "list") {
    stop("Please check your input modules! The input modules should be list object! \n")
  } else if (class(Modulelist[[1]]) == "list"){
    Modulelist <- lapply(seq(Modulelist), function(i) unique(unlist(Modulelist[[i]])))
    names(Modulelist) <- paste("miRSM", seq_along(Modulelist), sep=" ")
  }
    
    if (Analysis.type == "FEA") {
        Res <- moduleFEA(Modulelist, ont = GOont, KEGGorganism = KEGGorganism,
            Reactomeorganism = Reactomeorganism, OrgDb = OrgDb, padjustvaluecutoff = padjustvaluecutoff,
            padjustedmethod = padjustedmethod)
    } else if (Analysis.type == "DEA") {
        Res <- moduleDEA(Modulelist, OrgDb = OrgDb, padjustvaluecutoff = padjustvaluecutoff,
            padjustedmethod = padjustedmethod)
    }

    return(Res)
}

#' Cancer enrichment analysis of miRNA sponge modules using hypergeometric distribution test
#'
#' @title module_CEA
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param Cancergenes A SummarizedExperiment object: a list of cancer genes given.
#' @param Modulelist List object: a list of the identified miRNA sponge modules. 
#' @import SummarizedExperiment
#' @importFrom stats phyper
#' @export
#' @return Cancer enrichment significance p-values of the identified miRNA sponge modules
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM.CEA.pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, 
#'                               miRSM_WGCNA_SRVC_genes)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Johnson NL, Kotz S, Kemp AW (1992) 
#' "Univariate Discrete Distributions", Second Edition. New York: Wiley.
module_CEA <- function(ceRExp, 
                       mRExp = NULL, 
                       Cancergenes,
                       Modulelist) {
  
  if(class(Modulelist) != "list") {
    stop("Please check your input modules! The input modules should be list object! \n")
  } else if (class(Modulelist[[1]]) == "list"){
    Modulelist <- lapply(seq(Modulelist), function(i) unique(unlist(Modulelist[[i]])))
    names(Modulelist) <- paste("miRSM", seq_along(Modulelist), sep=" ")
  }
  
  if(is.null(mRExp)){
  ExpData <- assay(ceRExp)  
  } else {
  ExpData <- cbind(assay(ceRExp), assay(mRExp))      
  }
  
  B <- ncol(ExpData)
  N <- length(intersect(colnames(ExpData), as.matrix(assay(Cancergenes))))
  M <- unlist(lapply(seq_along(Modulelist), function(i) length(Modulelist[[i]])))
  x <- unlist(lapply(seq_along(Modulelist), function(i) 
      length(intersect(Modulelist[[i]], as.matrix(assay(Cancergenes))))))    
  p.value <- 1 - phyper(x - 1, N, B - N, M)
  
  names(p.value) <- names(Modulelist)
  return(p.value)
}

#' Validation of miRNA sponge interactions in each miRNA sponge module
#'
#' @title module_Validate
#' @param Modulelist List object: a list of the identified miRNA sponge modules. 
#' @param Groundtruth Matrix object: a list of experimentally validated miRNA sponge interactions.
#' @export
#' @return List object: a list of validated miRNA sponge interactions in each miRNA sponge module
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' Groundtruthcsv <- system.file("extdata", "Groundtruth_high.csv", package="miRSM")
#' Groundtruth <- read.csv(Groundtruthcsv, header=TRUE, sep=",") 
#' miRSM.Validate <- module_Validate(miRSM_WGCNA_SRVC_genes, Groundtruth)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
module_Validate <- function(Modulelist, 
                            Groundtruth) {
  
  if(class(Modulelist) != "list") {
    stop("Please check your input modules! The input modules should be list object! \n")
  } else if (class(Modulelist[[1]]) == "list"){
    Modulelist <- lapply(seq(Modulelist), function(i) unique(unlist(Modulelist[[i]])))
    names(Modulelist) <- paste("miRSM", seq_along(Modulelist), sep=" ")
  }
  
  validate_res <- lapply(seq(Modulelist), function(i) 
      Groundtruth[intersect(which(as.matrix(Groundtruth[, 1]) %in% 
      Modulelist[[i]]), which(as.matrix(Groundtruth[, 2]) %in% Modulelist[[i]])), ])
  names(validate_res) <- names(Modulelist)
  
  return(validate_res)
}


#' Co-expression analysis of each miRNA sponge module and its corresponding random miRNA sponge modules
#' 
#' @title module_Coexpress
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp NULL (default) or a SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs. 
#' @param Modulelist List object: a list of the identified miRNA sponge modules. 
#' @param resample The number of random miRNA sponge modules generated, and 1000 times in default.
#' @param method The method used to evaluate the co-expression level of each miRNA sponge module.
#' Users can select "mean" or "median" to calculate co-expression value of each miRNA sponge module
#' and its corresponding random miRNA sponge module.
#' @param test.method The method used to evaluate statistical significance p-value of 
#' co-expression level higher than random miRNA sponge modules.
#' Users can select "t.test" or "wilcox.test" to calculate statistical significance p-value of 
#' co-expression level higher than random miRNA sponge modules.
#' @import SummarizedExperiment
#' @importFrom WGCNA cor
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom stats t.test
#' @importFrom stats wilcox.test
#' @export
#' @return List object: co-expression values of miRNA sponge modules and their corresponding random miRNA sponge modules, 
#' and statistical significance p-value of co-expression level higher than random miRNA sponge modules.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_Coexpress <-  module_Coexpress(ceRExp, mRExp, 
#'                                            miRSM_WGCNA_SRVC_genes, 
#'                                            resample = 10, method = "mean",
#'                                            test.method = "t.test")
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
module_Coexpress <- function(ceRExp, 
                             mRExp = NULL, 
                             Modulelist, 
                             resample = 1000, 
                             method = c("mean", "median"),
                             test.method = c("t.test", "wilcox.test")) {
  
  if(class(Modulelist) != "list") {
    stop("Please check your input modules! The input modules should be list object! \n")
  } 
  
  if(is.null(mRExp)){
    ceRExp <- assay(ceRExp)
    module_ceRExp1 <- lapply(seq_along(Modulelist), function(i) 
      ceRExp[, which(colnames(ceRExp) %in% Modulelist[[i]][[1]])])
    module_ceRExp2 <- lapply(seq_along(Modulelist), function(i) 
      ceRExp[, which(colnames(ceRExp) %in% Modulelist[[i]][[2]])])
    
    if (method == "mean"){
      module_avg_cor <- unlist(lapply(seq_along(Modulelist), function(i) 
        mean(abs(cor(module_ceRExp1[[i]], module_ceRExp2[[i]])))))
    } else if (method == "median"){
      module_avg_cor <- unlist(lapply(seq_along(Modulelist), function(i) 
        median(abs(cor(module_ceRExp1[[i]], module_ceRExp2[[i]])))))
    }
    
    module_avg_cor_resample <- c()
    module_avg_cor_pvalue <- c()
    for (i in seq_along(Modulelist)){
      temp1 <- replicate(resample, sample(seq_len(ncol(ceRExp)), size = ncol(module_ceRExp1[[i]])))
      temp2 <- replicate(resample, sample(seq_len(ncol(ceRExp)), size = ncol(module_ceRExp2[[i]])))
      module_ceRExp1_resample <- lapply(seq_len(resample), function(i) ceRExp[, temp1[, i]])
      module_ceRExp2_resample <- lapply(seq_len(resample), function(i) ceRExp[, temp2[, i]])
      
      if (method == "mean"){
        temp <- unlist(lapply(seq_len(resample), function(i) 
          mean(na.omit(abs(cor(module_ceRExp1_resample[[i]], module_ceRExp2_resample[[i]]))))))
        module_avg_cor_resample[i] <- mean(temp)
        if (test.method == "t.test"){
          module_avg_cor_pvalue[i] <- t.test(temp, alternative = "less", mu = module_avg_cor[i])$p.value
        } 
        if (test.method == "wilcox.test"){
          module_avg_cor_pvalue[i] <- wilcox.test(temp, alternative = "less", mu = module_avg_cor[i])$p.value
        }
      } else if (method == "median"){
        temp <- unlist(lapply(seq_len(resample), function(i) 
          median(na.omit(abs(cor(module_ceRExp1_resample[[i]], module_ceRExp2_resample[[i]]))))))
        module_avg_cor_resample[i] <- median(temp)
        if (test.method == "t.test"){
          module_avg_cor_pvalue[i] <- t.test(temp, alternative = "less", mu = module_avg_cor[i])$p.value
        } 
        if (test.method == "wilcox.test"){
          module_avg_cor_pvalue[i] <- wilcox.test(temp, alternative = "less", mu = module_avg_cor[i])$p.value
        }
      }
    }
  } else {
  ceRExp <- assay(ceRExp)
  mRExp <- assay(mRExp)
  module_ceRExp <- lapply(seq_along(Modulelist), function(i) 
    ceRExp[, which(colnames(ceRExp) %in% Modulelist[[i]][[1]])])
  module_mRExp <- lapply(seq_along(Modulelist), function(i) 
    mRExp[, which(colnames(mRExp) %in% Modulelist[[i]][[2]])])
  
  if (method == "mean"){
    module_avg_cor <- unlist(lapply(seq_along(Modulelist), function(i) 
      mean(abs(cor(module_ceRExp[[i]], module_mRExp[[i]])))))
  } else if (method == "median"){
    module_avg_cor <- unlist(lapply(seq_along(Modulelist), function(i) 
      median(abs(cor(module_ceRExp[[i]], module_mRExp[[i]])))))
  }
  
  module_avg_cor_resample <- c()
  module_avg_cor_pvalue <- c()
  for (i in seq_along(Modulelist)){
    temp1 <- replicate(resample, sample(seq_len(ncol(ceRExp)), size = ncol(module_ceRExp[[i]])))
    temp2 <- replicate(resample, sample(seq_len(ncol(mRExp)), size = ncol(module_mRExp[[i]])))
    module_ceRExp_resample <- lapply(seq_len(resample), function(i) ceRExp[, temp1[, i]])
    module_mRExp_resample <- lapply(seq_len(resample), function(i) mRExp[, temp2[, i]])
    
      if (method == "mean"){
          temp <- unlist(lapply(seq_len(resample), function(i) 
            mean(na.omit(abs(cor(module_ceRExp_resample[[i]], module_mRExp_resample[[i]]))))))
          module_avg_cor_resample[i] <- mean(temp)
          if (test.method == "t.test"){
          module_avg_cor_pvalue[i] <- t.test(temp, alternative = "less", mu = module_avg_cor[i])$p.value
          } 
          if (test.method == "wilcox.test"){
          module_avg_cor_pvalue[i] <- wilcox.test(temp, alternative = "less", mu = module_avg_cor[i])$p.value
          }
      } else if (method == "median"){
        temp <- unlist(lapply(seq_len(resample), function(i) 
          median(na.omit(abs(cor(module_ceRExp_resample[[i]], module_mRExp_resample[[i]]))))))
        module_avg_cor_resample[i] <- median(temp)
        if (test.method == "t.test"){
          module_avg_cor_pvalue[i] <- t.test(temp, alternative = "less", mu = module_avg_cor[i])$p.value
        } 
        if (test.method == "wilcox.test"){
          module_avg_cor_pvalue[i] <- wilcox.test(temp, alternative = "less", mu = module_avg_cor[i])$p.value
        }
      }
  }
  }
  
  module_coexpress <- list(module_avg_cor, module_avg_cor_resample, module_avg_cor_pvalue)
  names(module_coexpress) <- c("Real miRNA sponge modules", "Random miRNA sponge modules", "Statistical significance p-value")
  return(module_coexpress)
}

#' Extract common miRNAs of each miRNA sponge module
#' 
#' @title share_miRs
#' @param miRExp NULL (default) or a SummarizedExperiment object. miRNA expression data: 
#' rows are samples and columns are miRNAs.
#' @param miRTarget A SummarizedExperiment object. Putative 
#' miRNA-target binding information.
#' @param Modulelist List object: a list of the identified miRNA sponge modules.
#' @import SummarizedExperiment
#' @export
#' @return List object: a list of common miRNAs of each miRNA sponge module.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, miRTarget, miRSM_WGCNA_SRVC_genes)
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
share_miRs <- function(miRExp = NULL,
                       miRTarget, 
                       Modulelist){
  
  if(class(Modulelist) != "list") {
    stop("Please check your input modules! The input modules should be list object! \n")
  } 
  
  if(is.null(miRExp)){
    miRTarget <- assay(miRTarget)
  } else {
    miRNames <- colnames(miRExp)
    miRTarget <- assay(miRTarget)
    miRTarget <- miRTarget[which(miRTarget[, 1] %in% miRNames), ]  
  }
    
    Res <- list()
    for (i in seq_along(Modulelist)){
      modulegenes <- Modulelist[[i]]
      tmp1 <- unique(miRTarget[which( miRTarget[, 2] %in% 
                                        modulegenes[[1]] ), 1])
      tmp2 <- unique(miRTarget[which( miRTarget[, 2] %in% 
                                        modulegenes[[2]] ), 1])
      tmp3 <- intersect( tmp1, tmp2 )
      Res[[i]] <- tmp3
    }
  
  names(Res) <- names(Modulelist)
  return(Res)
}

#' miRNA distribution analysis of sharing miRNAs by the identified miRNA sponge modules
#' 
#' @title module_miRdistribute
#' @param share_miRs List object: a list of common miRNAs of each miRNA sponge module 
#' generated by share_miRs function. 
#' @export
#' @return Matrix object: miRNA distribution in each miRNA sponge module.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, miRTarget, miRSM_WGCNA_SRVC_genes)
#' miRSM_WGCNA_miRdistribute <- module_miRdistribute(miRSM_WGCNA_share_miRs)
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
module_miRdistribute <- function(share_miRs) {
  
  miRs <- unique(unlist(share_miRs))
  res <- NULL
  interin <- NULL
  for (i in seq_along(miRs)) {
    for (j in seq_along(share_miRs)) {
      if (length(which(miRs[i] %in% share_miRs[[j]]) == 1)) {
        interin <- c(interin, names(share_miRs)[j])
      }
    }
    res1 <- paste(interin, collapse = ", ")        
    res2 <- length(interin)
    res <- rbind(res, c(miRs[i], res1, res2))
    interin <- NULL
  }
  colnames(res) <- c("miRNA", "Module ID", "Number of modules")
    
  return(res)
}

#' Extract miRNA-target interactions of each miRNA sponge module
#' 
#' @title module_miRtarget
#' @param share_miRs List object: a list of common miRNAs of each miRNA sponge module 
#' generated by share_miRs function. 
#' @param Modulelist List object: a list of the identified miRNA sponge modules.
#' @export
#' @return List object: miRNA-target interactions of each miRNA sponge module.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, miRTarget, miRSM_WGCNA_SRVC_genes)
#' miRSM_WGCNA_miRtarget <- module_miRtarget(miRSM_WGCNA_share_miRs, 
#'                                           miRSM_WGCNA_SRVC_genes)
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
module_miRtarget <- function(share_miRs, 
                             Modulelist){
  
  if(class(Modulelist) != "list") {
    stop("Please check your input modules! The input modules should be list object! \n")
  } else if (class(Modulelist[[1]]) == "list"){
    Modulelist <- lapply(seq(Modulelist), function(i) unique(unlist(Modulelist[[i]])))
    names(Modulelist) <- paste("miRSM", seq_along(Modulelist), sep=" ")
  }
  
  res_int <- list()
  for (k in seq(share_miRs)){
      CommonmiRs <- share_miRs[[k]]
      targets <- Modulelist[[k]]
      len_CommonmiRs <- length(CommonmiRs)
      len_targets <- length(targets)
      res_interin <- matrix(NA, len_CommonmiRs*len_targets, 2)
      for (i in seq_len(len_CommonmiRs)){
          for (j in seq_len(len_targets)){
              res_interin[(i-1)*len_targets+j, 1] <- CommonmiRs[i]
              res_interin[(i-1)*len_targets+j, 2] <- targets[j]
          }
      }
  res_int[[k]] <- res_interin
  }
  names(res_int) <- names(Modulelist)
  return(res_int)
}

#' Extract miRNA sponge interactions of each miRNA sponge module
#'
#' @title module_miRsponge
#' @param Modulelist List object: a list of the identified miRNA sponge modules.
#' @export
#' @return List object: miRNA sponge interactions of each miRNA sponge module.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_miRsponge <- module_miRsponge(miRSM_WGCNA_SRVC_genes)
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
module_miRsponge<- function(Modulelist){
  
    if(class(Modulelist) != "list") {
      stop("Please check your input modules! The input modules should be list object! \n")
    }
  
    res_int <- list()
    
      for (k in seq(Modulelist)){
        ceRNA1 <- Modulelist[[k]][[1]]
        ceRNA2 <- Modulelist[[k]][[2]]
        len_ceRNA1 <- length(ceRNA1)
        len_ceRNA2 <- length(ceRNA2)
        res_interin <- matrix(NA, len_ceRNA1*len_ceRNA2, 2)
        for (i in seq_len(len_ceRNA1)){
          for (j in seq_len(len_ceRNA2)){
            res_interin[(i-1)*len_ceRNA2+j, 1] <- ceRNA1[i]
            res_interin[(i-1)*len_ceRNA2+j, 2] <- ceRNA2[j]
          }
        }
        res_int[[k]] <- res_interin     
      }
    
    names(res_int) <- names(Modulelist)
  return(res_int)
}
