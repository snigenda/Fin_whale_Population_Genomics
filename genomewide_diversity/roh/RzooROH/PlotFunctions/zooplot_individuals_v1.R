#'Plot individual curves with proportion of the genome in each HBD class or cumulated
#'proportion in HBD classes with rates smaller than a threshold.
#'
#'For each individual, the function plots the mean percentage of the genome in different
#' HBD classes or the inbreeding coefficient obtained by summing autozygosity associated
#' with HBD classes with a rate lower or equal to a threshold (e.g., including all HBD
#'classes with longer and more recent HBD segments than a selected threshold).
#'
#'@param input a named list with one or several zres objects obtained after
#'  running zoorun. The zres objects are the output of the zoorun function. For
#'  instance, putting list(name1 = zres1, name2 = zres2). The function will then
#'  use the names in the plot (in case several zres objects are used).
#'
#'@param cumulative a logical indicating whether individual autozygosity is plotted
#'  per class (FALSE) or summed over all HBD class with a rate smaller than a
#'  value (these cumulated values are obtained for every rate defined in the
#'  model). By default, this value is TRUE. When FALSE, the percentages correspond
#'  to the individual genome-wide probabilities of belonging to each HBD-class
#'  or to the fraction of the genome in an autozygosity class. When TRUE, we obtain
#'  the probability of belonging to an HBD class with a rate smaller or equal than
#'  a threshold (here we use the pre-defined rates of the model as thresholds), averaged
#'  over the whole genome for each individual. This corresponds to report individual genomic
#'  inbreeding coefficients estimated with respect to different base populations obtained
#'  by selecting different thresholds T that determine which HBD classes are considered
#'  in the estimation of the genomic inbreeding coefficient (setting the base population
#'  approximately 0.5 * T generations ago).
#'
#'@param toplot A list of vectors indicating the zres@@ids to be plotted. This
#'  option can be used to select the individuals to plot. The list must contain
#'  one vector per population or zres object. By default, all individuals are
#'  plotted.
#'
#'@param ncols when several populations are plotted, ncols determines how many results (graphs)
#'are plotted per row.
#'
#'@return The function plots either the individual proportions of the genome associated with
#'different HBD classes or individual genomic inbreeding coefficients estimated with respect
#'to different base populations (from young to older). With both option, the average values are
#'plotted in red.
#'
#'@export

zooplot_individuals_v1 <- function (input,  cumulative=TRUE, toplot=NULL, ncols=2, col) {
  layout ( matrix (1:2, nrow=1), widths=c(0.90, 0.1))
  par (mar =c(1, 0, 4, 3))
  par (oma=c(6, 6, 0, 0))
  
  if (is (input, "list")) {
    if(any (lapply (input, class) != "zres")) {
      stop ("Some objects are NOT of class  \"zres\"\n")
    }
  }else {
    if(is (input,"zres")){
      input <- list(input)
    }
    else {
      stop ("input should be a list\n")
    }}
  
  if(length(names(input))==0 & length(input) > 1){
    warning("No names were provided for the input list!\n We will use capital letters.\n")
    names(input)=LETTERS[1:length(input)]}
  
  
  ks <- c()
  for (i in 1:length (input)) {
    myres <- input[[i]]@realized
    ks [i] <- ncol (myres)
    myres <- myres [, -c(ncol(myres))]
    myres <- data.frame (id=input[[i]]@ids, fullid=input[[i]]@sampleids, myres)
    
    
    myres$pop <- i
    if (i ==1) {
      allres <- myres
    }else {
      allres <- rbind(allres, myres)
    }
  }
  
  if (length (unique (ks)) >1) {
    stop  ("different models used for the data\n")
  }else {
    k <- unique (ks) -1
  }
  
  
  if (cumulative ==TRUE) {
    xlab=expression("Value of the threshold T used to estimate F"['G-T'] * " (using HBD classes with R"['K']<=" T)")
    ylab =expression ("Genomic inbreeding coefficient" ~ ( F [G-T]  ))
    ylim <- c(0, max(apply (allres [, 3:(k+2)], 1, cumsum))) * 1.06
    
  }else {
    xlab="Rate of the HBD class"
    ylab="Proportion of the genome in HBD class"
    ylim <- c(0,  max(allres [, 3:(k+2)])) *1.06
  }
  
  npop = length(input)
  if(npop>1 & npop%%2==1){npop=npop+1}
  mymat <- matrix (1:npop, ncol=ncols, byrow=TRUE)
  layout (mymat)
  
  for (j in 1:length (input)) {
    myres <- allres [allres$pop ==j, ]
    mymean <- apply (myres[, 3:(k+2)] , 2, mean)
    if (!is.null (toplot)) {
      if (length (toplot) != length (input)) {
        stop ('length of topot should match the length of input\n')
      }
      if (sum (myres$id %in% toplot [[j]]) == length (toplot [[j]])  ) {
        myres <- myres [myres$id %in% toplot [[j]], ]
        cat ('population ',  j,  ": ", nrow (myres),  "\n")
      }else {
        warning ("some ids for data ",  j,  " was not found\n")
      }
    }else {
      if (j ==1) {
        warning ("\nAll individuals are plotted; use toplot to select individuals\n")
      }
    }
    
    
    ##mymean <- apply (myres [, 3:(2+k)] , 2, mean)

    
    myres2 <- as.data.frame(myres)
    myres2$id <- NULL
    myres2$pop <- NULL
  
    mymean2 <- melt(cumsum(mymean))
    mymean2$krates <- input[[1]]@krates[1, 1:9]
    mymean3 <- cbind("Mean","mean", mymean2, "Mean", 0.9)
    colnames(mymean3) <- c("fullid", "variable", "value", "krates", "group", "size")
    
    if (cumulative==TRUE) {
      myres3 <- apply(myres2[2:10],1,cumsum)
      myres4 <- cbind(myres2[1],t(myres3))
      myres5 <- melt(myres4)
      myres5$krates <- rep(input[[1]]@krates[1, 1:9], each = length(myres2$fullid))
      myres5$group <- "Samples"
      myres5$size <- 0.8
      myres6 <- rbind(myres5, mymean3)
      p2 <- ggplot(myres6, aes(x= as.factor(krates), y= value, group=fullid, color = group, size = as.factor(size))) + geom_line() + 
        theme_light() +labs(title = "Proportion of the genome associated with different HBD classes", 
                                            x = xlab, y = ylab, color = "") + scale_color_manual(values = c(col,"#000000")) +
        guides(size = FALSE) + scale_size_manual(values= c(0.7,1.8))
      print(p2)
    }else {
      myres3 <- melt(myres2)
      myres3$krates <- rep(input[[1]]@krates[1, 1:9], each = length(myres2$fullid))
      p2 <- ggplot(myres3, aes(x= as.factor(krates), y= value, group=fullid)) + geom_line(size=1) + 
        geom_point() + theme_light() + labs(title = "Proportion of the genome associated with different HBD classes (individuals)", 
                                            x = xlab, y = ylab) 
      print(p2)
    }
  
  }
  
  
}