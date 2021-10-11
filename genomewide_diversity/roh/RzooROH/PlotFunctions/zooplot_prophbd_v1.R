#'Plot proportion of the genome associated with different HBD classes
#'
#'Plot the mean percentage of the genome in different HBD classes or the
#'inbreeding coefficient obtained by summing autozygosity associated with HBD
#'classes with a rate lower or equal to a threshold (e.g., including all HBD
#'classes with longer and more recent HBD segments than a selected threshold).
#'
#'@param input a named list with one or several zres objects obtained after
#'  running zoorun. The zres objects are the output of the zoorun function. For
#'  instance, putting list(name1 = zres1, name2 = zres2). The function will then
#'  use the names in the plot (in case several zres objects are used).
#'
#'@param cols a vector with the colors to be used for each population or zres
#'  object.
#'
#'@param style select "barplot", "lines" or "boxplot" for the graphic styles. Boxplot
#'can be used with a single zres file or population.
#'
#'@param cumulative a logical indicating whether mean autozygosity is estimated
#'  per class (FALSE) or summed over all HBD class with a rate smaller than a
#'  value (these cumulated values are obtained for every rate defined in the
#'  model). By default, this value is FALSE. When FALSE, the percentages correspond
#'  to the mean individual genome-wide probabilities of belonging to each HBD-class
#'  or to the fraction of the genome in an autozygosity class. When TRUE, we obtain
#'  the mean probability of belonging to an HBD class with a rate smaller or equal than
#'  a threshold (here we use the pre-defined rates of the model as thresholds), averaged
#'  over the whole genome and all individuals. This corresponds to report mean genomic
#'  inbreeding coefficients estimated with respect to different base populations obtained
#'  by selecting different thresholds T that determine which HBD classes are considered
#'  in the estimation of the genomic inbreeding coefficient (setting the base population
#'  approximately 0.5 * T generations ago).
#'
#'@return The function plots either the average proportion of the genome associated with
#'different HBD classes or the average genomic inbreeding coefficient estimated with respect
#'to different base populations (from young to older).
#'
#'@export

zooplot_prophbd_v1 <- function (input,  cols=NULL, style="barplot", cumulative = FALSE) {
  if(length(input) > 1){
    layout ( matrix (1:2, nrow=1), widths=c(0.85, 0.15))
    par (mar =c(6, 6, 2, 0))
  }
  if(length(input) == 1){par (mar =c(6, 6, 2, 2))}
  
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
  
  if (is.null (cols)) {
    if (length (input) >8) {
      stop ("Please provide ",  length (input), " colors,  There are only 8 default colors\n")
    }
    cols <- brewer.pal(8, "Dark2") [1:length (input)]
  }
  ks <- c()
  means <- list ()
  for (i in 1:length (input)) {
    myres <- input[[i]]@realized
    myres <- myres [, -c(ncol(myres))]
    myids <- input [[i]]@sampleids
    myres <- data.frame (id=myids, myres)
    means [[i]] <- apply (myres[, -1], 2, mean)
    ks [i] <- ncol (input [[1]]@krates)
    
    myres$pop <- i #if ind.results are used
    if (i ==1) {
      allres <- myres
    }else {
      allres <- rbind(allres, myres)
    }
    
  }
  
  if (style =="barplot") {
    if (length (unique (ks)) ==1  ) {
      k <- unique (ks) -1
      allmeans <- matrix (nrow=length (input), ncol=k) #means / cum.means
      for (j in 1:length (input)) {
        if (cumulative==TRUE) {
          xlab=expression("Value of the threshold T used to estimate F"['G-T'] * " (using HBD classes with R"['K']<=" T)")
          ylab =expression ("Genomic inbreeding coefficient" ~ ( F [G-T]  ))
          main = "Proportion of the genome associated with different HBD classes"
          sub ="Cumulative proportions"
          allmeans [j, ] <- cumsum(means [[j]]) [1:k]
        }else {
          allmeans [j, ] <-  means [[j]] [1:k]
          ylab="Proportion of the genome in HBD class"
          xlab=expression("Rate R"[k] * " of the HBD class")
          main ="Proportion of the genome associated with different HBD classes "
          sub = ""
        }
      }
    }else {
      stop ('Different models used for different data sets \n')
    }
    
    rownames (allmeans) <- names (input)
    colnames (allmeans) <- input[[1]]@krates[1, 1:k]
    
    
    ylim <-  c(0, (max(allmeans)*1.01))
    #mybp <- barplot (allmeans, beside=TRUE, col=cols,cex.lab=1, cex.axis=1, xlab=xlab, ylab=ylab, bty = "o",
     #                xaxs="i", yaxs="i",  ylim=ylim, width=4, xaxt="n", main = main,  family = "arial", sub = sub)
    
    
    allmeans2 <- melt(as.data.frame(allmeans))
    allmeans2$Groups <- names (input)
    ggplot(allmeans2, aes(x= variable, y = value, fill = Groups)) + geom_bar(stat="identity", position=position_dodge()) +
      scale_fill_brewer(palette="Dark2") + theme_light() + labs(title= main, subtitle = sub, y = ylab, x=xlab)
    ##abline (h=par ("usr") [3], col='black', lwd=2)
    #xats <- c(-2, apply (mybp, 2, mean))
    #axis (side =1, at=xats,  labels=c("", input[[1]]@krates [1, 1:k]), cex.axis=1, cex.lab=1)
    
    
    
  }else if (style =="lines") {
    if (length (unique (ks)) ==1  ) {
      k <- unique (ks) -1
      if (cumulative==TRUE) {
        xlab=expression("Value of the threshold T used to estimate F"['G-T'] * " (using HBD classes with R"['K']<=" T)")
        ylab =expression ("Genomic inbreeding coefficient" ~ ( F [G-T]  ))
        main = "Proportion of the genome associated with different HBD classes "
        sub = "Cumulative proportions"
        cummeans <- lapply (means, cumsum)
        ylim <- range(unlist(lapply (cummeans, range)))
      }else {
        xlab=expression("Rate R"[k] * " of the HBD class")
        ylab="Proportion of the genome in HBD class"
        ylim <- range(unlist(lapply (means, range)))
        main = "Proportion of the genome associated with different HBD classes"
        sub = ""
      }
      
      #plot (1,  type="n",xlim=c(1, k), ylim=ylim,  xaxt="n", cex.lab=2, cex.axis=1.5, bty="l", xlab="", ylab="", main = main)
      #mtext(side=1, xlab, line=4, cex=1)
      #mtext(side=2, ylab, line=4, cex=1)
      
      #axis (side=1, at=1:k, labels=input[[1]]@krates [1, 1:k], cex.axis=1.5, cex.lab=2)
      #for (i in 1:length (input)){
        #if (cumulative==TRUE) {
          #lines (1:k, cummeans [[i]], lwd=3, col=cols [i], type="b", pch=19)
        #}else {
          #lines (1:k, means [[i]], lwd=3, col=cols [i], type="b", pch=19)
        #}
      #}
      
      names(cummeans) <-  names (input)
      cummeans2 <- as.data.frame(cummeans)
      #cummeans2$Groups <- names (input)
      cummeans2$krates <- input[[1]]@krates[1, 1:k]
      cummeans3 <- melt(cummeans2, id = "krates")
      
      ggplot(cummeans3, aes(x= as.character(krates), y= value, group=variable)) + geom_line(aes(color=variable), size=1.2) + 
        geom_point(aes(color=variable)) + theme_light() + scale_color_brewer(palette="Dark2") +
        labs(title = main, subtitle = sub, x = xlab, y = ylab, color = "Groups") +
        scale_x_discrete(limits = c("2","4","8","16","32","64","128","256","512"))
      
      
      
    }else {
      stop ('Different models used for different data sets \n')
    }
  }else if (style =="boxplot"){
    if(length(input) > 1){stop('Style boxplot can be used only with one results file. \n')}
    if (length (unique (ks)) ==1  ) {
      k <- unique (ks) -1
      prophbd <- matrix (nrow=input[[1]]@nind, ncol=k) #means / cum.means
      if (cumulative==TRUE) {
        xlab=expression("Value of the threshold T used to estimate F"['G-T'] * " (using HBD classes with R"['K']<=" T)")
        ylab =expression ("Genomic inbreeding coefficient" ~ ( F [G-T]  ))
        
        prophbd <- t(apply(input[[1]]@realized[,(1:k)],1,cumsum))
      }else {
        prophbd <- input[[1]]@realized[,(1:k)]
        ylab="Proportion of the genome in HBD class"
        xlab=expression("Rate R"[k] * " of the HBD class")
      }
    }else {
      stop ('Different models used for different data sets \n')
    }
    
    colnames (prophbd) <- input[[1]]@krates[1, 1:k]
    
    ylim <-  c(0, (max(prophbd)*1.05))
    #boxplot (prophbd, col=cols,cex.lab=1.5, cex.axis=1.25, xlab=xlab, ylab=ylab, ylim=ylim,
             #names = input[[1]]@krates [1, 1:k])
    
    prophbd2 <- as.data.frame(prophbd)
    prophbd3 <- melt(prophbd2)
    
    p <- ggplot(prophbd3, aes(x=as.factor(variable), y=value)) + geom_boxplot(fill = cols,notch=TRUE) + 
      geom_jitter(shape=16, position=position_jitter(0.2)) + theme_light()+
      labs(title = "Proportion of the genome associated with different HBD classes", x = xlab, y = ylab)
    print(p)
  } else {
    stop ('style not recognized\n')
  }
  
  #if(length(input) > 1){
    #par (mar =c (0, 0, 0, 0) )
    #plot (1,  type="n", xlab="", ylab="", axes=FALSE)
    #if(!is.null(names(input))){
      #legend ('topright', legend=names(input), col=cols, pch=15, cex=1,
       #       pt.cex=1,  bty="n", xjust=0, y.intersp=1.5, text.font=1 )}
    #if(is.null(names(input))){
      #warning("You can add names to your list to improve the legend.")
      #legend ('topright', legend=1:length(input), col=cols, pch=15, cex=1,
        #      pt.cex=1,  bty="n", xjust=0, y.intersp=1.5, text.font=1)}
  #}
  
}