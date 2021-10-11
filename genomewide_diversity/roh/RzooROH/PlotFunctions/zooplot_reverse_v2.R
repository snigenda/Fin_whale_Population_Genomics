#'Plot the partitioning of the genome in different HBD classes for each individual
#'
#'@param input a named list with one or several zres objects obtained after
#'  running zoorun. The zres objects are the output of the zoorun function. For
#'  instance, putting list(name1 = zres1, name2 = zres2). The function will then
#'  use the names in the plot (in case several zres objects are used).
#'
#'@param cols A vector with the colors to be used for each class in the model.
#'
#'@param plotids A logical indicating whether the IDs of the individuals are
#'  plotted on the graph (TRUE by default).
#'
#'@param toplot A list of vectors indicating the zres@@ids to be plotted. This
#'  option can be used to select the individuals to plot. The list must contain
#'  one vector per population or zres object. By default, all individuals are
#'  plotted.
#'
#'@param randomids A logical indicating whether a randomset of individuals is
#'  plotted. This option allows to reduce the number of individuals in the plot.
#'  The option can not be used simultaneously with the toplot option. By
#'  default, randomids is FALSE.
#'
#'@param nrandom A vector indicating the number of individuals to be randomly
#'  sampled per population or per zres object when randomids is TRUE. By
#'  default, we select 10 individuals per zres object. This vector must have the
#'  same length as the input list.
#'
#'@param seed A value for the random seed used to sample individuals to plot
#'  (when the randomids option is TRUE).
#'
#'@param ylim The limits of the y-axis.
#'
#'@param border Whether a border is plotted around each block of the barplot
#'or not. When set to FALSE, it allows to get a less dense plot when many
#'individuals are plotted.

#'@param nonhbd Whether the a border is plotted around the non-hbd contribution.
#'When set to FALSE, it allows to get a less dense plot when many individuals
#'are plotted.
#'
#'@param vertical Whether the populations or zres labels are printed vertically
#' or not.
#'
#'@return Individuals are presented with stacked barplots. Each vertical stack of bars
#'represents one individual. Each class is represented with a bar of a different
#'color. The height of the bar represents the proportion associated
#'with the corresponding class. The total height of the stack is the total
#'autozygosity.
#'
#'@export


zooplot_partitioning_v1 <- function (input, cols=NULL, plotids=TRUE, toplot=NULL,
                                     randomids=FALSE, nrandom=NULL, seed=100, ylim=c(0, 1),
                                     border=TRUE, nonhbd=TRUE, vertical=FALSE){
  
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
    allcols <-  c("#11A4C8","#63C2C5","#1D4F9F","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#7E277C","#F7EC16","#F8941E")
  }else {
    allcols <- cols
  }
  set.seed=seed
  ns <- c()
  ks <- c()
  allids <- c()
  
  for (i in  1:length (input)) {
    myres <- input [[i]]@realized
    ks[i] <- ncol(input[[i]]@krates)
    
    if (randomids ==FALSE ) {
      if (!is.null (nrandom)) {
        stop ("Choose only one option, either toplot, either randomids.\n")
      }
      if ( !is.null (toplot)) {
        if (length(toplot)==length(input)) { ##
          mysample <- input [[i]]@ids %in% toplot [[i]]
          if (sum (mysample) == length (toplot [[i]])){
            myres <- myres [mysample, ]
            myids <- input [[i]]@sampleids [mysample]
          }else {
            stop ("Results not found for ids provided for data ",  i,  "\n")
          }
        }else {
          stop ("toplot must be a list with the same length as the input list.\n")
        }
      }else {
        ##all default
        myres <- myres
        myids <- input[[i]]@sampleids
      }
    }else if (randomids ==TRUE) {
      if (!is.null (toplot)) {stop ("Choose only one option, either toplot, either randomids.\n")}
      if(!is.null (nrandom)  & (length(nrandom) == length(input)   )  ){
        if (i ==1) {warning ("\nRandom seed ", seed,  " is used to sample individuals.\n")}
        if (length (input[[i]]@ids) >= nrandom [i]) {
          mysample <- sample (input[[i]]@ids, nrandom [i])
          myres <- myres [mysample, ]
          myids <- input [[i]]@sampleids[mysample]
        }else {
          stop ("nrandom for populaltion ",  i,  " larger than the number of samples\n")
        }
      }else {
        stop ("When randomids is TRUE, a VECTOR of length equal to the input list  must be provided with the option nrandom.\n")
      }
    }
    
    ns [i] <- length (myids)
    rownames (myres) <- myids
    
    if (nonhbd==FALSE) {
      myres <- myres [, -ncol(myres)]
    }
    myres <- t (myres)
    if (i==1) {
      allres <- myres
    }else {
      allres <- cbind (allres, myres)
    }
  } #FOR I
  
  
  if (length (input) >1  &  length (unique (ks)) >1  ) {
    stop ("different models used for the datas\n")
  }else {
    k <- unique (ks)
    if (k >14 | length (allcols) <k  ) {
      stop ('Please Provide ',  k,  'colors\n')
    }else {
      if (k==14) {
        cols <- allcols
      }else if(k < 8){
        cols <- allcols [3:(k+2)]
        cols[k] <- allcols[14]
      }else if (input[[1]]@krates[1,1] < 4) {
        cols <- allcols [1:k]
        cols[k] <- allcols[14]
      }else if (input[[1]]@krates[1,1] >= 4) {
        cols <- allcols [2:(k+1)]
        cols[k] <- allcols[14]
      }
    }
  }
  
  if (length (input) >1) {
    
    allres2 <- as.data.frame(allres)
    allres2$krates <- input[[1]]@krates[1, 1:9]
    allres3 <- melt(allres2, id = "krates")
    allres3$Groups <- rep(names(input), c(270,180))
    
    ggplot(allres3, aes(x = variable, y = value, fill = as.factor(krates))) + 
      geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ Groups, scales = "free") + theme_light() +
      theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      labs(title = "Partitioning individual genomes in different HBD classes", x = "Samples", y = "", fill = "HBD classes") +
      scale_fill_viridis(discrete = TRUE)
    
    
  }else {
    
    allres2 <- as.data.frame(allres)
    allres2$krates <- input[[1]]@krates[1, 1:9]
    allres3 <- melt(allres2, id = "krates")
    ggplot(allres3, aes(x = variable, y = value, fill = as.factor(krates))) + 
      geom_bar(stat = 'identity', position = 'stack') + theme_light()  + theme(axis.text.x= element_text(angle = 90)) +
      labs(title = "Partitioning individual genomes in different HBD classes", x = "Samples", y = "", fill = "HBD classes") +
      scale_fill_viridis(discrete = TRUE)
    
  }
  
  
}
