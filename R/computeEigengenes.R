computeEigengenes <- function(genExpr=NULL, eigenloci=NULL, netPath, geNames,
                                Labels, Label1, Label2, mus, combiningMu=NA,
                                survival, event="Dead", doIgnoreNas=FALSE, 
                                mu2modules, doWarn=TRUE, verbose=0){
    ## Computes eigengenes (features) for expr values and dnam levels.
    ## Input:
    ## genExpr: A matrix which rows are genes, cols are patients and values are gene expression values.
    ## eigenloci: A matrix which rows are patients, cols are genes and values are methylation levels 
    ##^of genes (computed using clustering and PCA).
    ## netPath:outPath beshe Path to the combinedNetwork object, e.g. networkPath[[tIndx]].
    ## geNames: Genes that 1) both their methylation level and gene expression values are available 2) either
    ##their expression value or metylation level highly correlate with the survival time of dead patients.
    ## Labels: A character names vector that determines the risk group for each patient.
    ## Names are patients' names and values are their related risk category. 
    ## mus: A single value or a vector of mu values used to construct the networks.
    ## combiningMu: A numeric vector determining the l value(s) used for eigengene computation,
    ## i.e., (1-l)*expr + l*dnam.
    ##^If set to NA, the same value that was used for network construction, will be
    ## be used for eigengene computation.
    ## survival: A matrix which rows are patinets' name and columns are  “Time” (to event; survival time) and
    ## event (which is "Dead" by default), where 1 means the corresponding patient was dead, 
    ##and 0 means was alive at the last time of contact.
    ## Output:
    ## pigengenes: pigengene object for each mu value. 
    ## eigengenes: Inferred values 
    ## eigenCor: A named vector of correlation of eigengenes with the survival time of dead patientspatients.
    result <- list()
    ## QC:
    if(is.null(genExpr)& is.null(eigenloci))
        stop("Both genExpr and eigenloci can't be null!")
    if(!file.exists(netPath))
        stop("netPath is not a file path!")
    for(l1 in c(Label1, Label2)){
        if(! l1 %in% Labels)
            stop("There is no patient labeled as: ", l1)
    }
    if(is.null(geNames))
        stop("geNames cannot be NULL!")

    ## lCombine: 3 different modes to combine expr and dnam data
    ## The eigengenes will be computed based on (1-lCombine)*expr+lCombine*dnam
    ## e.g., if lCombine=c("e"=0, "m"=1, "em"=0.6), the above formula will be computed 3 times.
    lCombine <- c()
    ## If in the next step of the analysis, the eigengenes will be combined in a linear model
    ## (common), then the em eigengene is redundant and not needed.
    
    ## Naming combiningMu
    names(combiningMu) <- paste0("em", combiningMu)
    names(combiningMu)[names(combiningMu)=="emNA"] <- "em" ##e.g., for when combiningMu=NA.
    eigenlociS <- NULL
    if(!is.null(eigenloci)){
        patients <- rownames(eigenloci)
        eigenlociS <- scale(eigenloci[patients, geNames])
        lCombine <- c(lCombine, "m"=1)
    }
    genExprS <- NULL
    if(!is.null(genExpr)){
        genExpr <- t(genExpr)
        patients <- rownames(genExpr)
        genExprS <- scale(genExpr[patients, geNames])
        lCombine <- c(lCombine, "e"=0)
    }
    ## If both eigenlociS and genExprS are provided,
    ##restrict every thing to the overlap. 
    if(!is.null(eigenlociS) & !is.null(genExprS)){
        patients <- intersect(rownames(eigenlociS), rownames(genExprS))
        eigenlociS <- eigenlociS[patients, geNames]
        genExprS <- genExprS[patients, geNames]
        lCombine <- c(lCombine, combiningMu)
        totalNum <- length(unique(union(rownames(eigenlociS), rownames(genExprS))))
        if(length(patients) < totalNum & doWarn)
            warning("Computing combined eigengene is possible for only ",
                    length(patients), " out of ", totalNum, 
                    "patients, as some cases do not have expr or dnam data!")
    }
    result[["patients"]] <- patients
    result[["lCombine"]] <- lCombine 
    patientLabel <- Labels[patients]
    ## QC
    if(length(patientLabel)==0){
        stop("No patients! Do names of Labels correspond to rownames of
             eigenLoci and genExpr?")
    }
    ## Main group of patients:
    mainGroup <- patientLabel[patientLabel %in% c(Label1, Label2)]

    ##DeadPatientsTime
    DeadPatientsTime <- survival[survival[, event]==1, "Time"]
    eigengeneFiles <- matrix("NA", nrow=length(mus), ncol=length(lCombine),
                             dimnames=list(as.character(mus), names(lCombine)))
    for(mu in mus){ ##mu values
        muChar <- as.character(mu)
        ## extract combinedNetwork, which determines the module assignment for each mu value
        combinedNetwork <- mu2modules[paste(mu), ]
        
        for(ind in seq_along(lCombine)){ ##e.g., lCombine=c("e"=0, "m"=1, "em0.6"=0.6)
            l1 <- lCombine[ind]
            eOrM <-  names(lCombine)[ind]
            if(is.na(l1))
                l1 <- mu
            ## Combining:
            message.if(verbose=verbose+1, paste("mu value:", mu, " lCombine:", l1))

            combinedValue <- combineValues(genExpr=genExprS, eigenloci=eigenlociS,
                                             l1=l1, verbose=verbose+1)
            ##QC
            if(any(names(combinedNetwork)!=colnames(combinedValue)))
                stop("Genes in the network and the data differ!")
            ## Out path:
            muFolder <- file.path(netPath, paste0("mu", mu), paste0("Pigen_", eOrM))
            dir.create(muFolder, recursive=TRUE, showWarnings=FALSE)
            pigengenesFile <- file.path(muFolder, "pigengene.RData")
            ##
            ## QC:
            if(!sum(combinedValue!=0, na.rm=TRUE))
                stop("combinedValue is all zero or NAs!! Cannot compute eigengenes!")
            Data1 <- combinedValue[names(mainGroup),]
            naGenes <- names(which(is.na(colSums(Data1))))##genes with >=1 NA value!
            if(length(naGenes) > 1){
                if(!doIgnoreNas){
                    stop("There are ", length(naGenes), " genes with NA values!")
                } else {
                    Data1 <- Data1[,!(colnames(Data1) %in% naGenes)]
                    if(doWarn)
                        warning(length(naGenes), " genes with NA values were removed.")
                }
            }
            pigengenes <- compute.pigengene(Data=Data1, Labels=mainGroup,
                                            modules=combinedNetwork[colnames(Data1)], 
                                            selectedModules="All", doPlot=TRUE, verbose=verbose-2,
                                            saveFile=pigengenesFile)
            ## Projection for all patients
            inferred <- project.eigen(Data=combinedValue, pigengene=pigengenes, 
                                      naTolerance=0.05, verbose=2, ignoreModules=c())
            eigengenes <- inferred$projected
            ##Adding e, m or em to the colnames
            colnames(eigengenes) <- paste0(colnames(eigengenes), "-", eOrM)
            ## Correlation of eigengenes with time:
            common <- intersect(names(DeadPatientsTime), rownames(eigengenes))
            eigenCor <- stats::cor(eigengenes[common, ], as.numeric(DeadPatientsTime[common]))
            eigengeneFiles[muChar, eOrM] <- file.path(muFolder, "eigengenes.RData")
            ## Save:
            save.if(eigengenes, file=eigengeneFiles[muChar, eOrM], verbose=verbose-2)
            save.if(eigenCor, file=file.path(muFolder, "eigenCor.RData"), verbose=verbose-3)
        }##^ for(l1 in lCombine)
        gc()
    }##^for(mu in mus)
    result[["eigengeneFiles"]] <- eigengeneFiles
    return(result)
}
