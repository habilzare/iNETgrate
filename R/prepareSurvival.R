prepareSurvival <- function(clinical, patientIDCol="bcr_patient_barcode",
                             eventCol="vital_status", event="Dead",
                             timeCol="days_to_last_followup",
                             riskFactorCol, riskCatCol=NULL, riskHigh=NULL,
                             riskLow=NULL, otherLabel=NULL, verbose=0){
    ## clinical: a matrix where row names are cases
    ## riskCatCol: The column in the clinical that has the risk levels.
    ## riskFactorCol: The column in the clinical that has description of the risk (e.g., cytogenetic
    ## abnormalities in AML).
    ## otherLabel: a vector with names the same as row names of clinical.
   
    ## QC
    if(all(c(is.null(riskCatCol), is.null(riskHigh), is.null(riskLow))))
        stop("All riskCatCol, riskHigh, and riskLow cannot be NULL!")
    if(is.null(riskCatCol))
        if(any(c(is.null(riskHigh), is.null(riskLow))))
            stop("riskHigh and/or riskLow cannot be NULL, when riskCatCol is NULL")
    if(is.null(riskHigh) != is.null(riskLow))
        stop("riskHigh and riskLow are not defined appropriately!")

    ## Removing duplicate patient clinical data  
    clinUnique <- clinical[!duplicated(clinical[, patientIDCol]), ]
    rownames(clinUnique) <- as.character(clinUnique[, patientIDCol])

    ## Preparing survival data
    ## If riskCatCol is available        
    if(!is.null(riskCatCol)){
        survival <- clinUnique[,c(patientIDCol, eventCol, timeCol,
                                  riskFactorCol, riskCatCol)]
        survival[, riskCatCol] <- as.character(survival[, riskCatCol])
        riskCategories <- unique(survival[, riskCatCol])
        if(!(any(c(is.null(riskHigh), is.null(riskLow))))){
            if(!(riskHigh %in% riskCategories && riskLow %in% riskCategories)){
                stop("riskHigh and riskLow values don't match riskCategories")
            }
            message.if(paste("Risk categories will be reassigned as High, Int and Low based on",
                            riskHigh, " and ", riskLow, "\n"), verbose=verbose-1)
            survival[which(survival[, riskCatCol] %in% riskHigh), riskCatCol] <- "High"
            survival[which(survival[, riskCatCol] %in% riskLow), riskCatCol] <- "Low"
            survival[which(!survival[, riskCatCol] %in% c("High", "Low")),
                     riskCatCol] <- "Int"
        }  ## End (if(!(any(c(is.null(riskHigh), is.null(riskLow)))))))
    } else{ ## If riskCatCol is NULL
        riskFactors <- unique(as.character(clinUnique[, riskFactorCol]))
        if(!(riskHigh %in% riskFactors && riskLow %in% riskFactors)){
            stop("riskHigh and riskLow values don't match riskFactors")
        } 
        survival <- clinUnique[,c(patientIDCol, eventCol, timeCol, riskFactorCol)]
        survival[which(as.character(survival[, riskFactorCol]) %in% riskHigh), 
                         "Risk1"] <- "High"
        survival[which(as.character(survival[, riskFactorCol]) %in% riskLow), 
                         "Risk1"] <- "Low"
        survival[which(!survival[, "Risk1"] %in% c("High", "Low")), 
                         "Risk1"] <- "Int"
    }  ## End else{}
    
    colnames(survival) <- c("PatientID", event, "Time", "Abnormality", "Risk1")

    ## Setting event col to 0 and 1
    survival[, event] <- as.factor(survival[, event])
    levels(survival[, event])[levels(survival[, event])!=event] <- 0
    levels(survival[, event])[levels(survival[, event])==event] <- 1
    
    DeadPatients <- rownames(survival[survival[, event]==1, ])

    if("days_to_death" %in% colnames(clinUnique)){
        survival[DeadPatients, "Time"] <- clinUnique[which(clinUnique[, patientIDCol] %in% DeadPatients), 
                                                     "days_to_death"]
    }

    ## Add more data:
    if(!is.null(otherLabel)){
        survival[, "Risk2"] <- otherLabel[rownames(survival)]
    }

    ## Checking for NA "Time"
    toKeep <- rownames(survival)[which(is.na(survival[, "Time"])==FALSE)]
    if(length(toKeep) < nrow(survival)){
        paste(nrow(survival) - length(toKeep), 
                 "cases are removed because of NA in the Time column. \n")
    }
    survival <- survival[toKeep,]

    ## Checking for NA event
    toKeep <- rownames(survival)[which(is.na(survival[,event])==FALSE)]
    if(length(toKeep) < nrow(survival)){
        paste(nrow(survival) - length(toKeep), 
                      "cases are removed because of NA in the event column. \n")
    }

    survival <- survival[toKeep,]
    survival[, event] <- as.numeric(as.character(survival[, event]))
    survival[, "Time"] <- as.numeric(survival[, "Time"])
    return(survival)
}
