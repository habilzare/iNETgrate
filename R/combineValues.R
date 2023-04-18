combineValues <- function(genExpr, eigenloci, l1, verbose=0){  
    if(is.null(genExpr)){
        combinedValue <- l1*eigenloci ## will add doCorMethod later
        message.if(verbose=verbose, "Using only dnam")
    } else { ##genExpr is not NULL,
        combinedValue <- (1-l1)*genExpr
        message.if(verbose=verbose, "Using expr")
        if(!is.null(eigenloci)){
            message.if(verbose=verbose, "Add dnam")
            combinedValue <- combinedValue + l1*eigenloci
        }
    }
    
    return(combinedValue)
}
