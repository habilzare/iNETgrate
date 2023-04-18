\name{prepareSurvival}
\alias{prepareSurvival}

\title{Prepare survival data}
\description{
    Processes and cleans clinical data to gather the required survival information.}

\usage{
prepareSurvival(clinical, patientIDCol="bcr_patient_barcode", eventCol="vital_status", 
                 event="Dead", timeCol="days_to_last_followup", riskFactorCol, 
                 riskCatCol=NULL, riskHigh=NULL, riskLow=NULL, otherLabel=NULL,
                 verbose=0)}

\arguments{
  \item{clinical}{
    A matrix or a data frame of clinical data with a minimum four columns that
    provide information about patient IDs, event occurance, survival or follow up time, 
    and risk factors. The rows are patients.}

  \item{patientIDCol}{
    A string determining the column in clinical data that lists patient IDs
    or case IDs, e.g., "bcr_patient_barcode" in the sample AML data.}

  \item{eventCol}{
    A string determining the column in clinical data that lists the
    \code{event}. The Event column must be a binary column
    E.g. "vital_status" column in sample AML data  where two possible
    events are "Dead" and "Alive".}

  \item{event}{
    A string determining the event of interest from \code{eventCol} in clinical
    data e.g., "Dead" is the event of interest in the sample data.}

  \item{timeCol}{
    A string determining the column in clinical data storing survival time
    information e.g., "days_to_last_followup" in sample AML data}

  \item{riskFactorCol}{
    A string determining the column in clinical data that describes the risk factors
    associated with risk levels.}

  \item{riskCatCol}{
    A string determining the column in clinical data that lists risk levels of 
    the disease.}

  \item{riskHigh}{
    A character vector that describes high-risk group of patients in clinical data. 
    See section "Details".}

  \item{riskLow}{
    A character vector that describes low-risk group of patients in clinical data. 
    See section "Details".}

  \item{otherLabel}{
    A vector of other risk levels where names are the same as row names of clinical 
    data.}

  \item{verbose}{
    An integer that sets the level of verbosity. 0 means silent, and higher values produce 
    more details of the computation.}
}

\details{
    This function removes any cases in clinical data for which vital status or 
    survival time is missing. When \code{riskCatCol} is \code{NULL},
    \code{riskHigh} and \code{riskLow} are expected to have values that belong
    to \code{riskFactorCol}. When \code{riskCatCol} is defined, \code{riskHigh} 
    and \code{riskLow} are expected to have values that belong to
    \code{riskCatCol}. All \code{riskCatCol},\code{riskHigh} and \code{riskLow} 
    cannot be \code{NULL}.}

\value{
    A data frame with a minimum of five columns namely: "PatientID", \code{event}, 
    "Time", "Risk1", and "Abnormality". The row names are the same as "PatientID".}

\seealso{
    \code{\link{downloadData}}, \code{\link{cleanAllData}},
    \code{\link{analyzeSurvival}}
}

\examples{
data(toyRawAml)

riskCatCol = "acute_myeloid_leukemia_calgb_cytogenetics_risk_category"
riskFactorCol = "cytogenetic_abnormalities"
survival <- prepareSurvival(clinical=toyRawAml$clinical,
                             patientIDCol="bcr_patient_barcode", eventCol="vital_status", 
                             event="Dead", timeCol="days_to_last_followup", 
                             riskCatCol=riskCatCol, riskFactorCol=riskFactorCol, 
                             riskHigh="Favorable", riskLow="Poor",
                             otherLabel=NULL)
}

\keyword{manip}
