\name{myeloma}
\alias{myeloma}
\title{Multiple Myeloma Data}
\description{
Krall, Uthoff, and Harley (1975) analyzed data from a study
on multiple myeloma in which researchers treated 65 patients with
alkylating agents.  Of those patients, 48 died during the study
and 17 survived.  In the data set MYELOMA, the variable TIME
represents the survival time in months from diagnosis. The variable
VSTATUS consists of two values, 0 and 1, indicating whether the
patient was alive or dead, respectively, at the of end the study.
If the value of VSTATUS is 0, the corresponding value of TIME is
censored.  The variables thought to be related to survival are
LOGBUN (log BUN at diagnosis), HGB (hemoglobin at diagnosis),
PLATELET (platelets at diagnosis: 0=abnormal, 1=normal), 
AGE (age at diagnosis in years), LOGWBC (log WBC at diagnosis), FRAC
(fractures at diagnosis: 0=none, 1=present), LOGPBM (log percentage
of plasma cells in bone marrow), PROTEIN (proteinuria at 
diagnosis), and SCALC (serum calcium at diagnosis).

Data are from 
\code{http://ftp.sas.com/techsup/download/sample/samp_lib/}

\code{statsampExamples_of_Coxs_Model.html}
}
\usage{data(myeloma)}
\format{A data frame containing 65 observations on 11 variables:
\tabular{rl}{
        [,1] \tab "time"\cr
        [,2] \tab "vstatus"\cr
        [,3] \tab "logBUN"\cr
        [,4] \tab "HGB"\cr
        [,5] \tab "platelet"\cr
        [,6] \tab "age"\cr
        [,7] \tab "logWBC"\cr
        [,8] \tab "FRAC"\cr
        [,9] \tab "logPBM"\cr
        [,10] \tab "protein"\cr
        [,11] \tab "SCALC"}

}
\references{
Krall, J.M., Uthoff, V.A., and Harley, J. B. (1975), 
A Step-up Procedure for Selecting Variables Associated with Survival, 
{\em Biometrics}, {\bf 31}, 49 -57.
}
\keyword{datasets}
