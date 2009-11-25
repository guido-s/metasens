\name{Crowther2003}
\docType{data}
\alias{Crowther2003}
\title{Aspirin after Myocardial Infarction}

\description{Meta-analysis on phenobarbital prior to preterm birth for
  preventing neonatal periventricular haemorrhage}

\usage{data(Crowther2003)}

\format{A data frame with the following columns:
  \describe{
  \item{study}{study label}
  \item{event.e}{number of periventricular haemorrhages in experimental group}
  \item{total.e}{number of observations in experimental group}
  \item{event.c}{number of periventricular haemorrhages in control group}
  \item{total.c}{number of observations in control group}
  }
}

\source{
  Crowther CA, Henderson-Smart DJ (2003),
  Phenobarbital prior to preterm birth for preventing neonatal
  periventricular haemorrhage.
  \emph{Cochrane Database of Systematic Reviews}, Issue 3. Art. No.:
  CD000164. DOI: 10.1002/14651858.CD000164
}

\examples{
data(Crowther2003)
metabin(event.e, total.e, event.c, total.c,
        data=Crowther2003,
        studlab=study)
}

\keyword{datasets}
