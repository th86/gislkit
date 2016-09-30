\name{getMI}
\title{Mutual information of two numeric vectors}
\usage{
getMI(x,y)
getMI(x,y, normalize=TRUE, negateMI=FALSE)
}
\arguments{
	\item{x}{A numeric vector} 
	\item{y}{A numeric vector}
	\item{bin}{An integer specifying the bin number}
	\item{so}{An integer specifying spline order}
	\item{normalize}{A logical indicating whether you want a normalized mutual information}
	\item{negateMI}{A logical indicating whether you want a signed mutual information}
}
\value{
	The mutual information of \code{x} and \code{y}
}

\description{
	Evaluates the mutual information of two aligned numeric vectors
}

\details{
	Evaluates the mututal information using B-spline functions.
}

\examples{
	getMI(c(0.1,0.5,0.9,10),c(0,0.2,0.7,0.9))
}

\references{
  Carsten O Daub, Ralf Steuer, Joachim Selbig and Sebastian Kloska,
  Estimating mutual information using B-spline functions - an improved similarity measure for 
  analysing gene expression data. BMC Bioinformatics 2004, 5:118.

  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
}

\author{Wei-Yi Cheng}