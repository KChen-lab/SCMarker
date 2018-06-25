\name{SCMarker-package}
\Rdversion{1.0}
\alias{SCMarker-package}
\alias{SCMarker}
\docType{package}
\title{
\packageTitle{SCMarker}
}
\description{
\packageDescription{SCMarker}
}
\details{
\tabular{ll}{
Package: \tab SCMarker\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2018-06-08\cr
License: \tab GPL(>=2)\cr
}
}
\author{
Fang Wang

Maintainer: Fang Wang <fwang9@mdanderson.org>

}
\references{

}

\keyword{ SCMarker }
\seealso{

}
\examples{
~data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
filterres=ModalFilter(data=melanoma1,geneK=10,cellK=10,width=2,cutoff=2)
filterres=GeneFilter(filterres=filterres,maxexp=dim(filterres$data)[1]*0.8)
filterres=getMarker(filterres=filterdata,MNN=200,MNNIndex=20)
}
