\name{genecount}
\alias{genecount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count the number of non-zero
}
\description{
Count the number of non-zero value for a vector
}
\usage{
genecount(k, data, index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
The kth row or column
}
  \item{data}{
A matrix
}
  \item{index}{
This must match one of "row" or "column" to extract the kth row or column vector
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
genecount(k=6,data=melanoma1,index="row")

## The function is currently defined as
function (k, data, index) 
{
    if (index == "row") {
        x = data[k, ]
    }
    else {
        x = data[, k]
    }
    return(sum(x != 0))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
