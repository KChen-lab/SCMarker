\name{Bina}
\alias{Bina}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discretized matrix
}
\description{
Discretizing a matrix of continuous value to binary value 
}
\usage{
Bina(data, cutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A matrix which needs to be binary
}
  \item{cutoff}{
A numeric threshold. If the value in data is less than threshold, it is set as 0, otherwise 1.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Binary matrix}
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
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
Bina(melanoma1)
## The function is currently defined as
function (data, cutoff) 
{
    data[data < cutoff] = 0
    data[data >= cutoff] = 1
    return(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
