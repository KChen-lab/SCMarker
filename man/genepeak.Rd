\name{genepeak}
\alias{genepeak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating modal distribution
}
\description{
Estimating the number of peak corresponding to non-zero value from probability density distribution
}
\usage{
genepeak(x, width = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector
}
  \item{width}{
The numeric width which is used to estimate probability density distribution function
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
peak
}
\examples{
x=c(rnorm(50,4,0.5),rep(0,length=20))
genepeak(x)
genepeak(x,width=2)

## The function is currently defined as
function (x, width = 1) 
{
    if (length(x[x != 0]) > 1) {
        pp = density(x[x != 0], width = width)
        Hpeak = peak(pp)
        return(dim(Hpeak)[1])
    }
    else {
        return(0)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
