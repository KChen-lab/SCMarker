\name{peak}
\alias{peak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of bi/multimodal distribution
}
\description{
Estimating peak number of probabiliy density distribution
}
\usage{
peak(pp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pp}{
The output of density function
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
The value corresponding to each peak
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
genepeak, ModalFilter
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
x=c(rnorm(50,4,0.5),rep(0,length=20))
pp=density(x)
peak(pp)
## The function is currently defined as
function (pp) 
{
    y = pp$y
    x = pp$x
    Hpeak = c()
    for (i in 2:(length(y) - 1)) {
        if (y[i] >= y[i - 1] & y[i] > y[i + 1]) {
            Hpeak = rbind(Hpeak, c(x[i], y[i]))
        }
    }
    if (length(Hpeak) == 0) {
        Hpeak = rbind(Hpeak, c(0, 0))
    }
    Hpeak = as.data.frame(Hpeak)
    names(Hpeak) = c("H", "Density")
    return(Hpeak)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
