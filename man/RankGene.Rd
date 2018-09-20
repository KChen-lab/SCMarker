\name{RankGene}
\alias{RankGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gene pair order
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RankGene(kk, k, HamD, geneName, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kk}{
%%     ~~Describe \code{k} here~~
}
  \item{k}{
%%     ~~Describe \code{MNN} here~~
}
  \item{HamD}{
%%     ~~Describe \code{HamD} here~~
}
  \item{geneName}{
%%     ~~Describe \code{geneName} here~~
}
  \item{n}{
%%     ~~Describe \code{MNNIndex} here~~
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
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (kk, k, HamD, geneName, n)
{
    x = HamD[k, ]
    xrank = rank(-x)
    MNNgene = geneName[xrank <= k & x > n]
    return(MNNgene)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
