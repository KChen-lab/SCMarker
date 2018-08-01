\name{MNNpair}
\alias{MNNpair}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene pairs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MNNpair(k, MNNgene, geneName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{MNNgene}{
%%     ~~Describe \code{MNNgene} here~~
}
  \item{geneName}{
%%     ~~Describe \code{geneName} here~~
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
function (k, MNNgene, geneName) 
{
    subgene = MNNgene[[k]]
    index = match(subgene, geneName)
    PP <- function(i, index, MNNgene, k, geneName) {
        if (geneName[k] \%in\% MNNgene[[index[i]]]) 
            return(geneName[index[i]])
    }
    if (length(index[!is.na(index)]) > 0) {
        a = do.call(cbind, lapply(1:length(index), PP, index = index, 
            MNNgene = MNNgene, k = k, geneName = geneName))
        if (!is.null(a)) {
            return(a)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
