\name{findDFC}
\alias{findDFC}
\title{
Finds the double first cousins in a pedigree}
\description{
Given a pedigree, all pairs of individuals that are double first cousins are returned.}

\usage{
findDFC(pedigree, exact = FALSE, parallel = FALSE,
	ncores = getOption("mc.cores", 2L))
}

\arguments{
  \item{pedigree }{
A pedigree with columns organized: ID, Dam, Sire
}
  \item{exact }{
A logical statement indicating if individuals who are exactly double first cousins are to be identified
}
  \item{parallel }{
A logical statement indicating if parallelization should be attempted.  Note, only reliable for Mac and Linux operating systems.
}
  \item{ncores }{
Number of cpus to use, default is maximum available
}
}

\details{
When exact = TRUE, only those individuals whose grandparents are completely unrelated will be identified as double first cousins.  When exact = FALSE, as long as the parents of individuals i and j are two sets of siblings (i.e., either sires full brothers/dams full sisters or two pairs of opposite sex full sibs) then i and j will be considered double first cousins.  In the event where the grandparents of i and j are also related, exact = FALSE will still consider i and j full sibs, even though genetically they will be more related than exact = TRUE double first cousins.

\code{parallel} = TRUE should only be used on Linux or Mac OSes (i.e., not Windows).
}

\value{
  \item{PedPositionList }{gives the list of row numbers for all the pairs of indidivuals that are related as double first cousins
}
  \item{DFC }{gives the list of IDs, as characters, for all the pairs of individuals that are related as double first cousins
}
  \item{FamilyCnt }{If two individuals, i and j, are double first cousins, then i's siblings will also be double first cousins with j's siblings.  Therefore, this is the total number of family pairs where offspring are related as double first cousins.
}
}
\author{\email{matthewwolak@gmail.com}
}

