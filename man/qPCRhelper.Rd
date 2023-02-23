\name{qPCRhelper}
\alias{qPCRhelper}
%-
\title{qPCRhelper
%%  Calculates gene expression from qPCR Ct values
}
\description{
Computes nomalized Ct values (delta Ct) from raw qPCR Ct values and conducts test of significance using t.test. Plots expression values based from log2(2^(-1*delta delta Ct)) across groups per gene of interest.
}
\usage{
qPCRhelper(data.dir = NULL, ref.gene = NULL, ref.group = NULL,
plot.ref.group = NULL, plot.nrow = 1,
plot.title = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.dir}{file path.
}
  \item{ref.gene}{string: value should be one of column names in input table. Used for delta Ct computation.
}
  \item{ref.group}{string: value should be one of 'Group' values in input table. Used for delta delta Ct computation.
}
  \item{plot.ref.group}{string: value should be one of 'Group' values in input table. Used to set reference in plotting.
}
  \item{plot.nrow}{numeric: optional. Number of rows for plotting n number of plots corresponding to n number of genes.
}
  \item{plot.title}{plot title: optional.
}
}
\value{
A dataframe with columns for normalized Ct values (dCt), and gene expression (log2RelExp).
}
\note{
If gene names start with a number, e.g. 18S, please precede the gene name with
'X' without space, e.g.:18S -> X18S
}
\author{Wilson Jr. Aala
%%
}
\examples{
## Create sample table with expected 'Sample', 'Group', and gene Ct columns
Sample <- c("C1", "C2", "T1", "T2") #required column
Group <- c("C", "C", "T", "T") #required column
# Gene Ct values, at least two columns: one reporter, one target gene
GAPDH <- c(18.1,18.2,18.1,18.2) #reporter, ref.gene
IL4 <- c(30.1,30.5,20.1,20.2) #target
a <- data.frame(Sample,Group,GAPDH,IL4) #export using write.table(a,...)

## Write the data frame to a file in a temporary directory
temp_file <- file.path(tempdir(), "a.txt")
write.table(a, file = temp_file, sep = "\t")

## Run qPCRhelper directly on the file
library(qPCRhelper)
b <- qPCRhelper(data.dir=temp_file,
                ref.gene="GAPDH",
                ref.group="C",
                plot.ref.group="C",
                plot.nrow=1,
                plot.title="My cool qPCR data")
}
