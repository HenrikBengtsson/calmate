%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fitCalMaTeCNprobes.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{fitCalMaTeCNprobes.matrix}
\alias{fitCalMaTeCNprobes.matrix}
\alias{fitCalMaTeCNprobes.matrix}

\alias{fitCalMaTeCNprobes}

 \title{Normalizes non-polymorphic copy number loci according to the CalMaTe method}

 \description{
  Normalizes non-polymorphic copy number loci according to the CalMaTe method.
 }

 \usage{\method{fitCalMaTeCNprobes}{matrix}(T, references, ...)}

 \arguments{
  \item{T}{A JxI \code{\link[base]{numeric}} \code{\link[base]{matrix}}, where J is the number of loci
                      and I is the number of samples.}
  \item{references}{A \code{\link[base]{logical}} or \code{\link[base]{numeric}} \code{\link[base]{vector}} specifying which
     samples should be used as the reference set.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J.
 }


\keyword{methods}
\keyword{internal}
