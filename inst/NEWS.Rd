\name{NEWS}
\title{News for \R Package \pkg{iNETgrate}}

\section{Changes in version 0.99.102 (2023-06-12)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item In \code{prepareSurvival} and \code{cleanAllData} function,
      the default values changes as follows: \code{riskLow="Low"} and
      \code{riskHigh="High"}, respectively. 
    }
  }
}

\section{Changes in version 0.99.100 (2023-06-09)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{Labels} argument is removed from \code{analyzeSurvival}.
      and the \code{riskCol} argument is added.
    }
  }
}

\section{Changes in version 0.99.96 (2023-06-08)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{abnormalityCol} argument is removed from \code{accelFailAnalysis}.
      \item The \code{otherLabel} argument is removed from
      \code{prepareSurvival} and \code{cleanAllData} functions.
    }
  }
}

\section{Changes in version 0.99.94 (2023-06-08)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{computeEigengenes}, \code{computeEigenloci}, and
      \code{sample2atient} functions were renamed to
      \code{computEigengenes}, \code{computEigenloci}, and
      \code{sample2pat}, respectively.
      \item The \code{makeNetwork} function now saves the result for
      each mu in a separate folder.
    }
  }
}

\section{Changes in version 0.99.74 (2023-05-24)}{
  \subsection{General}{
    \itemize{
      \item Data were compressed using "xz" to save 20\% space.
    }
  }
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{downloaData} function is not exported anymore
      because it is specific to the legacy TCGA data, which is no longer
      available.
    }
  }
}

\section{Changes in version 0.99.70 (2023-05-23)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Habil renamed \code{downloadData} to \code{downloaData}
      and updated the \code{data.category} argument of \code{GDCquery}. 
      \item Habil removed the legacy option from the \code{donwlowData}
      function because it is no longer supported by TCGAbiolinks.
    }
  }
}

\section{Changes in version 0.99.64 (2023-05-11)}{
  \subsection{General}{
    \itemize{
      \item Habil commented the installation command line in the
      vignette.
    }
  }
}

\section{Changes in version 0.99.50 (2023-04-14)}{
  \subsection{General}{
    \itemize{
      \item Habil changed function names to camelCase i.e, from
      \code{function.name} to \code{functionName}.
    }
  }
}

\section{Changes in version 0.99.46 (2023-04-11)}{
  \subsection{General}{
    \itemize{
      \item Habil changed \code{getwd()} to \code{tempdr()} in the docs.
      \item Habil changed most of occurrences of \code{cat()} to \code{paste()}.
      \item Habil changed most of occurrences of \code{print()} to
      \code{message()} or \code{Pigengene::message.if()}.
      \item Habil removed \code{paste()} inside all \code{stop()} to \code{warning()}.
      \item Habil replaced \code{dontrun} with \code{donttest}.
    }
  }
}

\section{Changes in version 0.99.30 (2023-04-06)}{
  \subsection{General}{
    \itemize{
      \item Tested on R 4.3.
      \item TCGAbiolinks (<= 2.24.3) because of issues with newer versions.
    }
  }
}

\section{Changes in version 0.5.42 (2021-10-31)}{
  \subsection{General}{
    \itemize{
      \item Being prepared for submission to Bioconductor.
      \item Created.
    }
  }
}

%% \subsection{General}
%% \subsection{New functions}{
%% \subsection{Changes in existing functions} %% functionality, signature, output
%% \subsection{Bug Fixes}
