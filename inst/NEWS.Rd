\name{NEWS}
\title{News for \R Package \pkg{iNETgrate}}

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
