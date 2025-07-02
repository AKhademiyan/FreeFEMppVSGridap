\documentclass[11pt]{article}
\usepackage[a4paper,margin=2.5cm]{geometry}
\usepackage{hyperref}
\usepackage{booktabs}

\title{FEM vs Gridap: Benchmarking Stokes Flow in a Wavy Channel}
\author{Alireza Khademiyan}
\date{\today}

\begin{document}
\maketitle

This project benchmarks two finite element solvers—\textbf{FreeFEM} and \textbf{Julia/Gridap}—for solving the two‐dimensional Stokes equations in a sinusoidal (wavy‐walled) channel.  Performance is compared in terms of execution time and memory usage across a range of mesh sizes.

\section*{Repository Contents}

\\begin{tabular}{@{}ll@{}}
\\toprule
\\textbf{File} & \\textbf{Description} \\\\
\\midrule
\texttt{stokes.edp} & FreeFEM script for the Stokes problem \\\\
\texttt{stokes1.jl} & Julia/Gridap script (same geometry/physics) \\\\
\texttt{wavy.geo} & Gmsh geometry script defining the channel \\\\
\texttt{wavy.msh} & Pre‐generated mesh used by both solvers \\\\
\texttt{data.ipynb} & Jupyter notebook for benchmarking/plots \\\\
\texttt{performance\_comparison.png} & Combined timing\,+\,memory plot \\\\
\texttt{timing\_vs\_mesh.png} & Execution time vs.\ mesh size \\\\
\texttt{memory\_vs\_mesh.png} & Peak memory vs.\ mesh size \\\\
\texttt{first.png} & Legacy combined figure \\\\
\texttt{report.pdf} & Scientific report with full analysis \\\\
\texttt{.DS\_Store} & macOS system file (ignore) \\\\
\\bottomrule
\\end{tabular}

\section*{Methods}

\\begin{itemize}
  \\item Taylor–Hood elements (\(\\mathrm{P2}\\)–\(\\mathrm{P1}\\)) for two‐dimensional Stokes flow.
  \\item Sinusoidal upper and lower walls; pressure BCs at inlet/outlet, no‐slip on walls.
\\end{itemize}

\section*{Performance Comparison}

\\begin{itemize}
  \\item \\textbf{Execution time}: measured with \texttt{@benchmark} (Julia) and \texttt{clock()} (FreeFEM).
  \\item \\textbf{Memory usage}: captured via \\texttt{@benchmark.memory} (Julia) and \\texttt{storageused()} (FreeFEM).
\\end{itemize}

\\textbf{Key findings}
\\begin{itemize}
  \\item FreeFEM shows lower memory overhead and faster runtimes for small–medium meshes.
  \\item Gridap matches performance at large scale but uses more memory.
  \\item Gridap offers modern tooling (GPU readiness, composability).
\\end{itemize}

\section*{How to Run}

\\subsection*{FreeFEM}
\\begin{verbatim}
FreeFEM++ stokes.edp
\\end{verbatim}

\\subsection*{Julia}
Install required packages
\\begin{verbatim}
using Pkg
Pkg.add(\"Gridap\")
Pkg.add(\"GridapGmsh\")
Pkg.add(\"BenchmarkTools\")
\\end{verbatim}
Then run
\\begin{verbatim}
julia stokes1.jl
\\end{verbatim}

\\subsection*{Jupyter Notebook}
\\begin{verbatim}
jupyter notebook data.ipynb
\\end{verbatim}

\section*{Report}
See \\href{./report.pdf}{\\texttt{report.pdf}} for detailed methodology, plots, and discussion.

\section*{License}
MIT License.  Feel free to use, modify, and share.

\section*{Author}
Alireza Khademiyan\\\\
\\url{https://github.com/AKhademiyan}

\section*{Future Improvements}
\\begin{itemize}
  \\item Add \\texttt{.tex} source for the report.
  \\item Automate benchmarking over mesh sizes.
  \\item Add unit tests for solver output validation.
\\end{itemize}

\\end{document}
