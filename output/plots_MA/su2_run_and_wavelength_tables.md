# SU(2) Run And Wavelength Tables

The first table lists the SU(2) ensembles used in the main MA plots, excluding the dedicated $160\times80^3$ comparison runs.  The second table lists those $160\times80^3$ runs separately.  The third table summarizes the smoothed wavelength and peak-to-trough amplitude diagnostics.

\begin{table}[htbp]
\centering
\scriptsize
\caption{SU(2) ensemble parameters used in the main analysis.  The $160\times80^3$ runs are excluded and listed separately.}
\label{tab:su2-main-ensembles}
\begin{tabular}{lccccccc}
\hline
BC & $\beta$ & $a~[\mathrm{fm}]$ & lattice & cut & $N_{\rm run}$ & $N_{\rm MC}$/run & seeds \\
\hline
PBC & 2.500 & 0.0774 & $16\times 16^3$ & 0 & 2 & 400 & 1902689425, 2018065041 \\
OBC & 2.500 & 0.0774 & $22\times 16^3$ & 3 & 2 & 400 & 101631635, 2133559185 \\
PBC & 2.600 & 0.0560 & $22\times 22^3$ & 0 & 2 & 400 & 1488952329, 1604860937 \\
OBC & 2.600 & 0.0560 & $30\times 22^3$ & 4 & 2 & 400 & 1720270345, 1835695113 \\
PBC & 2.700 & 0.0408 & $30\times 30^3$ & 0 & 2 & 400 & 1231330499, 1346372035 \\
OBC & 2.700 & 0.0408 & $42\times 30^3$ & 6 & 2 & 400 & 1624488965, 1739959301 \\
PBC & 2.800 & 0.0299 & $41\times 41^3$ & 0 & 2 & 400 & 144114837, 259336341 \\
OBC & 2.800 & 0.0299 & $57\times 41^3$ & 8 & 2 & 400 & 1882700213, 1998022069 \\
PBC & 2.950 & 0.0191 & $65\times 65^3$ & 0 & 2 & 400 & 560401647, 1348947599 \\
OBC & 2.950 & 0.0191 & $89\times 65^3$ & 12 & 2 & 400 & 469074389, 1483978919 \\
PBC & 3.150 & 0.0107 & $115\times 108^3$ & 0 & 1 & 50 & 275974907 \\
OBC & 3.150 & 0.0107 & $159\times 108^3$ & 22 & 2 & 50 & 999060221, 1208018241 \\
\hline
\end{tabular}
\end{table}

\begin{table}[htbp]
\centering
\scriptsize
\caption{Dedicated SU(2) $160\times80^3$ comparison ensembles.}
\label{tab:su2-t160-ensembles}
\begin{tabular}{lccccccc}
\hline
BC & $\beta$ & $a~[\mathrm{fm}]$ & lattice & cut & $N_{\rm run}$ & $N_{\rm MC}$/run & seeds \\
\hline
PBC & 3.150 & 0.0107 & $160\times 80^3$ & 40 & 1 & 50 & 785195877 \\
OBC & 3.150 & 0.0107 & $160\times 80^3$ & 40 & 1 & 50 & 825666065 \\
\hline
\end{tabular}
\end{table}

\begin{table}[htbp]
\centering
\scriptsize
\caption{Smoothed wavelength and peak-to-trough amplitude diagnostics.  The status columns label the reliability checks described in the text.}
\label{tab:su2-wavelength-diagnostics}
\begin{tabular}{lccccccccc}
\hline
BC & $a~[\mathrm{fm}]$ & $\sigma$ & seg. & $\lambda_{\rm peak}$ & $\lambda_{\rm FFT}$ & $\langle\Delta Q\rangle$ & edge & peak status & FFT status \\
\hline
PBC & 0.0299 & 10 & 2 & 182 & 170 & 4.42 & no & too\_few\_cycles & ok \\
PBC & 0.0299 & 10 & 1 & 138 & 170 & 2.16 & no & ok & ok \\
OBC & 0.0299 & 10 & 2 & 203 & 340 & 3.24 & no & too\_few\_cycles & too\_few\_cycles \\
OBC & 0.0299 & 10 & 1 & 160 & 340 & 2.31 & no & ok & too\_few\_cycles \\
PBC & 0.0191 & 20 & 2 & 209 & 140 & 0.69 & no & too\_few\_cycles & ok \\
PBC & 0.0191 & 20 & 1 & 168 & 140 & 1.60 & no & too\_few\_cycles & ok \\
OBC & 0.0191 & 20 & 2 & 191 & 280 & 1.90 & no & too\_few\_cycles & too\_few\_cycles \\
OBC & 0.0191 & 20 & 1 & 272 & 280 & 3.47 & yes & too\_few\_cycles & too\_few\_cycles \\
\hline
\end{tabular}
\end{table}

