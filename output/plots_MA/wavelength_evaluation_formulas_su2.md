# Smoothed Wavelength Diagnostics

This note documents the exploratory wavelength estimates used for the smoothed
topological charge histories.  The quantity plotted and analyzed is the
Gaussian-smoothed history of each continuous run segment.  Estimates are made
per segment and then shown together in the summary plot.

## Gaussian Smoothing

For samples \(Q_i = Q(t_i)\), the smoothed curve at sample \(t_i\) is

\[
\widetilde Q(t_i)
=
\frac{\sum_j Q(t_j)
\exp\!\left[-\frac{(t_i-t_j)^2}{2\sigma^2}\right]}
{\sum_j
\exp\!\left[-\frac{(t_i-t_j)^2}{2\sigma^2}\right]} .
\]

For the thesis-style PBC/OBC comparison, the same smoothing width is applied to
PBC and OBC at fixed lattice spacing.  The display widths are chosen to increase
toward smaller lattice spacing:

\[
\sigma =
\begin{cases}
10, & a = 0.0299~\mathrm{fm},\\
20, & a = 0.0191~\mathrm{fm},\\
40, & a = 0.0107~\mathrm{fm}.
\end{cases}
\]

This is a visual and diagnostic convention, not a fitted physical parameter.
The separate stability scan below is used to check whether wavelength estimates
are stable under changes of \(\sigma\).

For the wavelength diagnostics the unrounded smoothed curve is used.  Rounding
is useful visually, but it can flatten extrema and bias peak finding.

## Artifact Checks

Gaussian smoothing can create visually wave-like structure even in stochastic
histories.  The extracted \(\lambda\) values are therefore treated as diagnostic
time scales, not as proof of a physical oscillation.

The first and last \(3\sigma\) samples define the conservative usable region:

\[
i \in [3\sigma, N-3\sigma].
\]

The usable length is

\[
N_{\mathrm{usable}} = N - 6\sigma.
\]

FFT estimates are computed on this usable region.  The peak/extrema diagnostic
is allowed to use the full broad-smoothed segment, because otherwise visible
troughs close to the boundary can be removed.  Such estimates are marked with
`edge_sensitive = yes` in the CSV.

An estimate is marked as robust only if the wavelength is larger than the
smoothing width,

\[
\lambda > 3\sigma,
\]

and if the usable region contains at least two wavelengths,

\[
N_{\mathrm{usable}} \geq 2\lambda.
\]

The first condition rejects scales directly tied to the smoothing kernel.  It is
not sufficient by itself, because smoothed random histories can also generate
pseudo-periodic extrema.  The second condition rejects estimates where the
segment is too short to contain repeated cycles.

The diagnostic CSV files use the following status labels and flags:

- `ok`: the estimate passes the edge, smoothing-scale, and usable-cycle checks.
- `edge_too_short`: the edge mask leaves too little usable data.
- `lambda_too_close_to_sigma`: the estimated wavelength is too close to the
  smoothing width.
- `too_few_cycles`: the usable segment does not contain at least two
  wavelengths.
- `no_estimate`: the method did not return a usable estimate.
- `edge_sensitive`: a separate `yes`/`no` flag.  If it is `yes`, at least one
  peak/trough used by the peak diagnostic lies inside the \(3\sigma\) boundary
  region.

## Detrending

Before estimating periods, the smoothed curve is linearly detrended:

\[
y_i = \widetilde Q(t_i) - (\alpha t_i + \beta),
\]

where \(\alpha\) and \(\beta\) are obtained from a least-squares linear fit.
The mean is then subtracted,

\[
\bar y_i = y_i - \frac{1}{N}\sum_{k=1}^{N} y_k .
\]

This reduces the chance that a slow drift is mistaken for a periodic signal.

## Peak Spacing

For the peak diagnostic, the already smoothed curve is smoothed once more with
the same \(\sigma\).  This gives a broader extrema curve and prevents small
local staircase wiggles from being counted as physical peak-to-peak structure.
Local maxima of \(\bar y_i\) and local minima of \(-\bar y_i\) are then found
with a minimum distance of approximately \(2\sigma\) samples and a prominence
threshold of \(0.18\) times the range of the broad detrended curve.

If repeated maxima or repeated minima are available, the peak-spacing estimate
uses the median of the distances between neighboring extrema of the same type:

\[
\lambda_{\mathrm{peaks}}
=
\mathrm{median}
\left(
\{p_{m+1}-p_m\}
\cup
\{q_{m+1}-q_m\}
\right),
\]

where \(p_m\) are peak positions and \(q_m\) are trough positions.

If a segment contains only one visible maximum and one visible minimum, the
method reports the corresponding peak-to-trough diagnostic scale,

\[
\lambda_{\mathrm{peaks}}
=
2\,\mathrm{median}\left(\{|e_{m+1}-e_m|\}\right),
\]

where \(e_m\) and \(e_{m+1}\) are neighboring extrema of opposite type.  Such
values are useful for checking a visible wave-like feature, but they can still
receive the `too_few_cycles` status if the segment is too short to contain two
full wavelengths.

## Peak-To-Trough Amplitude

For every neighboring maximum/minimum pair used by the peak diagnostic, the
peak-to-trough height is

\[
\Delta Q_m = |e_{m+1}-e_m| .
\]

The reported amplitude is the average over these visible peak-to-trough pairs,

\[
\langle \Delta Q\rangle_{\mathrm{peak-trough}}
=
\frac{1}{M}\sum_{m=1}^{M}\Delta Q_m .
\]

The file `amplitudes_su2.png` shows this quantity per segment.  It is a
diagnostic size of the visible wave-like structure, not a fitted physical
amplitude.

## Fourier Spectrum

A Hann window is applied before taking the real Fourier transform:

\[
Y_k
=
\mathcal{F}
\left[
w_i \bar y_i
\right]_k,
\qquad
w_i =
\frac{1}{2}
\left(
1-\cos\!\left(\frac{2\pi i}{N-1}\right)
\right).
\]

The power spectrum is

\[
P(f_k) = |Y_k|^2.
\]

The dominant period is

\[
\lambda_{\mathrm{FFT}}
=
\frac{1}{f_{k^\star}},
\qquad
k^\star =
\operatorname*{arg\,max}_{k} P(f_k),
\]

using only frequencies whose periods satisfy

\[
2\sigma \leq \frac{1}{f_k} \leq N.
\]

The file `fft_spectra_su2.png` shows the normalized Fourier power,
\(P_{\mathrm{FFT}}/P_{\mathrm{max}}\), as a function of
\(\lambda_{\mathrm{FFT}} = 1/f\).  Each continuous run segment is normalized to
unit maximum power, so the plot compares the location and shape of the spectral
peaks rather than the absolute power normalization.

## Sigma Stability Scan

The file `wavelength_sigma_scan_su2.png` shows
\(\lambda(\sigma)\) for all non-\(160\times80^3\) ensembles.  The accompanying
tables are:

- `wavelength_sigma_scan_detailed_su2.csv`: one row per segment, method, and
  \(\sigma\).
- `wavelength_sigma_scan_medians_su2.csv`: segment medians at fixed
  \(\sigma\).
- `wavelength_sigma_suggestions_su2.csv`: a suggested \(\sigma\) per ensemble.

The suggested value is not chosen by making PBC and OBC agree.  It is chosen as
the smallest \(\sigma\) in a stable region.  For every ensemble, several values
of \(\sigma\) are scanned.  At each \(\sigma\), the history is smoothed and the
wavelength \(\lambda(\sigma)\) is estimated with the peak-spacing and Fourier
methods.  A method contributes to the scan only if its estimate passes the
artifact checks above.

The stability score uses the symmetric relative change between neighboring
scanned smoothing widths:

\[
\delta_\sigma \lambda_m
=
\frac{
2|\lambda_m(\sigma_{r+1})-\lambda_m(\sigma_r)|
}{
\lambda_m(\sigma_{r+1})+\lambda_m(\sigma_r)
}.
\]

A candidate is considered stable if the average neighboring change satisfies

\[
\delta_\sigma \lambda_m \lesssim 0.05
\]

for at least two valid methods.  Thus the suggested value means:

\[
\sigma_{\mathrm{suggested}}
=
\min\left\{
\sigma :
\delta_\sigma \lambda_m \lesssim 0.05
\ \mathrm{for\ at\ least\ two\ methods}
\right\}.
\]

In words, it is the first smoothing width where the wavelength estimate has
become approximately stable.  Small \(\sigma\) values can leave too many local
wiggles, while very large \(\sigma\) values can oversmooth the history.  If
valid estimates exist but no stable region is found, the ensemble is marked as
`unstable_sigma` and no suggested \(\sigma\) is reported.  If no valid estimates
exist, it is marked as `no_stable_candidate`.  These suggestions are therefore
diagnostic choices, not physical fit parameters.  In the main PBC/OBC
comparison, common display values of \(\sigma\) are used at fixed lattice
spacing so that PBC and OBC are filtered in the same way.

## Reading The Summary Plot

In `wavelengths_su2.png`:

- `o` markers are PBC estimates.
- `x` markers are OBC estimates.
- Blue points are PBC estimates.
- Orange points are OBC estimates.
- Each point is one continuous run segment.
- Each panel is one lattice spacing.  The PBC and OBC estimates are shown in
  the same panel with a small horizontal offset.
- The upper row lists the wavelength-estimation method:
  `peak` and `FFT`.
- The plotted summary shows every finite wavelength returned by the methods.
  The reliability interpretation is stored in the CSV status column.  In
  particular, values marked `too_few_cycles` are diagnostic scales only: the
  segment is too short to confirm that the scale repeats.
- Missing estimates are left blank; lattice spacings with no plotted estimates
  are omitted from the summary plot.

In `amplitudes_su2.png`, the same PBC/OBC marker convention is used for the
peak-to-trough amplitude diagnostic.
