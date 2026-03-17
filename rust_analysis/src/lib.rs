/// qcd_analysis — fast analysis routines for lattice QCD topological charge.
///
/// Exposed Python functions (all accept numpy arrays as Vec<f64>):
///
///   find_optimal_alpha(q, alpha_min, alpha_max, tol) -> float
///     Parallel golden-section search for α minimising ⟨(αQ − round(αQ))²⟩.
///
///   jackknife_error(data) -> float
///     O(N) delete-one jackknife error.  Much faster than the O(N²) loop.
///
///   bootstrap_alpha(q, n_boot, alpha_min, alpha_max, tol) -> (float, float)
///     Parallel bootstrap → (alpha_central, alpha_std).
///
///   bootstrap_Q2(q, n_boot, alpha_min, alpha_max, tol) -> (float, float)
///     Parallel bootstrap on ⟨Q_re²⟩ → (Q2_central, Q2_std).
///
/// Build (from rust_analysis/):
///   maturin develop --release

use pyo3::prelude::*;
use rayon::prelude::*;
use rand::prelude::*;
use rand::rngs::SmallRng;

// ---------------------------------------------------------------------------
// Core maths (plain Rust, no Python types)
// ---------------------------------------------------------------------------

/// ⟨(αQ − round(αQ))²⟩ evaluated in parallel over the full Q array.
#[inline]
fn alpha_objective(q: &[f64], alpha: f64) -> f64 {
    let n = q.len();
    if n == 0 {
        return f64::INFINITY;
    }
    let sum: f64 = q
        .par_iter()
        .map(|&qi| {
            let d = alpha * qi - (alpha * qi).round();
            d * d
        })
        .sum();
    sum / n as f64
}

/// Golden-section minimisation of alpha_objective on [alpha_min, alpha_max].
fn golden_section(q: &[f64], alpha_min: f64, alpha_max: f64, tol: f64) -> f64 {
    let phi = (5.0_f64.sqrt() - 1.0) / 2.0;
    let (mut a, mut b) = (alpha_min, alpha_max);
    let (mut c, mut d) = (b - phi * (b - a), a + phi * (b - a));

    while (b - a).abs() > tol {
        if alpha_objective(q, c) < alpha_objective(q, d) {
            b = d;
            d = c;
            c = b - phi * (b - a);
        } else {
            a = c;
            c = d;
            d = a + phi * (b - a);
        }
    }
    (a + b) / 2.0
}

// ---------------------------------------------------------------------------
// Python bindings
// Note: Vec<f64> as parameter type lets pyo3 accept numpy arrays directly
// via their buffer protocol — no separate numpy crate needed.
// ---------------------------------------------------------------------------

/// Find α minimising ⟨(αQ − round(αQ))²⟩ via golden-section search.
///
/// Parameters
/// ----------
/// q         : numpy float64 array of raw topological charge values
/// alpha_min : lower search bound (default 0.8)
/// alpha_max : upper search bound (default 1.2)
/// tol       : convergence tolerance (default 1e-6)
#[pyfunction]
#[pyo3(signature = (q, alpha_min=0.8, alpha_max=1.2, tol=1e-6))]
fn find_optimal_alpha(q: Vec<f64>, alpha_min: f64, alpha_max: f64, tol: f64) -> f64 {
    golden_section(&q, alpha_min, alpha_max, tol)
}

/// O(N) jackknife standard error using the delete-one mean identity.
///
/// Replaces the O(N²) Python loop; identical result.
///
/// Parameters
/// ----------
/// data : numpy float64 array
#[pyfunction]
fn jackknife_error(data: Vec<f64>) -> f64 {
    let n = data.len();
    if n < 2 {
        return 0.0;
    }
    let nf = n as f64;
    let total: f64 = data.iter().sum();
    // θ̂_i = (total − x_i) / (n−1);  mean(θ̂) = total/n
    let jack_mean = total / nf;
    let var_jack: f64 = data
        .par_iter()
        .map(|&xi| {
            let ji = (total - xi) / (nf - 1.0);
            (ji - jack_mean) * (ji - jack_mean)
        })
        .sum();
    ((nf - 1.0) / nf * var_jack).sqrt()
}

/// Parallel bootstrap standard error on α.
///
/// Parameters
/// ----------
/// q         : numpy float64 array of raw Q values
/// n_boot    : bootstrap samples (default 500)
/// alpha_min, alpha_max, tol : search bounds (defaults 0.8, 1.2, 1e-6)
///
/// Returns (alpha_central, alpha_std).
#[pyfunction]
#[pyo3(signature = (q, n_boot=500, alpha_min=0.8, alpha_max=1.2, tol=1e-6))]
fn bootstrap_alpha(
    q: Vec<f64>,
    n_boot: usize,
    alpha_min: f64,
    alpha_max: f64,
    tol: f64,
) -> (f64, f64) {
    let n = q.len();
    let alpha_central = golden_section(&q, alpha_min, alpha_max, tol);
    if n_boot < 2 {
        return (alpha_central, 0.0);
    }

    let alphas: Vec<f64> = (0..n_boot)
        .into_par_iter()
        .map(|seed| {
            let mut rng = SmallRng::seed_from_u64(seed as u64);
            let sample: Vec<f64> = (0..n).map(|_| q[rng.gen_range(0..n)]).collect();
            golden_section(&sample, alpha_min, alpha_max, tol)
        })
        .collect();

    let mean = alphas.iter().sum::<f64>() / n_boot as f64;
    let var = alphas.iter().map(|&a| (a - mean) * (a - mean)).sum::<f64>()
        / (n_boot - 1) as f64;
    (alpha_central, var.sqrt())
}

/// Parallel bootstrap standard error on ⟨Q_re²⟩.
///
/// Computes α on the full sample, builds Q_re = round(α·Q), then bootstraps
/// the mean of Q_re².
///
/// Parameters
/// ----------
/// q      : numpy float64 array of raw Q values
/// n_boot : bootstrap samples (default 500)
///
/// Returns (Q2_central, Q2_std).
#[pyfunction]
#[pyo3(name = "bootstrap_Q2", signature = (q, n_boot=500, alpha_min=0.8, alpha_max=1.2, tol=1e-6))]
fn bootstrap_q2(
    q: Vec<f64>,
    n_boot: usize,
    alpha_min: f64,
    alpha_max: f64,
    tol: f64,
) -> (f64, f64) {
    let n = q.len();
    let alpha = golden_section(&q, alpha_min, alpha_max, tol);

    let q_re: Vec<f64> = q.iter().map(|&qi| (alpha * qi).round()).collect();
    let q2_central = q_re.iter().map(|&qi| qi * qi).sum::<f64>() / n as f64;

    if n_boot < 2 {
        return (q2_central, 0.0);
    }

    let q2_boots: Vec<f64> = (0..n_boot)
        .into_par_iter()
        .map(|seed| {
            let mut rng = SmallRng::seed_from_u64(seed as u64 + 99991);
            (0..n)
                .map(|_| {
                    let qi = q_re[rng.gen_range(0..n)];
                    qi * qi
                })
                .sum::<f64>()
                / n as f64
        })
        .collect();

    let mean = q2_boots.iter().sum::<f64>() / n_boot as f64;
    let var = q2_boots.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>()
        / (n_boot - 1) as f64;
    (q2_central, var.sqrt())
}

// ---------------------------------------------------------------------------
// Module
// ---------------------------------------------------------------------------

#[pymodule]
fn qcd_analysis(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_optimal_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(jackknife_error, m)?)?;
    m.add_function(wrap_pyfunction!(bootstrap_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(bootstrap_q2, m)?)?;
    Ok(())
}
