mod hydrogen;
mod coeffs;
mod complex;

pub use coeffs::ProtonLayout;
use ndarray::Array2;
use ndarray_linalg::*;

pub fn hamiltonian<const PROTONS: usize, const N_MAX: usize>(layout: &ProtonLayout<PROTONS>) -> Array2<f64> {
    let qs = coeffs::make_q::<PROTONS, N_MAX>(layout);
    let ss = coeffs::make_s::<PROTONS, N_MAX>(layout);
    let dim = PROTONS * N_MAX * (N_MAX + 1) * (2 * N_MAX + 1) / 6;
    let mut array = Array2::zeros((dim, dim));

    let mut i = 0;
    for ap in 0..PROTONS {
        for np in 0..N_MAX {
            for lp in 0..np {
                for mp in -(lp as isize)..(lp+1) as isize {
                    let mut j = 0;
                    for a in 0..PROTONS {
                        for n in 0..N_MAX {
                            for l in 0..n {
                                for m in -(l as isize)..(l+1) as isize {
                                    let mut val = 0.0;
                                    // Calculate Hamiltonian
                                    array[[i,j]] = val;
                                    j += 1;
                                }
                            }
                        }
                    }
                    i += 1;
                }
            }
        }
    }
    array
}

pub fn energies<const PROTONS: usize, const N_MAX: usize>(layout: &ProtonLayout<PROTONS>) -> Vec<f64> {
    let h = hamiltonian::<PROTONS, N_MAX>(layout);
    let (eigs, _) = h.eig().unwrap();
    eigs.to_vec().iter().map(|c| {println!("{}", c.im()); c.re()}).collect::<Vec<_>>()
}