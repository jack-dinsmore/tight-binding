mod hydrogen;
mod coeffs;
mod complex;

pub use coeffs::ProtonLayout;
use ndarray::Array2;

pub fn hamiltonian<const PROTONS: usize, const N_MAX: usize>(layout: &ProtonLayout<PROTONS>) -> Array2<f64> {
    let qs = coeffs::make_q(layout);
    let ss = coeffs::make_s(layout);
    
}

pub fn energies<const PROTONS: usize, const N_MAX: usize>(layout: &ProtonLayout<PROTONS>) -> Vec<f64> {
    let h = hamiltonian::<PROTONS, N_MAX>(layout);
    let (eigs, _) = h.eig().unwrap();
    eigs
}