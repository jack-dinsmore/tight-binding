use super::complex::*;

const FINE_STRUCTURE: f64 = 0.0072973525;

pub fn ylm(l: u32, m: u32) -> impl Fn((f64, f64)) -> Complex {
    unimplemented!()
}

pub fn legendre(l: u32, m: u32) -> impl Fn((f64, f64)) -> f64 {
    unimplemented!()
}

pub fn nlm(n: u32, l: u32, m: u32) -> impl Fn((f64, f64, f64)) -> Complex {
    #[cfg(debug)]
    {
        assert_eq!(m <= l && -m <= l, true);
    }

    // Position is in units of a_0

    let norm = ((2.0 / n as f64).pow(3) * (n - l - 1).factorial() as f64 / (2 * n * (n + l).factorial()) as f64).sqrt();
    let rho_scale = 2.0 / n as f64;
    let L = legendre(n - l - 1, 2 * l + 1);
    let Y = ylm(l, m);

    move |x: (f64, f64, f64) | -> f64 {
        let r = (x.0 * x.0 + x.1 * x.1 + x.2 * x.2).sqrt();
        let r_asym = (x.0 * x.0 + x.1 * x.1).sqrt();
        let cos_theta = x.2 / r;
        let cos_phi = x.0 / r_asym;
        let sin_phi = x.1 / r_asym;

        let rho = r * rho_scale;

        norm * rho.pow(l) * (-rho / 2.0).exp() * L(rho) * Y(cos_theta, cos_phi, sin_phi)
    }
}

pub fn energy(n: u32) -> f64 {
    // Units where position is a_0, mass is m_e, and t such that c = 1.
    - FINE_STRUCTURE * FINE_STRUCTURE / 2.0 / n as f64 / n as f64
}