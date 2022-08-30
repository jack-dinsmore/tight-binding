use super::complex::*;
use std::f64::consts::PI;

const FINE_STRUCTURE: f64 = 0.0072973525;

fn parity(l: i32) -> f64 {
    if l % 2 == 0 { 1.0 } else { -1.0 }
}

fn factorial(num: u32) -> u32 {
    match num {
        0  => 1,
        1  => 1,
        2  => 2,
        3  => 6,
        4  => 24,
        5  => 120,
        6.. => (1..num+1).product(),
    }
}

fn choose(a: u32, b: u32) -> u32 {
    factorial(a) / (factorial(a - b) * factorial(b))
}

fn ylm(l: u32, m: i32) -> impl Fn(f64, f64, f64) -> Complex {
    let norm = parity(m) * ((2 * l + 1) as f64 / (4.0 * PI) * factorial((l as i32 - m) as u32) as f64 / factorial((l as i32 + m) as u32) as f64).sqrt();
    let P = legendre(l, m);
    move |cos_theta: f64, cos_phi: f64, sin_phi: f64| {
        Complex::new(cos_phi, sin_phi) * norm * P(cos_theta)
    }
}

fn legendre(l: u32, m: i32) -> impl Fn(f64) -> f64 {
    const LEGENDRE_NORM: [f64; 25] = [
        1.0,

        0.5,
        1.0,
        -1.0,

        3.0/24.0,
        3.0/6.0,
        0.5,
        -3.0,
        3.0,

        15.0/720.0,
        15.0/120.0,
        1.5/12.0,
        0.5,
        1.5,
        15.0,
        -15.0,

        105.0/40320.0,
        105.0/5040.0,
        15.0/2.0 * 1.0/360.0,
        2.5/20.0,
        1.0/8.0,
        -2.5,
        15.0/2.0,
        -105.0,
        105.0,
    ];
    const LEGENDRE_COEFF: [[f64; 4]; 15] = [
        [0.0, 0.0, 0.0, 0.0],

        [1.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],

        [-1.0, 3.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],

        [-3.0, 5.0, 0.0, 0.0],
        [1.0, -5.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],

        [3.0, -30.0, 35.0, 0.0],
        [-3.0, 7.0, 0.0, 0.0],
        [-1.0, 7.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],

    ];
    let norm_index = ((l * (l + 1)) as i32 + m) as usize;
    let coeff_index = ((l * (l + 1)) / 2 + m.abs() as u32) as usize;
    let base = if (l - m.abs() as u32) % 2 == 0 {0} else {1};
    move |x: f64| {
        let mut value = 0.0;
        for power in 0..(l - m.abs() as u32) / 2 {
            value += x.powi((2 * power + base) as i32) * LEGENDRE_COEFF[coeff_index][power as usize];
        }
        LEGENDRE_NORM[norm_index] * value * (1.0 - x * x).powf(m.abs() as f64 / 2.0)
    }
    
}

fn laguerre(n: u32, alpha: u32) -> impl Fn(f64) -> f64 {
    move |x: f64| {
        let mut value = 0.0;
        for i in 0..(n+1) {
            value += parity(i as i32) * choose(n+alpha, n-i) as f64 * x.powi(i as i32) / factorial(i) as f64;
        }
        value
    }
}

pub fn nlm(n: u32, l: u32, m: i32) -> impl Fn((f64, f64, f64)) -> Complex {
    #[cfg(debug)]
    {
        assert_eq!(m <= l && -m <= l, true);
    }

    // Position is in units of a_0

    let norm = ((2.0 / n as f64).powi(3) * factorial(n - l - 1) as f64 / (2 * n * factorial(n + l)) as f64).sqrt();
    let rho_scale = 2.0 / n as f64;
    let L = laguerre(n - l - 1, 2 * l + 1);
    let Y = ylm(l, m);

    move |x: (f64, f64, f64) | -> Complex {
        let r = (x.0 * x.0 + x.1 * x.1 + x.2 * x.2).sqrt();
        let r_asym = (x.0 * x.0 + x.1 * x.1).sqrt();
        let cos_theta = x.2 / r;
        let cos_phi = x.0 / r_asym;
        let sin_phi = x.1 / r_asym;

        let rho = r * rho_scale;

        Y(cos_theta, cos_phi, sin_phi) * norm * rho.powi(l as i32) * (-rho / 2.0).exp() * L(rho)
    }
}

pub fn energy(n: u32) -> f64 {
    // Units where position is a_0, mass is m_e, and t such that c = 1.
    - FINE_STRUCTURE * FINE_STRUCTURE / 2.0 / n as f64 / n as f64
}