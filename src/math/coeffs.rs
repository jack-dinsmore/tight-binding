use super::complex::*;
use super::hydrogen::*;
use std::collections::HashMap;

const GRID_SIZE: f64 = 0.1;
const GRID_LIMIT: f64 = 5.0;

pub struct ProtonLayout<const PROTONS: usize> {
    zs: [u32; PROTONS],
    poses: [(f64, f64, f64); PROTONS],
}

impl<const PROTONS: usize> ProtonLayout<PROTONS> {
    pub fn new(zs: [u32; PROTONS], poses: [(f64, f64, f64); PROTONS]) -> Self {
        Self {
            zs,
            poses,
        }
    }

    pub fn single(z: u32) -> ProtonLayout<1> {
        ProtonLayout::<1> {
            zs: [z],
            poses: [(0.0, 0.0, 0.0)]
        }
    }

    pub fn zs(&self, i: usize) -> u32 {
        self.zs[i]
    }

    pub fn pos(&self, i: usize) -> (f64, f64, f64) {
        self.poses[i]
    }

    pub fn count(&self) -> u32 {
        let mut val = 0;
        for z in &self.zs {
            val += z;
        }
        val
    }
}

fn sub(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (a.0 - b.0, a.1 - b.1, a.2 - b.2)
}

fn add(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (a.0 + b.0, a.1 + b.1, a.2 + b.2)
}

fn norm(a: (f64, f64, f64)) -> f64 {
    (a.0 * a.0 + a.1 * a.1 + a.2 * a.2).sqrt()
}

fn integrate(f: impl Fn((f64, f64, f64)) -> Complex) -> Complex {
    let mut value = Complex::ZERO;
    let mut x = -GRID_LIMIT;
    while x < GRID_LIMIT {
        let mut y = -GRID_LIMIT;
        while y < GRID_LIMIT {
            let mut z = -GRID_LIMIT;
            while z < GRID_LIMIT {
                value += f((x, y, z));
                z += GRID_SIZE;
            }
            y += GRID_SIZE;
        }
        x += GRID_SIZE;
    }
    value * GRID_SIZE * GRID_SIZE * GRID_SIZE
}

pub fn make_s<const PROTONS: usize, const N_MAX: usize>(layout: &ProtonLayout<PROTONS>) -> HashMap<(u32, u32, u32, u32, u32, u32, u32, u32), Complex> {
    let mut map = HashMap::new();
    for a in 0..PROTONS {
        for ap in 0..PROTONS {
            for n in 1..N_MAX {
                for l in 0..n {
                    for m in (-(l as isize))..(l+1) as isize {
                        let wf = nlm(n, l, m);
                        for np in 1..N_MAX {
                            for lp in 0..np {
                                for mp in (-(lp as isize))..(lp+1) as isize {
                                    let wfp = nlm(np, lp, mp);
                                    let key = (a as u32, ap as u32, n as u32, l as u32, m as u32, np as u32, lp as u32, mp as u32);
                                    if map.contains_key(&key) {continue;}
                                    let value = integrate(|r: (f64, f64, f64)| {
                                        wf(r) * wfp(add(sub(r, layout.pos(a)), layout.pos(ap)))
                                    });
                                    map.insert(key, value);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    map
}

pub fn make_q<const PROTONS: usize, const N_MAX: usize>(layout: &ProtonLayout<PROTONS>) -> HashMap<(u32, u32, u32, u32, u32, u32, u32, u32, u32, u32, u32, u32), Complex> {
    let mut hmap = HashMap::new();
    for na in 1..N_MAX {
        for la in 0..na {
            for ma in (-(la as isize))..(la+1) as isize {
                let wfa = nlm(na as u32, la as u32, ma as u32);
                for nap in 1..N_MAX {
                    for lap in 0..nap {
                        for map in (-(lap as isize))..(lap+1) as isize {
                            let wfap = nlm(nap as u32, lap as u32, map as u32);
                            for nb in 1..N_MAX {
                                for lb in 0..nb {
                                    for mb in (-(lb as isize))..(lb+1) as isize {
                                        let wfb = nlm(nb as u32, lb as u32, mb as u32);
                                        for nbp in 1..N_MAX {
                                            for lbp in 0..nbp {
                                                for mbp in (-(lbp as isize))..(lbp+1) as isize {
                                                    let wfbp = nlm(nbp as u32, lbp as u32, mbp as u32);
                                                    let key = (na as u32, la as u32, ma as u32, nap as u32, lap as u32, map as u32,
                                                               nb as u32, lb as u32, mb as u32, nbp as u32, lbp as u32, mbp as u32);
                                                    if hmap.contains_key(&key) {continue;}
                                                    let value = integrate(|ra: (f64, f64, f64)| {
                                                        let ra_scale = ;
                                                        integrate(|rb:(f64, f64, f64)| {
                                                            let rb_scale = ;
                                                            wfa(ra) * wfap(ra) * wfb(rb) * wfbp(rb) * 1 / norm(sub(ra_scale, rb_scale))
                                                        })
                                                    });
                                                    hmap.insert(key, value);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    hmap
}