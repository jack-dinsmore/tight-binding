use std::ops::{Add, Mul, Sub, AddAssign, SubAssign, MulAssign};
pub struct Complex {
    r: f64,
    i: f64,
}

impl Complex {
    pub const I: Self = Self {r: 0.0, i: 1.0};
    pub const ZERO: Self = Self {r: 0.0, i:0.0};

    pub fn real(&self) -> f64 {
        self.r
    }
    pub fn imag(&self) -> f64 {
        self.i
    }
    pub fn conj(&self) -> Self {
        Self {r: self.r, i: -self.i}
    }
    pub fn norm(&self) -> f64 {
        self.r * self.r + self.i * self.i
    }
}

impl Add for Complex {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            r: self.r + other.r,
            i: self.i + other.i,
        }
    }
}

impl AddAssign for Complex {
    fn add_assign(&mut self, other: Self) {
        self.r += other.r;
        self.i += other.i;
    }
}

impl Sub for Complex {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            r: self.r - other.r,
            i: self.i - other.i,
        }
    }
}

impl SubAssign for Complex {
    fn sub_assign(&mut self, other: Self) {
        self.r -= other.r;
        self.i -= other.i;
    }
}

impl Mul for Complex {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self {
            r: self.r * other.r - self.i * other.i,
            i: self.r * other.i + self.i * other.r,
        }
    }
}

impl MulAssign for Complex {
    fn mul_assign(&mut self, other: Self) {
        *self = Self {
            r: self.r * other.r - self.i * other.i,
            i: self.r * other.i + self.i * other.r,
        };
    }
}

impl Mul<f64> for Complex {
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        Self {
            r: self.r * other,
            i: self.i * other,
        }
    }
}

impl MulAssign<f64> for Complex {
    fn mul_assign(&mut self, other: f64) {
        self.r *= other;
        self.i *= other;
    }
}