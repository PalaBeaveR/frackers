use std::{
    fmt::Display,
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign},
};

use num_bigint::{BigInt, Sign, ToBigInt};
use primes::{PrimeSet, Sieve};

#[derive(Debug, Clone)]
pub struct Fracker {
    pub numerator: BigInt,
    pub denominator: BigInt,
}

pub trait ToFracker {
    fn fr(self) -> Fracker;
}

impl ToFracker for (BigInt, BigInt) {
    fn fr(self) -> Fracker {
        Fracker {
            numerator: self.0,
            denominator: self.1,
        }
    }
}

impl ToFracker for (i32, i32) {
    fn fr(self) -> Fracker {
        Fracker {
            numerator: self.0.into(),
            denominator: self.1.into(),
        }
    }
}

impl ToFracker for (i64, i64) {
    fn fr(self) -> Fracker {
        Fracker {
            numerator: self.0.into(),
            denominator: self.1.into(),
        }
    }
}

impl ToFracker for i32 {
    fn fr(self) -> Fracker {
        Fracker {
            numerator: self.into(),
            denominator: 1.into(),
        }
    }
}

impl MulAssign for Fracker {
    fn mul_assign(&mut self, rhs: Self) {
        self.numerator *= rhs.numerator;
        self.denominator *= rhs.denominator;
    }
}

impl Mul for Fracker {
    type Output = Fracker;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl Mul<Fracker> for i32 {
    type Output = Fracker;

    fn mul(self, mut rhs: Fracker) -> Self::Output {
        rhs.numerator *= self;
        rhs
    }
}

impl Div for Fracker {
    type Output = Fracker;

    fn div(mut self, rhs: Self) -> Self::Output {
        self.numerator *= rhs.denominator;
        self.denominator *= rhs.numerator;
        self
    }
}

impl Div<BigInt> for Fracker {
    type Output = Fracker;

    fn div(mut self, rhs: BigInt) -> Self::Output {
        self.denominator *= rhs;
        self
    }
}

impl Div<&Fracker> for &Fracker {
    type Output = Fracker;

    fn div(self, rhs: &Fracker) -> Self::Output {
        Fracker {
            numerator: &self.numerator * &rhs.denominator,
            denominator: &self.denominator * &rhs.numerator,
        }
    }
}

impl Add for Fracker {
    type Output = Fracker;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl AddAssign for Fracker {
    fn add_assign(&mut self, rhs: Self) {
        self.numerator = self.numerator.clone() * rhs.denominator.clone()
            + rhs.numerator * self.denominator.clone();
        self.denominator *= rhs.denominator;
    }
}

impl Fracker {
    pub fn abs(mut self) -> Fracker {
        if self.numerator.sign() == Sign::Minus {
            self.numerator = -self.numerator
        }

        if self.denominator.sign() == Sign::Minus {
            self.denominator = -self.denominator
        }

        self
    }

    pub fn simplify(mut self) -> Fracker {
        for n in Sieve::new().iter().take(2000) {
            loop {
                if self.numerator.clone() % n == 0.to_bigint().unwrap()
                    && self.denominator.clone() % n == 0.to_bigint().unwrap()
                {
                    self.numerator /= n;
                    self.denominator /= n;
                    //println!("Reduced by {}", n);
                } else {
                    break;
                }
            }
        }

        self
    }

    pub fn gcd(mut a: BigInt, mut b: BigInt) -> BigInt {
        loop {
            if b == 0.into() {
                return a;
            }
            std::mem::swap(&mut a, &mut b);
            b %= a.clone();
        }
    }

    pub fn lcm(a: BigInt, b: BigInt) -> BigInt {
        if a > b {
            a.clone() / Fracker::gcd(a, b.clone()) * b
        } else {
            b.clone() / Fracker::gcd(a.clone(), b) * a
        }
    }

    /// Inverts probability
    /// Equivalent to 1 - fraction
    pub fn inv_prob(mut self) -> Fracker {
        self.numerator = self.denominator.clone() - self.numerator;
        self
    }

    pub fn pow(mut self, power: u32) -> Fracker {
        self.numerator = self.numerator.pow(power);
        self.denominator = self.denominator.pow(power);
        self
    }

    pub fn sqrt_fast(mut self) -> Fracker {
        self.numerator = self.numerator.sqrt();
        self.denominator = self.denominator.sqrt();
        self
    }

    pub fn halve(mut self) -> Fracker {
        self.denominator *= 2;
        self
    }

    pub fn halve_lossy(mut self) -> Fracker {
        self.numerator = self.numerator / 2;
        self
    }

    pub fn add_lcm(mut self, other: Fracker) -> Fracker {
        let lcm = Fracker::lcm(self.denominator.clone(), other.denominator.clone());
        self.numerator = &self.numerator * (&lcm / &self.denominator)
            + other.numerator * (&lcm / other.denominator);
        self.denominator = lcm;
        self
    }

    pub fn sub_lcm(mut self, other: Fracker) -> Fracker {
        let lcm = Fracker::lcm(self.denominator.clone(), other.denominator.clone());
        self.numerator = &self.numerator * (&lcm / &self.denominator)
            - other.numerator * (&lcm / other.denominator);
        self.denominator = lcm;
        self
    }

    pub fn sqrt_iter(self, iter: u32) -> Fracker {
        let mut estimate = self.clone().sqrt_fast();
        estimate.numerator += estimate.denominator.clone();

        let mut x = estimate.clone();

        for _ in 0..iter {
            let tempx = x.clone() + (&self / &x);
            x = tempx.halve_lossy();
        }

        x
    }

    pub fn sqrt_iter_lcm(self, iter: u32) -> Fracker {
        let mut estimate = self.clone().sqrt_fast();
        estimate.numerator += estimate.denominator.clone();

        let mut x = estimate.clone();

        for _ in 0..iter {
            let tempx = x.clone().add_lcm(&self / &x);
            x = tempx.halve_lossy();
        }

        x
    }

    pub fn flip(mut self) -> Fracker {
        std::mem::swap(&mut self.numerator, &mut self.denominator);
        self
    }

    pub fn is_negative(&self) -> bool {
        (self.numerator.sign() == Sign::Minus) ^ (self.denominator.sign() == Sign::Minus)
    }
}

impl Display for Fracker {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut clone = self.clone();
        if self.is_negative() {
            write!(f, "-")?;
            if self.numerator.sign() == Sign::Minus {
                clone.numerator *= -1;
            } else {
                clone.denominator *= -1
            }
        }
        let digit = &clone.numerator / &clone.denominator;

        write!(f, "{digit}")?;

        let mut remainder = &clone.numerator - digit * &clone.denominator;

        if remainder == 0.to_bigint().unwrap() {
            return writeln!(f);
        }

        write!(f, ".")?;

        for _ in 0..f.precision().unwrap_or(5) {
            remainder *= 10;
            let digit = &remainder / &clone.denominator;

            write!(f, "{digit}")?;

            remainder = remainder - digit * &clone.denominator;
        }

        writeln!(f)
    }
}

impl SubAssign for Fracker {
    fn sub_assign(&mut self, rhs: Self) {
        self.numerator = &self.numerator * &rhs.denominator - &rhs.numerator * &self.denominator;
        self.denominator *= rhs.denominator;
    }
}

impl SubAssign<&Fracker> for Fracker {
    fn sub_assign(&mut self, rhs: &Self) {
        self.numerator = &self.numerator * &rhs.denominator - &rhs.numerator * &self.denominator;
        self.denominator *= rhs.denominator.clone();
    }
}

impl Sub for Fracker {
    type Output = Fracker;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl Sub<&Fracker> for Fracker {
    type Output = Fracker;

    fn sub(mut self, rhs: &Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl Neg for Fracker {
    type Output = Fracker;

    fn neg(mut self) -> Self::Output {
        self.numerator = -self.numerator;
        self
    }
}
