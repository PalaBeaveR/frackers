use std::time::Instant;

use cached::proc_macro::cached;
use frackers::{Fracker, ToFracker};
use num_bigint::{BigInt, ToBigInt};

fn main() {
    let time = Instant::now();
    let sqrt_of_two = 2.fr().sqrt_iter_lcm(10);
    println!("Finished part 1");
    let pi = calculate_pi(sqrt_of_two.clone(), 10);
    println!("Finished part 1");
    let sqrt_pi = pi.clone().sqrt_iter(10);
    println!("Finished part 1");
    //let one_over_sqrt_of_two_pi = (2 * pi).sqrt_iter(10).flip();
    let two_over_sqrt_pi = 2 * sqrt_pi.clone().flip();
    println!("Finished part 1");
    let one_over_sqrt_of_two_pi = (sqrt_of_two.clone() * sqrt_pi.clone()).flip();
    println!("Finished part 1");
    let sqrt_of_pi_over_two = sqrt_pi.clone() / sqrt_of_two.clone();
    println!("Finished part 1");

    //let fr = (624 * 13, 402 * 13).fr();
    //println!("{:.50}", fr);
    //println!("{:.50}", fr.simplify());
    //return;

    let n = 1000;
    let p = (1, 10).fr();
    let q = p.clone().inv_prob();
    let k1 = 301;
    let k2 = 308;
    println!("Finished part 2");

    let np = n * p.clone();
    let sqrt_npq = (np.clone() * q).sqrt_iter(10);
    let k1_top = k1.fr() - np.clone();
    let k2_top = k2.fr() - np.clone();
    let x1 = k1_top / sqrt_npq.clone();
    let x2 = k2_top / sqrt_npq.clone();
    println!("Finished part 3");
    println!("{x1:.10}");
    println!("{x2:.10}");

    println!(
        "Approximate {:.150}",
        double_laplace(
            one_over_sqrt_of_two_pi,
            sqrt_of_pi_over_two,
            sqrt_of_two,
            two_over_sqrt_pi,
            x1,
            x2
        )
    );

    println!("Accurate: {:.200}", chance_accurate(k1, k2, n, p));

    println!("{} seconds", time.elapsed().as_secs_f64());

    //println!("{:?}", combinations(4, 10));
    //println!("{:.50}", chance_accurate(4, 5, 6, (1, 10).fr()));
    //println!("{:.50}", (100000, 3000).fr().sqrt_fast());
    //println!("{:.50}", (100000, 3000).fr().sqrt_iter(10));
    //println!("{:.50}", calculate_pi(2.fr().sqrt_iter(10), 10));
    //println!("{:.50}", erf(two_over_sqrt_pi, 5.fr(), 100))
}

pub fn chance_accurate(k_start: i32, k_end: i32, n: i32, p: Fracker) -> Fracker {
    let q = p.clone().inv_prob();

    let mut probability = 0.fr();

    for k in k_start..=k_end {
        probability += combinations(k, n)
            * p.clone().pow(k.try_into().unwrap())
            * q.clone().pow((n - k).try_into().unwrap())
    }

    probability
}
pub fn double_laplace(
    one_over_sqrt_of_two_pi: Fracker,
    sqrt_of_pi_over_two: Fracker,
    sqrt_of_two: Fracker,
    two_over_sqrt_pi: Fracker,
    x1: Fracker,
    x2: Fracker,
) -> Fracker {
    let zero_erf = erf(0.fr(), 100);
    println!("Zero erf: {zero_erf:.50}");
    let erf1 = erf(x1 / sqrt_of_two.clone(), 100);
    println!("erf1: {erf1:.50}");
    let erf2 = erf(x2 / sqrt_of_two, 100);
    println!("erf2: {erf2:.50}");

    let integral_part =
        sqrt_of_pi_over_two * two_over_sqrt_pi * (erf2 - zero_erf.clone() - (erf1 - zero_erf));
    println!("integral part: {integral_part:.50}");
    one_over_sqrt_of_two_pi * integral_part
}

pub fn erf(x: Fracker, iters: u32) -> Fracker {
    let mut approx = 0.fr();
    let mut pow = x.clone();
    let mut pow_increment = x.clone() * x.clone();
    for n in 0..iters {
        println!("{n}, {approx:.50}");
        let two_n_plus_one = 2 * n + 1;
        let x_power_two_n_plus_one = pow.clone();
        pow *= pow_increment.clone();

        let n_factorial = factorial(n.try_into().unwrap());
        let frac = x_power_two_n_plus_one / (two_n_plus_one * n_factorial);
        if n % 2 == 0 {
            approx = approx.add_lcm(frac);
        } else {
            approx = approx.sub_lcm(frac);
        }

        //println!("{:.50}", approx)
    }

    approx
}

pub fn erf_geometric(x: Fracker, iters: u32) -> Fracker {
    let mut approx = 0.fr();
    let mut pow = x.clone();
    let mut pow_increment = x.clone() * x.clone();
    for n in (0..iters).map(|n| 2u32.pow(n)) {
        println!("{n}, {approx:.50}");
        let two_n_plus_one = 2 * n + 1;
        let x_power_two_n_plus_one = x.clone().pow(two_n_plus_one);
        //let x_power_two_n_plus_one = pow.clone();
        //pow *= pow_increment.clone();

        let n_factorial = factorial(n.try_into().unwrap());
        let frac = x_power_two_n_plus_one / (two_n_plus_one * n_factorial);
        if n % 2 == 0 {
            approx = approx.add_lcm(frac);
        } else {
            approx = approx.sub_lcm(frac);
        }

        //println!("{:.50}", approx)
    }

    approx
}

pub fn combinations(this: i32, out_of: i32) -> Fracker {
    let largest = this.max(out_of - this);

    Fracker {
        numerator: ((largest + 1)..=out_of).map(BigInt::from).product(),
        denominator: (1..=(out_of - largest)).map(BigInt::from).product(),
    }
}

pub fn calculate_pi(sqrt_of_two: Fracker, iterations: u32) -> Fracker {
    let mut pi = (2, 9801).fr() * sqrt_of_two;
    let mut buff = 1103.to_bigint().unwrap();
    let mut acc = 0.fr();
    for i in 0..iterations {
        //if i != 0 {
        //    buff += 26390.to_bigint().unwrap();
        //}

        let a = factorial((4 * i).try_into().unwrap());
        let b = factorial(i.try_into().unwrap()).pow(4);
        let c = 396.to_bigint().unwrap().pow(4 * i);
        let iter = (a * (buff.clone() + 26390.to_bigint().unwrap() * i), b * c).fr();
        acc = acc.add_lcm(iter);
    }
    (pi * acc).flip()
}

#[cached]
fn factorial(n: i32) -> BigInt {
    (1..=n).map(BigInt::from).product()
}

fn erf_approx(mut x: Fracker) -> Fracker {
    println!("x: {x:.50}");
    let negative = if x.is_negative() {
        x = x.abs();
        true
    } else {
        false
    };
    println!("x: {x:.50}");

    let a1 = (705230784i64, 10000000000i64).fr();
    println!("{a1:.50}");
    let a2 = (422820123i64, 10000000000i64).fr();
    println!("{a2:.50}");
    let a3 = (92705272i64, 10000000000i64).fr();
    println!("{a3:.50}");
    let a4 = (1520143i64, 10000000000i64).fr();
    println!("{a4:.50}");
    let a5 = (2765672i64, 10000000000i64).fr();
    println!("{a5:.50}");
    let a6 = (430638i64, 10000000000i64).fr();
    println!("{a6:.50}");

    println!(
        "POS ERF: {:.50}",
        (1.fr()
            + a1.clone() * x.clone()
            + a2.clone() * x.clone().pow(2)
            + a3.clone() * x.clone().pow(3)
            + a4.clone() * x.clone().pow(4)
            + a5.clone() * x.clone().pow(5)
            + a6.clone() * x.clone().pow(6))
        .pow(16)
        .flip()
    );

    let out = 1.fr()
        - (1.fr()
            + a1 * x.clone()
            + a2 * x.clone().pow(2)
            + a3 * x.clone().pow(3)
            + a4 * x.clone().pow(4)
            + a5 * x.clone().pow(5)
            + a6 * x.pow(6))
        .pow(16)
        .flip();

    if negative {
        -out
    } else {
        out
    }
}
