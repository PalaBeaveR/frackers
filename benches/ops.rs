use criterion::{black_box, criterion_group, criterion_main, Criterion};
use frackers::{Fracker, ToFracker};
use num_bigint::{BigInt, ToBigInt};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("add", |b| {
        b.iter(|| {
            black_box(0..10000)
                .map(ToFracker::fr)
                .map(|f| f * (1, 523).fr())
                .fold(0.fr(), |acc, el| acc + el);
        })
    });
    c.bench_function("add_lcm", |b| {
        b.iter(|| {
            black_box(0..10000)
                .map(ToFracker::fr)
                .map(|f| f * (1, 523).fr())
                .fold(0.fr(), |acc, el| acc.add_lcm(el));
        })
    });

    c.bench_function("sqrt_10", |b| {
        b.iter(|| {
            black_box(10000000..10000010)
                .map(ToFracker::fr)
                .map(|f| f * (1, 523).fr())
                .fold(0.fr(), |acc, el| acc + el.sqrt_iter(10));
        })
    });

    c.bench_function("sqrt_12", |b| {
        b.iter(|| {
            black_box(10000000..10000010)
                .map(ToFracker::fr)
                .map(|f| f * (1, 523).fr())
                .fold(0.fr(), |acc, el| acc + el.sqrt_iter(12));
        })
    });

    c.bench_function("erf", |b| {
        let sqrt_of_two = 2.fr().sqrt_iter_lcm(10);
        let pi = calculate_pi(sqrt_of_two.clone(), 1);
        let sqrt_pi = pi.clone().sqrt_iter_lcm(10);
        //let one_over_sqrt_of_two_pi = (2 * pi).sqrt_iter(10).flip();
        let two_over_sqrt_pi = 2 * sqrt_pi.clone().flip();
        b.iter(|| {
            black_box(0..2)
                .map(ToFracker::fr)
                .map(|f| f * (1, 523).fr())
                .fold(0.fr(), |acc, el| {
                    acc.add_lcm(erf(two_over_sqrt_pi.clone(), el, 60))
                });
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

pub fn erf(two_over_sqrt_pi: Fracker, x: Fracker, iters: u32) -> Fracker {
    let mut approx = 0.fr();
    for n in 0..iters {
        let two_n_plus_one = 2 * n + 1;
        let x_power_two_n_plus_one = x.clone().pow(two_n_plus_one);

        let n_factorial = factorial(n.try_into().unwrap());
        let frac = x_power_two_n_plus_one / (two_n_plus_one * n_factorial);
        if n % 2 == 0 {
            approx = approx.add_lcm(frac);
        } else {
            approx = approx.sub_lcm(frac);
        }

        //println!("{:.50}", approx)
    }

    two_over_sqrt_pi * approx
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
        acc += iter;
    }
    (pi * acc).flip()
}

fn factorial(n: i32) -> BigInt {
    (1..=n).map(BigInt::from).product()
}
