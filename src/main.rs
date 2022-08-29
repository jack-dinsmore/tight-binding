mod math;

use math::{ProtonLayout, energies};

fn main() {
    const N_MAX: usize = 1;
    let hydrogen = ProtonLayout::single(1);
    let es = energies::<1, N_MAX>(&hydrogen);
    
    println!("For hydrogen with N_MAX={}", N_MAX);
    for e in es {
        println!("{}", e);
    }
}