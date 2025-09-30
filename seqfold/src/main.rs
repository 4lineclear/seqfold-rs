use clap::Parser;

/// Predict the minimum free energy (kcal/mol) of a nucleic acid sequence
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[arg()]
    seq: String,

    /// The temperature in celsius
    #[arg(short, long)]
    temp: Option<f64>,
}

fn main() {
    let args = Args::parse();
    let values = seqfold_rs::fold(args.seq.as_bytes(), args.temp);
    let mfe: f64 = values.iter().map(|v| v.e).sum();
    println!("{mfe:.2}");
}
