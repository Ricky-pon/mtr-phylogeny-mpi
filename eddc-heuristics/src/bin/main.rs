use eddc_heuristics::fasta;
use std::vec;
use string_decomposer;
fn main() -> std::io::Result<()> {
    let mut reads = vec![];
    let mut units = vec![];

    let args: Vec<_> = std::env::args().collect();
    let records = fasta::parse_into_vec(&args[1])?;
    for record in records.iter() {
        if record.desc().unwrap().trim().starts_with('(') {
            reads.push(record.seq());
        } else {
            units.push(record.seq());
        }
    }

    let params = string_decomposer::StringDecompParameters::new(1, -1, -1, -1);
    let encoded_reads: Vec<Vec<_>> = reads
        .iter()
        .map(|read| {
            string_decomposer::DecomposedSeq::new(read, &units, &params)
                .encoding()
                .iter()
                .map(|enc| enc.id())
                .collect()
        })
        .collect();
    let score_matrix = eddc_heuristics::eddc_heuristic(&encoded_reads, &units);
    for v in &score_matrix {
        for &x in v {
            print!("{} ", x);
        }
        println!();
    }
    Ok(())
}
