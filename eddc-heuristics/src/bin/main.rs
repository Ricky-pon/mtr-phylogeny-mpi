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

    let mut encoded_reads = vec![vec![]; reads.len()];

    let params = string_decomposer::StringDecompParameters::new(1, -1, -1, -1);

    for (i, &read) in reads.iter().enumerate() {
        let decomposed = string_decomposer::DecomposedSeq::new(read, &units, &params);

        for e in decomposed.encoding() {
            encoded_reads[i].push(e.id());
        }
    }

    let score_matrix = eddc_heuristics::eddc_heuristic(&encoded_reads, &units);
    for v in &score_matrix {
        for &x in v {
            print!("{} ", x);
        }
        println!();
    }
    Ok(())
}
