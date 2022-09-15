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
    let score_matrix = eddc_heuristics::eddc_heuristic_parallel(&encoded_reads, &units);
    for (i, v) in score_matrix.iter().enumerate() {
        let scores: Vec<_> = v
            .iter()
            .enumerate()
            .map(|(j, &score)| {
                let len = ((reads[i].len() * reads[j].len()) as f64).sqrt();
                // Avoid zero div.
                let normalized_score = if 0.000001 < len { score / len } else { score };
                format!("{normalized_score}")
            })
            .collect();
        println!("{}", scores.join(" "));
    }
    Ok(())
}
