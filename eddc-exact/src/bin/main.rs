use eddc_exact::fasta;
use std::vec;
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

    println!("{} {}", reads.len(), units.len());

    let score_matrix = eddc_exact::eddc_exact_parallel(&reads, &units);

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
