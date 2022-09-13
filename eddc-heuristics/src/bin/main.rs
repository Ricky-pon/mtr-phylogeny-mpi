use eddc_heuristics::fasta;
use std::{io, vec};
use string_decomposer;
fn main() {
    let mut reads = vec![];
    let mut units = vec![];
    {
        loop {
            let mut seq_summary = String::new();
            let mut seq = String::new();
            match io::stdin().read_line(&mut seq_summary) {
                Ok(_) => match io::stdin().read_line(&mut seq) {
                    Ok(length) => {
                        if length == 0 {
                            break;
                        }
                        seq.retain(|c| c != '\n');
                        if seq_summary.into_bytes()[2] == '(' as u8 {
                            reads.push(seq);
                        } else {
                            units.push(seq);
                        }
                    }
                    Err(error) => println!("error: {error}"),
                },
                Err(_) => {
                    break;
                }
            }
        }
    }

    let mut encoded_reads = vec![vec![]; reads.len()];

    let reads_u8: Vec<&[u8]> = reads.iter().map(|seq| (&*seq).as_bytes()).collect();
    let units_u8: Vec<&[u8]> = units.iter().map(|seq| (&*seq).as_bytes()).collect();
    let params = string_decomposer::StringDecompParameters::new(1, -1, -1, -1);

    for (i, &read) in reads_u8.iter().enumerate() {
        let decomposed = string_decomposer::DecomposedSeq::new(read, &units_u8, &params);

        for e in decomposed.encoding() {
            encoded_reads[i].push(e.id());
        }
        // println!("{:?}", decomposed);
    }

    let score_matrix = eddc_heuristics::eddc_heuristic(&encoded_reads, &units_u8);
    for v in &score_matrix {
        for &x in v {
            print!("{} ", x);
        }
        println!();
    }
}
