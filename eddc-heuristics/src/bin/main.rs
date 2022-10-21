use eddc_heuristics::fasta;
use std::collections::{BTreeMap, BTreeSet};
use std::vec;
use string_decomposer::{self, Op, StringDecompParameters};
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

    let mut unit_variants = vec![];
    enum_variants(&reads, &units, &params, &mut unit_variants);
    for variant in &unit_variants {
        units.push(variant);
    }

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

fn enum_variants(
    reads: &[&[u8]],
    units: &[&[u8]],
    params: &StringDecompParameters,
    unit_variants: &mut Vec<Vec<u8>>,
) {
    let mut count_variant: BTreeMap<Vec<u8>, usize> = BTreeMap::new();
    let mut read_num: BTreeMap<Vec<u8>, usize> = BTreeMap::new(); // variant v -> # of reads containing v
    for read in reads {
        let mut variant_set: BTreeSet<Vec<u8>> = BTreeSet::new();
        for enc in string_decomposer::DecomposedSeq::new(read, units, params).encoding() {
            if enc.ops().iter().all(|op| matches!(op, Op::Match(_, _))) {
                continue;
            }

            let mut variant = vec![];
            for op in enc.ops() {
                if let Op::Match(x, _) = op {
                    variant.push(*x);
                }
                if let Op::Ins(x) = op {
                    variant.push(*x);
                }
            }

            count_variant.insert(
                variant.clone(),
                match count_variant.get(&variant) {
                    Some(val) => val + 1,
                    None => 1,
                },
            );
            variant_set.insert(variant.clone());
        }

        for variant in &variant_set {
            read_num.insert(
                variant.clone(),
                match read_num.get(variant) {
                    Some(val) => val + 1,
                    None => 1,
                },
            );
        }
    }

    let mut variant_vec = Vec::from_iter(count_variant.iter());
    variant_vec.sort_by(|a, b| b.1.cmp(a.1));
    variant_vec.truncate(5);
    for (variant, &num) in variant_vec {
        if num > 1 {
            unit_variants.push(variant.clone());
        }
    }
}
