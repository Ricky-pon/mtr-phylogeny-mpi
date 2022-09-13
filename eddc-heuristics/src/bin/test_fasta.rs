fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let records = eddc_heuristics::fasta::parse_into_vec(&args[1])?;
    for record in records.iter() {
        println!("{}\t{}", record.id(), record.seq().len());
    }
    Ok(())
}
