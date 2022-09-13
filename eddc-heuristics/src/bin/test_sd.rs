fn main() -> std::io::Result<()> {
    let read = b"TTTACGACGACGTTT";
    let units = vec![b"ACG".as_slice(), b"CCC".as_slice()];
    use string_decomposer::*;
    let param = StringDecompParameters::new(1, -1, -1, -1);
    println!("{:?}", DecomposedSeq::new(read, &units, &param));
    Ok(())
}
