//! String decomposer -- clone
//! # Example
//! ```rust
//! let seq = b"TTTACTACTACTGCGCGCTTT";
//! let units = vec![b"ACT", b"GC"];
//! let param = StrindDecompParameters::new(1,-1,-1,-1);
//! let decomposed = DecomposedSeq::new(seq, &units, &param);
//! ```
//! # Reference
//! - original paper
//! - Morishita-sensei's implementation
//!
//!

#[derive(Debug, Clone)]
pub struct DecomposedSeq {}

pub struct StringDecompParameters {
    match_score: i64,
    mismatch_score: i64,
    // Insertion to the *unit*. Should be negative.
    ins_score: i64,
    // Deletion from the *unit*. Should be negative.
    del_score: i64,
}

impl StringDecompParameters {
    pub const fn new(mat: i64, mism: i64, ins: i64, del: i64) -> Self {
        Self {
            match_score: mat,
            mismatch_score: mism,
            ins_score: ins,
            del_score: del,
        }
    }
}
impl DecomposedSeq {
    pub fn new<U: std::borrow::Borrow<[u8]>>(
        seq: &[u8],
        units: &[U],
        param: &StringDecompParameters,
    ) -> Self {
        string_decomposer(seq, units, param);
        todo!()
    }
}

// The first element, each "row" of the DP-table
type Row = (i64, Vec<Vec<i64>>);
#[derive(Debug, Clone)]
struct StrDecompDP {
    inner: Vec<Row>,
}

impl StrDecompDP {
    fn new(seq: &[u8], units: &[&[u8]]) -> Self {
        let row: Vec<_> = units.iter().map(|u| vec![0; u.len()]).collect();
        let inner = vec![(0, row.clone()); seq.len() + 1];
        Self { inner }
    }
    fn get_row_mut(&mut self, i: usize) -> (&mut i64, &mut [Vec<i64>]) {
        self.inner
            .get_mut(i)
            .map(|(x, y)| (x, y.as_mut_slice()))
            .unwrap()
    }
    fn get_mut(&mut self, read_pos: usize, unit_id: usize, unit_pos: usize) -> &mut i64 {
        match unit_pos {
            0 => &mut self.inner[read_pos].0,
            _ => &mut (self.inner[read_pos].1)[unit_id][unit_pos - 1],
        }
    }
    fn get(&self, read_pos: usize, unit_id: usize, unit_pos: usize) -> i64 {
        match unit_pos {
            0 => self.inner[read_pos].0,
            _ => self.inner[read_pos].1[unit_id][unit_pos - 1],
        }
    }
    fn last_max_of(&self, r_pos: usize) -> (i64, usize, usize) {
        self.inner[r_pos]
            .1
            .iter()
            .enumerate()
            .map(|(u_id, row)| (u_id, row.iter().enumerate().last()))
            .filter_map(|(u_id, x)| x.map(|(u_pos, score)| (*score, u_id, u_pos + 1)))
            .max_by_key(|x| x.0)
            .unwrap()
    }
}

// Do global  alignment...
fn string_decomposer<U: std::borrow::Borrow<[u8]>>(
    seq: &[u8],
    units: &[U],
    param: &StringDecompParameters,
) -> (i64, Vec<(usize, usize, usize)>) {
    let units: Vec<_> = units.iter().map(|u| u.borrow()).collect();
    // Initialization.
    let mut dp_table = StrDecompDP::new(seq, &units);
    {
        let (_, rows) = dp_table.get_row_mut(0);
        for row in rows.iter_mut() {
            for (i, slot) in row.iter_mut().enumerate() {
                *slot = param.del_score * (i as i64 + 1);
            }
        }
    }
    // DP
    for (i, s_base) in seq.iter().enumerate().map(|(i, b)| (i + 1, b)) {
        for (u_id, unit) in units.iter().enumerate() {
            for (j, u_base) in unit.iter().enumerate().map(|(j, b)| (j + 1, b)) {
                let mat_score = match s_base == u_base {
                    true => param.match_score,
                    false => param.mismatch_score,
                };
                let mat_score = dp_table.get(i - 1, u_id, j - 1) + mat_score;
                let del_score = dp_table.get(i, u_id, j - 1) + param.del_score;
                let ins_score = dp_table.get(i - 1, u_id, j) + param.ins_score;
                *dp_table.get_mut(i, u_id, j) = mat_score.max(del_score).max(ins_score);
            }
        }
        let (wrap, _u_id, _u_pos) = dp_table.last_max_of(i);
        *dp_table.get_mut(i, 0, 0) = wrap.max(0);
    }
    // Traceback
    let mut alnpath = vec![];
    let (mut score, mut r_pos, mut u_id, mut u_pos) = (0..seq.len() + 1)
        .map(|r_pos| {
            let (score, u_id, u_pos) = dp_table.last_max_of(r_pos);
            (score, r_pos, u_id, u_pos)
        })
        .max_by_key(|x| x.0)
        .unwrap();
    eprintln!("{score}");
    let opt_score = score;
    while !(score == 0 && u_pos == 0) {
        eprintln!("{r_pos},{u_id},{u_pos}");
        alnpath.push((r_pos - 1, u_id, u_pos - 1));
        let u_base = units[u_id][u_pos - 1];
        let r_base = seq[r_pos - 1];
        let mat_score = if u_base == r_base {
            param.match_score
        } else {
            param.mismatch_score
        };
        let mat_score = dp_table.get(r_pos - 1, u_id, u_pos - 1) + mat_score;
        let del_score = dp_table.get(r_pos, u_id, u_pos - 1) + param.del_score;
        let ins_score = dp_table.get(r_pos - 1, u_id, u_pos) + param.ins_score;
        if score == mat_score {
            r_pos -= 1;
            u_pos -= 1;
        } else if score == del_score {
            u_pos -= 1;
        } else if score == ins_score {
            r_pos -= 1;
        } else {
            unreachable!()
        }
        score = dp_table.get(r_pos, u_id, u_pos);
        // If reached the start of the unit, wrap around.
        if u_pos == 0 && r_pos != 0 && score != 0 {
            let (prev_max, next_u_id, next_u_pos) = dp_table.last_max_of(r_pos);
            assert_eq!(prev_max, score);
            u_id = next_u_id;
            u_pos = next_u_pos;
        }
    }
    alnpath.reverse();
    (opt_score, alnpath)
}

#[cfg(test)]
mod tests {
    use super::*;
    const PARAM: StringDecompParameters = StringDecompParameters::new(1, -1, -1, -1);
    #[test]
    fn string_decomposer_test() {
        let read = b"ACGACGACG";
        let units = vec![b"ACG".as_slice()];
        let (score, ops) = string_decomposer(read, &units, &PARAM);
        assert_eq!(score, 9);
        let path = vec![
            (0, 0, 0),
            (1, 0, 1),
            (2, 0, 2),
            (3, 0, 0),
            (4, 0, 1),
            (5, 0, 2),
            (6, 0, 0),
            (7, 0, 1),
            (8, 0, 2),
        ];
        assert_eq!(ops, path);
    }
    #[test]
    fn string_decomposer_test_2() {
        let read = b"TTTACGCCGACGTTT";
        let units = vec![b"ACG".as_slice(), b"CCG".as_slice()];
        let (score, ops) = string_decomposer(read, &units, &PARAM);
        assert_eq!(score, 9);
        let path = vec![
            (3, 0, 0),
            (4, 0, 1),
            (5, 0, 2),
            (6, 1, 0),
            (7, 1, 1),
            (8, 1, 2),
            (9, 0, 0),
            (10, 0, 1),
            (11, 0, 2),
        ];
        assert_eq!(path, ops);
    }
    #[test]
    fn string_decomposer_test_3() {}
    #[test]
    fn string_decomposer_test_4() {}
}
