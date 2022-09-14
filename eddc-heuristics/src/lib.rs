use std::vec;

pub mod fasta;
pub struct EDDCParameters {
    match_score: Vec<Vec<f64>>,
    // Insertion to the *unit*. Should be negative.
    ins_score: Vec<f64>,
    // Deletion from the *unit*. Should be negative.
    del_score: Vec<f64>,
    // Match score of unit vs concat_of_units
    multiunit_match_score: Vec<Vec<Vec<Vec<f64>>>>,
}

impl EDDCParameters {
    pub fn multiunit_match_score(
        &self,
        unit_id: usize,
        seq_id: usize,
        start_pos: usize,
        end_pos: usize,
    ) -> f64 {
        self.multiunit_match_score[unit_id][seq_id][start_pos][end_pos]
    }
}

#[allow(dead_code)]
fn alignment_dp(s: &Vec<usize>, t: &Vec<usize>, params: &EDDCParameters) -> f64 {
    let n = s.len();
    let m = t.len();
    let mut dp: Vec<Vec<Option<f64>>> = vec![vec![None; m + 1]; n + 1];

    dp[0][0] = Some(0.0);
    for i in 0..n + 1 {
        for j in 0..m + 1 {
            if i + 1 <= n {
                let tmp_val = dp[i][j].unwrap() + params.del_score[s[i]];
                dp[i + 1][j] = match dp[i + 1][j] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
            if j + 1 <= m {
                let tmp_val = dp[i][j].unwrap() + params.ins_score[t[j]];
                dp[i][j + 1] = match dp[i][j + 1] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
            if i + 1 <= n && j + 1 <= m {
                let tmp_val = dp[i][j].unwrap() + params.match_score[s[i]][t[j]];
                dp[i + 1][j + 1] = match dp[i + 1][j + 1] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
        }
    }

    dp[n][m].unwrap()
}

#[allow(dead_code)]
fn alignment_dp_multiunit(
    s_idx: usize,
    s: &Vec<usize>,
    t_idx: usize,
    t: &Vec<usize>,
    units: &Vec<&[u8]>,
    params: &EDDCParameters,
) -> f64 {
    let n = s.len();
    let m = t.len();

    let mut dp: Vec<Vec<Option<f64>>> = vec![vec![None; m + 1]; n + 1];
    dp[0][0] = Some(0.0);
    for i in 0..n {
        dp[i + 1][0] = Some(dp[i][0].unwrap() + params.ins_score[s[i]]);
    }
    for j in 0..m {
        dp[0][j + 1] = Some(dp[0][j].unwrap() + params.ins_score[t[j]]);
    }
    for i in 1..n + 1 {
        for j in 1..m + 1 {
            {
                let tmp_val = dp[i - 1][j].unwrap() + params.del_score[s[i - 1]];
                dp[i][j] = match dp[i][j] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
            {
                let tmp_val = dp[i][j - 1].unwrap() + params.ins_score[t[j - 1]];
                dp[i][j] = match dp[i][j] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
            // s[i] vs t[k..j]
            {
                let mut k = j;
                let mut len = 0;
                while k > 0 && len < 2 * units[s[i - 1]].len() {
                    k -= 1;
                    len += units[t[k]].len();
                    let tmp_val =
                        dp[i - 1][k].unwrap() + params.multiunit_match_score(s[i - 1], t_idx, k, j);
                    dp[i][j] = match dp[i][j] {
                        Some(val) => Some(val.min(tmp_val)),
                        None => Some(tmp_val),
                    }
                }
            }
            // s[k..i] vs t[j]
            {
                let mut k = i;
                let mut len = 0;
                while k > 0 && len < 2 * units[t[j - 1]].len() {
                    k -= 1;
                    len += units[s[k]].len();
                    let tmp_val =
                        dp[k][j - 1].unwrap() + params.multiunit_match_score(t[j - 1], s_idx, k, i);
                    dp[i][j] = match dp[i][j] {
                        Some(val) => Some(val.min(tmp_val)),
                        None => Some(tmp_val),
                    }
                }
            }
        }
    }
    // println!("{:?}", dp);

    dp[n][m].unwrap()
}

fn unit_match_score(u: &[u8], v: &[u8]) -> f64 {
    let n = u.len();
    let m = v.len();
    let mut dp: Vec<Vec<Option<f64>>> = vec![vec![None; m + 1]; n + 1];

    dp[0][0] = Some(0.0);
    for i in 0..n + 1 {
        for j in 0..m + 1 {
            if i + 1 <= n {
                let tmp_val = dp[i][j].unwrap() + 1.0;
                dp[i + 1][j] = match dp[i + 1][j] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
            if j + 1 <= m {
                let tmp_val = dp[i][j].unwrap() + 1.0;
                dp[i][j + 1] = match dp[i][j + 1] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
            if i + 1 <= n && j + 1 <= m {
                let tmp_val = dp[i][j].unwrap() + (if u[i] != v[j] { 1.0 } else { 0.0 });
                dp[i + 1][j + 1] = match dp[i + 1][j + 1] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
        }
    }

    dp[n][m].unwrap()
}

fn unit_ins_score(unit: &[u8]) -> f64 {
    unit.len() as f64 * 0.5
}

fn unit_del_score(unit: &[u8]) -> f64 {
    unit.len() as f64 * 0.5
}

fn calc_multiunit_match_score(u: &[u8], s: &[usize], units: &Vec<&[u8]>) -> Vec<Vec<f64>> {
    let n = s.len();
    let m = u.len();
    let mut res = vec![vec![0.0; n + 1]; n + 1];
    for start_pos in 0..n + 1 {
        let mut dp = vec![None; m + 1];
        dp[0] = Some(0.0);
        for read_unit_pos in start_pos..n {
            let read_unit_len = units[s[read_unit_pos]].len();
            for i in 0..read_unit_len {
                let mut new_dp: Vec<Option<f64>> = vec![None; m + 1];
                for j in 0..m + 1 {
                    if i + 1 <= read_unit_len {
                        let tmp_val = dp[j].unwrap() + 1.0;
                        new_dp[j] = match new_dp[j] {
                            Some(val) => Some(val.min(tmp_val)),
                            None => Some(tmp_val),
                        }
                    }
                    if j + 1 <= m {
                        let tmp_val = dp[j].unwrap() + 1.0;
                        dp[j + 1] = match dp[j + 1] {
                            Some(val) => Some(val.min(tmp_val)),
                            None => Some(tmp_val),
                        }
                    }
                    if i + 1 <= read_unit_len && j + 1 <= m {
                        let tmp_val = dp[j].unwrap()
                            + (if units[s[read_unit_pos]][i] != u[j] {
                                1.0
                            } else {
                                0.0
                            });
                        new_dp[j + 1] = match new_dp[j + 1] {
                            Some(val) => Some(val.min(tmp_val)),
                            None => Some(tmp_val),
                        }
                    }
                }
                dp = new_dp;
                // println!("{:?}", dp);
            }
            for j in 0..m {
                let tmp_val = dp[j].unwrap() + 1.0;
                dp[j + 1] = match dp[j + 1] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
            // println!("{:?}", dp);
            res[start_pos][read_unit_pos + 1] = dp[m].unwrap();
        }
    }
    // println!("unit: {:?}", u);
    // println!("seq: {:?}", s);
    // println!("score: {:?}", res);
    res
}

fn set_params(sequences: &Vec<Vec<usize>>, units: &Vec<&[u8]>) -> EDDCParameters {
    let unit_num = units.len();

    let mut match_score = vec![vec![0.0; unit_num]; unit_num];
    for i in 0..unit_num {
        for j in i + 1..unit_num {
            let tmp = unit_match_score(&units[i], &units[j]);
            match_score[i][j] = tmp;
            match_score[j][i] = tmp;
        }
    }
    let ins_score = units.iter().map(|unit| unit_ins_score(unit)).collect();
    let del_score = units.iter().map(|unit| unit_del_score(unit)).collect();

    let mut multiunit_match_score: Vec<Vec<Vec<Vec<f64>>>> =
        vec![vec![vec![vec![0.0]]; sequences.len()]; units.len()];

    for unit_id in 0..units.len() {
        for seq_id in 0..sequences.len() {
            multiunit_match_score[unit_id][seq_id] =
                calc_multiunit_match_score(&units[unit_id], &sequences[seq_id], &units);
        }
    }
    // println!("{:?}", multiunit_match_score);

    EDDCParameters {
        match_score,
        ins_score,
        del_score,
        multiunit_match_score,
    }
}

pub fn eddc_heuristic(sequences: &Vec<Vec<usize>>, units: &Vec<&[u8]>) -> Vec<Vec<f64>> {
    let params = set_params(sequences, units);

    let seq_num = sequences.len();
    let mut scores = vec![vec![0.0; seq_num]; seq_num];

    for i in 0..seq_num {
        for j in i + 1..seq_num {
            let score = alignment_dp_multiunit(i, &sequences[i], j, &sequences[j], units, &params);
            scores[i][j] = score;
            scores[j][i] = score;
        }
    }

    scores
}
/// Parallel version of the `eddc_heuristics`. The number of the thread is determined automatically.
/// To handle the # of the threads to be used, see [here](https://docs.rs/rayon/latest/rayon/struct.ThreadPoolBuilder.html#method.build_global).
pub fn eddc_heuristic_parallel(sequences: &Vec<Vec<usize>>, units: &Vec<&[u8]>) -> Vec<Vec<f64>> {
    let params = set_params(sequences, units);
    let seq_num = sequences.len();
    let scores: Vec<_> = sequences
        .iter()
        .enumerate()
        .flat_map(|(i, seq)| {
            sequences
                .iter()
                .enumerate()
                .skip(i + 1)
                .map(|(j, seq_aligned)| {
                    let score = alignment_dp_multiunit(i, seq, j, seq_aligned, units, &params);
                    (i, j, score)
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let mut score_matrix = vec![vec![0.0; seq_num]; seq_num];
    for (i, j, score) in scores {
        score_matrix[i][j] = score;
        score_matrix[j][i] = score;
    }
    score_matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn eddc_heuristic_test() {
        let sequences = vec![vec![0, 0, 0, 0, 0], vec![1, 1, 1, 1, 1]];
        let units = vec![b"AT".as_slice(), b"AAT".as_slice()];
        let correct_scores = vec![vec![0.0, 5.0], vec![5.0, 0.0]];
        let scores = eddc_heuristic(&sequences, &units);
        assert_eq!(correct_scores, scores);
    }
    #[test]
    fn eddc_heuristic_test2() {
        let sequences = vec![vec![0], vec![1, 1, 1]];
        let units = vec![b"ATTATTATC".as_slice(), b"ATT".as_slice()];
        let correct_scores = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let scores = eddc_heuristic(&sequences, &units);
        assert_eq!(correct_scores, scores);
    }
}
