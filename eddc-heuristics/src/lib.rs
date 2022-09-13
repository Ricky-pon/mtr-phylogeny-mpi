pub struct EDDCParameters {
    match_score: Vec<Vec<f64>>,
    // Insertion to the *unit*. Should be negative.
    ins_score: Vec<f64>,
    // Deletion from the *unit*. Should be negative.
    del_score: Vec<f64>,
}

fn f64_min(a: f64, b: f64) -> f64 {
    match a <= b {
        true => a,
        _ => b,
    }
}

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
                    Some(val) => Some(f64_min(val, tmp_val)),
                    None => Some(tmp_val),
                }
            }
            if j + 1 <= m {
                let tmp_val = dp[i][j].unwrap() + params.ins_score[t[j]];
                dp[i][j + 1] = match dp[i][j + 1] {
                    Some(val) => Some(f64_min(val, tmp_val)),
                    None => Some(tmp_val),
                }
            }
            if i + 1 <= n && j + 1 <= m {
                let tmp_val = dp[i][j].unwrap() + params.match_score[s[i]][t[j]];
                dp[i + 1][j + 1] = match dp[i + 1][j + 1] {
                    Some(val) => Some(f64_min(val, tmp_val)),
                    None => Some(tmp_val),
                }
            }
        }
    }

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
                    Some(val) => Some(f64_min(val, tmp_val)),
                    None => Some(tmp_val),
                }
            }
            if j + 1 <= m {
                let tmp_val = dp[i][j].unwrap() + 1.0;
                dp[i][j + 1] = match dp[i][j + 1] {
                    Some(val) => Some(f64_min(val, tmp_val)),
                    None => Some(tmp_val),
                }
            }
            if i + 1 <= n && j + 1 <= m {
                let tmp_val = dp[i][j].unwrap() + (if u[i] != v[j] { 1.0 } else { 0.0 });
                dp[i + 1][j + 1] = match dp[i + 1][j + 1] {
                    Some(val) => Some(f64_min(val, tmp_val)),
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

fn set_params(units: &Vec<&[u8]>) -> EDDCParameters {
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
    // println!("{:?}", match_score);

    EDDCParameters {
        match_score,
        ins_score,
        del_score,
    }
}

pub fn eddc_heuristic(sequences: &Vec<Vec<usize>>, units: &Vec<&[u8]>) -> Vec<Vec<f64>> {
    let params = set_params(units);

    let seq_num = sequences.len();
    let mut scores = vec![vec![0.0; seq_num]; seq_num];

    for i in 0..seq_num {
        for j in i + 1..seq_num {
            let score = alignment_dp(&sequences[i], &sequences[j], &params);
            scores[i][j] = score;
            scores[j][i] = score;
        }
    }

    scores
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
}
