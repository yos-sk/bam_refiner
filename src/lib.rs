use rust_htslib::bam;
use std::collections::HashMap;

pub fn reverse_complement(sequence: &str) -> String {
    // complement
    let complement = sequence
        .chars()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => c,
        })
        .collect::<String>();

    // reverse
    let rev_comp = complement.chars().rev().collect::<String>();

    rev_comp
}

#[allow(non_snake_case)]
pub fn convert_u82String(query: &[u8]) -> String {
    let mut converted_query = String::new();

    for item in query {
        let ch = *item as char;
        converted_query.push(ch)
    }
    converted_query
}

pub fn get_cigartuples(record: &bam::Record) -> Vec<(usize, u32)> {
    let mut cigartuples: Vec<(usize, u32)> = vec![];

    for op in record.cigar().iter() {
        match op {
            bam::record::Cigar::Match(len) => {
                cigartuples.push((0, *len));
            }
            bam::record::Cigar::Ins(len) => {
                cigartuples.push((1, *len));
            }
            bam::record::Cigar::Del(len) => {
                cigartuples.push((2, *len));
            }
            bam::record::Cigar::RefSkip(len) => {
                cigartuples.push((3, *len));
            }
            bam::record::Cigar::SoftClip(len) => {
                cigartuples.push((4, *len));
            }
            bam::record::Cigar::HardClip(len) => {
                cigartuples.push((5, *len));
            }
            bam::record::Cigar::Pad(len) => {
                cigartuples.push((6, *len));
            }
            bam::record::Cigar::Equal(len) => {
                cigartuples.push((7, *len));
            }
            bam::record::Cigar::Diff(len) => {
                cigartuples.push((8, *len));
            } /*
              bam::record::Cigar::Back(len) => {
                  cigaråtuples.push((9, len));
                  eprintln!("The backward operation exists.");
              }

              _ => {
                  eprintln!("Unepected cigar.");
              },
              */
              //}
        }
    }
    cigartuples
}

pub fn get_read_position(cigartuples: &Vec<(usize, u32)>) -> (u32, u32, u32) {
    let mut read_start: u32 = 0;
    let mut read_end: u32 = 0;
    let mut read_length: u32 = 0;

    for (i, (op, len)) in cigartuples.iter().enumerate() {
        if i == 0 {
            match op {
                4 | 5 => read_start += len,
                _ => (),
            }
        }
        match op {
            0 | 1 | 7 | 8 => read_length += len,
            _ => (),
        }
        if i == cigartuples.len() - 1 {
            read_end = read_start + read_length;
            match op {
                4 | 5 => read_length += len,
                _ => (),
            }
        }
    }
    (read_start, read_end, read_start + read_length)
}

pub fn get_current_ref_pos(
    cigartuples: &Vec<(usize, u32)>,
    ref_start: i64,
    ref_end: i64,
    it_start: usize,
    it_end: usize,
    strand: String,
) -> (u32, u32) {
    let mut out_start: u32 = ref_start as u32;
    let mut out_end: u32 = out_start;
    let r_start: u32 = it_start as u32;
    let r_end: u32 = it_end as u32;

    let mut read_length: u32 = 0;
    let mut ref_length: u32 = 0;

    let mut f_start = false;
    let mut f_end = false;

    for (op, len) in cigartuples.iter() {
        match op {
            0 | 7 | 8 => {
                read_length += len;
                ref_length += len;
            }
            1 | 4 | 5 => {
                read_length += len;
            }
            2 | 3 => {
                ref_length += len;
            }
            _ => (),
        }

        if read_length >= r_start {
            if !f_start {
                if *op == 0 || *op == 7 || *op == 8 {
                    out_start += ref_length - (read_length - r_start);
                    f_start = true;
                } else {
                    return (0, 0);
                }
            }
        }
        if read_length >= r_end {
            if !f_end {
                if *op == 0 || *op == 7 || *op == 8 {
                    out_end += ref_length - (read_length - r_end);
                } else {
                    return (0, 0);
                }
                f_end = true;
            }
        }
    }
    (out_start, out_end)
}

pub fn filter(
    alignments: &mut HashMap<
        String,
        Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>,
    >,
) -> HashMap<
    String,
    Vec<(
        String,
        i64,
        i64,
        u32,
        u32,
        String,
        usize,
        usize,
        usize,
        usize,
    )>,
> {
    let mut new_results: HashMap<
        String,
        Vec<(
            String,
            i64,
            i64,
            u32,
            u32,
            String,
            usize,
            usize,
            usize,
            usize,
        )>,
    > = HashMap::new();
    for (key, value) in alignments.iter_mut() {
        value.sort_by_key(|tuple| tuple.6);
        let mut prim_info: (
            String,
            i64,
            i64,
            u32,
            u32,
            String,
            usize,
            usize,
            usize,
            usize,
        ) = (
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
        );
        let mut supp_info: Vec<(
            String,
            i64,
            i64,
            u32,
            u32,
            String,
            usize,
            usize,
            usize,
            usize,
        )> = vec![];
        for tuple in value {
            if tuple.6 == 1 {
                let mut supp_id: isize = -1;
                let mut supp_dist: isize = -1;
                let mut supp_cnt: isize = -1;
                let mut prim_dist: isize = -1;
                let mut prim_cnt: isize = -1;
                if prim_info.0 != String::default() {
                    let diff_start = prim_info.3 as isize - tuple.3 as isize;
                    let diff_end = prim_info.4 as isize - tuple.4 as isize;
                    prim_dist = diff_start.abs() + diff_end.abs();
                    prim_cnt = prim_info.8 as isize;
                }
                if !supp_info.is_empty() {
                    for (i, t_supp_info) in supp_info.iter().enumerate() {
                        let s_diff_start = t_supp_info.3 as isize - tuple.3 as isize;
                        let s_diff_end = t_supp_info.4 as isize - tuple.4 as isize;
                        let t_supp_dist = s_diff_start.abs() + s_diff_end.abs();
                        if supp_dist == -1 {
                            supp_dist = t_supp_dist;
                            supp_id = i as isize;
                            supp_cnt = t_supp_info.8 as isize;
                        } else if supp_dist > t_supp_dist {
                            supp_dist = t_supp_dist;
                            supp_id = i as isize;
                            supp_cnt = t_supp_info.8 as isize;
                        }
                    }
                    /*
                    if key == "m64288_220501_014302/165611145/ccs" {
                        eprintln!("Secondary dist: {} {}", prim_dist, supp_dist);
                    }*/
                    if prim_dist <= supp_dist {
                        if prim_cnt < tuple.8.try_into().unwrap() {
                            let info: (
                                String,
                                i64,
                                i64,
                                u32,
                                u32,
                                String,
                                usize,
                                usize,
                                usize,
                                usize,
                            ) = (
                                tuple.0.clone(),
                                tuple.1,
                                tuple.2,
                                tuple.3,
                                tuple.4,
                                tuple.5.clone(),
                                0,
                                tuple.7,
                                tuple.8,
                                1,
                            );
                            prim_info = info;
                        } else if prim_cnt == tuple.8.try_into().unwrap() {
                            let info: (
                                String,
                                i64,
                                i64,
                                u32,
                                u32,
                                String,
                                usize,
                                usize,
                                usize,
                                usize,
                            ) = (
                                prim_info.0.clone(),
                                prim_info.1,
                                prim_info.2,
                                prim_info.3,
                                prim_info.4,
                                prim_info.5.clone(),
                                0,
                                prim_info.7,
                                prim_info.8,
                                0,
                            );
                            prim_info = info;
                        }
                    } else {
                        if supp_cnt < tuple.8.try_into().unwrap() {
                            let info: (
                                String,
                                i64,
                                i64,
                                u32,
                                u32,
                                String,
                                usize,
                                usize,
                                usize,
                                usize,
                            ) = (
                                tuple.0.clone(),
                                tuple.1,
                                tuple.2,
                                tuple.3,
                                tuple.4,
                                tuple.5.clone(),
                                0,
                                1,
                                tuple.8,
                                1,
                            );
                            let id = supp_id as usize;
                            supp_info[id] = info;
                        } else if supp_cnt == tuple.8.try_into().unwrap() {
                            let id = supp_id as usize;
                            let info: (
                                String,
                                i64,
                                i64,
                                u32,
                                u32,
                                String,
                                usize,
                                usize,
                                usize,
                                usize,
                            ) = (
                                supp_info[id].0.clone(),
                                supp_info[id].1,
                                supp_info[id].2,
                                supp_info[id].3,
                                supp_info[id].4,
                                supp_info[id].5.clone(),
                                0,
                                1,
                                supp_info[id].8,
                                0,
                            );
                            supp_info[id] = info;
                        }
                    }
                } else {
                    if prim_cnt < tuple.8.try_into().unwrap() {
                        let info: (
                            String,
                            i64,
                            i64,
                            u32,
                            u32,
                            String,
                            usize,
                            usize,
                            usize,
                            usize,
                        ) = (
                            tuple.0.clone(),
                            tuple.1,
                            tuple.2,
                            tuple.3,
                            tuple.4,
                            tuple.5.clone(),
                            0,
                            tuple.7,
                            tuple.8,
                            1,
                        );
                        prim_info = info;
                    }
                }
            } else if tuple.7 == 1 {
                if tuple.8 == 0 {
                    let info: (
                        String,
                        i64,
                        i64,
                        u32,
                        u32,
                        String,
                        usize,
                        usize,
                        usize,
                        usize,
                    ) = (
                        tuple.0.clone(),
                        tuple.1,
                        tuple.2,
                        tuple.3,
                        tuple.4,
                        tuple.5.clone(),
                        0,
                        1,
                        tuple.8,
                        0,
                    );
                    supp_info.push(info);
                } else {
                    let info: (
                        String,
                        i64,
                        i64,
                        u32,
                        u32,
                        String,
                        usize,
                        usize,
                        usize,
                        usize,
                    ) = (
                        tuple.0.clone(),
                        tuple.1,
                        tuple.2,
                        tuple.3,
                        tuple.4,
                        tuple.5.clone(),
                        0,
                        1,
                        tuple.8,
                        1,
                    );
                    supp_info.push(info);
                }
            } else {
                if tuple.8 == 0 {
                    let info: (
                        String,
                        i64,
                        i64,
                        u32,
                        u32,
                        String,
                        usize,
                        usize,
                        usize,
                        usize,
                    ) = (
                        tuple.0.clone(),
                        tuple.1,
                        tuple.2,
                        tuple.3,
                        tuple.4,
                        tuple.5.clone(),
                        0,
                        tuple.7,
                        tuple.8,
                        0,
                    );
                    prim_info = info;
                } else {
                    let info: (
                        String,
                        i64,
                        i64,
                        u32,
                        u32,
                        String,
                        usize,
                        usize,
                        usize,
                        usize,
                    ) = (
                        tuple.0.clone(),
                        tuple.1,
                        tuple.2,
                        tuple.3,
                        tuple.4,
                        tuple.5.clone(),
                        0,
                        tuple.7,
                        tuple.8,
                        1,
                    );
                    prim_info = info;
                }
            }
        }
        supp_info.push(prim_info);
        new_results.insert(key.to_string(), supp_info);
    }

    for (key, value) in new_results.iter_mut() {
        value.sort_by(|t1, t2| {
            let key_a = (t1.3 as isize, -(t1.4 as isize));
            let key_b = (t2.3 as isize, -(t2.4 as isize));
            key_a.cmp(&key_b)
        });

        print!("{}\t", key);
        for (i, tuple) in value.iter().enumerate() {
            if i != value.len() - 1 {
                print!(
                    "{},{},{},{},{},{},{},{},{},{}\t",
                    tuple.0,
                    tuple.1,
                    tuple.2,
                    tuple.3,
                    tuple.4,
                    tuple.5,
                    tuple.6,
                    tuple.7,
                    tuple.8,
                    tuple.9
                );
            } else {
                println!(
                    "{},{},{},{},{},{},{},{},{},{}",
                    tuple.0,
                    tuple.1,
                    tuple.2,
                    tuple.3,
                    tuple.4,
                    tuple.5,
                    tuple.6,
                    tuple.7,
                    tuple.8,
                    tuple.9
                );
            }
        }
    }
    new_results
}