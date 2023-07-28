use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::tbx;
use rust_htslib::tbx::Read as tbx_read;
use rust_htslib::tpool::Error;
use std::collections::HashMap;
use std::error::Error as stdError;

use marker_filter::convert_u82String;
use marker_filter::get_cigartuples;
use marker_filter::get_current_ref_pos;
use marker_filter::get_read_position;
use marker_filter::reverse_complement;

mod write_bam;

pub fn run(
    input_bam: &str,
    output_bam: &str,
    hap1_tabix: &str,
    hap2_tabix: &str,
    kmer_size: u32,
) -> Result<(), Box<dyn stdError>> {
    let mut alignments = cal_count_marker(input_bam, hap1_tabix, hap2_tabix, kmer_size);
    let filtered_alignments = filter(&mut alignments);
    write_bam::process_write_bam(input_bam, output_bam, &filtered_alignments);
    Ok(())
}

fn cal_count_marker(
    bamfile: &str,
    hap1_tabix: &str,
    hap2_tabix: &str,
    kmer_size: u32,
) -> HashMap<String, Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>> {
    let mut hap1_tbx_reader =
        tbx::Reader::from_path(hap1_tabix).expect(&format!("Could not open {}", hap1_tabix));

    let mut hap2_tbx_reader =
        tbx::Reader::from_path(hap2_tabix).expect(&format!("Could not open {}", hap2_tabix));

    let mut bam = bam::Reader::from_path(bamfile).expect(&format!("Could not open {}", bamfile));

    let header = bam.header().target_names();

    let mut headers: HashMap<u32, String> = HashMap::new();
    for name in header {
        let r_tid = bam.header().tid(name).unwrap();
        let r_string = convert_u82String(name);
        headers.insert(r_tid, r_string);
    }

    let mut alignments: HashMap<
        String,
        Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>,
    > = HashMap::new();

    let mut read_alignments: Vec<bam::record::Record> = Vec::new();
    let mut prev_read_id = String::new();

    let mut line_num = 0;
    for r_record in bam.records() {
        line_num += 1;
        let record = r_record.unwrap();
        let read_id = convert_u82String(record.qname());
        eprintln!("Processing line {}, {}", line_num, &read_id);
        if record.is_unmapped() {
            continue;
        }

        if prev_read_id.len() == 0 {
            prev_read_id = read_id;
            read_alignments.push(record);
            continue;
        }

        if read_id != prev_read_id {
            let t_alignments = process_read_alignments(
                &read_alignments,
                &headers,
                &mut hap1_tbx_reader,
                &mut hap2_tbx_reader,
                kmer_size,
            );
            alignments.insert(prev_read_id, t_alignments);
            read_alignments = Vec::new();
            read_alignments.push(record);
            prev_read_id = read_id;
        } else {
            read_alignments.push(record);
        }
    }
    let t_alignments = process_read_alignments(
        &read_alignments,
        &headers,
        &mut hap1_tbx_reader,
        &mut hap2_tbx_reader,
        kmer_size,
    );
    alignments.insert(prev_read_id, t_alignments);

    for (key, value) in alignments.iter_mut() {
        value.sort_by(|t1, t2| {
            let key_a = (t1.3 as isize, -(t1.4 as isize));
            let key_b = (t2.3 as isize, -(t2.4 as isize));
            key_a.cmp(&key_b)
        });

        eprint!("{}\t", key);
        for (i, tuple) in value.iter().enumerate() {
            if i != value.len() - 1 {
                eprint!(
                    "{},{},{},{},{},{},{},{},{}\t",
                    tuple.0, tuple.1, tuple.2, tuple.3, tuple.4, tuple.5, tuple.6, tuple.7, tuple.8
                );
            } else {
                eprintln!(
                    "{},{},{},{},{},{},{},{},{}",
                    tuple.0, tuple.1, tuple.2, tuple.3, tuple.4, tuple.5, tuple.6, tuple.7, tuple.8
                );
            }
        }
    }
    alignments
}

fn process_read_alignments(
    read_alignments: &Vec<bam::record::Record>,
    headers: &HashMap<u32, String>,
    hap1_tbx_reader: &mut tbx::Reader,
    hap2_tbx_reader: &mut tbx::Reader,
    kmer_size: u32,
) -> Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)> {
    let mut records = read_alignments.clone();
    records.sort_by_key(|t| t.flags());

    let mut counted_alignments: Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)> =
        Vec::new();
    let mut sequence = String::new();
    for (i, record) in records.iter().enumerate() {
        if i == 0 {
            if record.is_reverse() {
                sequence = reverse_complement(&convert_u82String(&record.seq().as_bytes()));
            } else {
                sequence = convert_u82String(&record.seq().as_bytes());
            }
        }
        let reference_id = record.tid() as u32;
        let reference_name = headers.get(&reference_id).unwrap();
        let ref_start = record.pos();
        let ref_end = record.cigar().end_pos();
        let read_strand = if record.is_reverse() { "-" } else { "+" };
        let mut read_seq = String::new();
        if record.is_reverse() {
            if record.is_supplementary() | record.is_secondary() {
                if sequence.len() != 0 {
                    read_seq = reverse_complement(&sequence);
                } else {
                    eprintln!("Bad alignment: BamFile needs to be sorted");
                }
            } else {
                read_seq = reverse_complement(&sequence);
            }
        } else {
            if record.is_supplementary() | record.is_secondary() {
                if sequence.len() != 0 {
                    read_seq = sequence.clone();
                } else {
                    eprintln!("Bad alignment: BamFile needs to be sorted");
                }
            } else {
                read_seq = sequence.clone();
            }
        }

        let is_supp: usize = if record.is_supplementary() { 1 } else { 0 };

        let is_sec: usize = if record.is_secondary() { 1 } else { 0 };

        let cigartuples = get_cigartuples(&record);
        let read_pos: (u32, u32) = get_read_position(&cigartuples);
        let read_start = read_pos.0;
        let read_end = read_pos.1;
        let read_length = read_end - read_start;

        let mut kmer_cnt: usize = 0;

        let delimiter: u8 = 9; // '\t' for ASCII code
        let mut tbx_reader = &mut *hap1_tbx_reader;
        if &reference_name[0..2] == "h1" {
            let tid = match tbx_reader.tid(reference_name) {
                Ok(tid) => tid,
                Err(_) => {
                    let info: (String, i64, i64, u32, u32, String, usize, usize, usize) = (
                        reference_name.to_string(),
                        ref_start,
                        ref_end,
                        read_start,
                        read_end,
                        read_strand.to_string(),
                        is_sec,
                        is_supp,
                        kmer_cnt,
                    );
                    counted_alignments.push(info);
                    continue;
                }
            };
            let result: Result<(), Error> =
                tbx_reader.fetch(tid as u64, ref_start as u64, ref_end as u64);
            match result {
                Ok(_) => {
                    ();
                }
                Err(_) => {
                    let info: (String, i64, i64, u32, u32, String, usize, usize, usize) = (
                        reference_name.to_string(),
                        ref_start,
                        ref_end,
                        read_start,
                        read_end,
                        read_strand.to_string(),
                        is_sec,
                        is_supp,
                        kmer_cnt,
                    );
                    counted_alignments.push(info);
                    continue;
                }
            }
        } else {
            tbx_reader = hap2_tbx_reader;
            let tid = match tbx_reader.tid(reference_name) {
                Ok(tid) => tid,
                Err(_) => {
                    let info: (String, i64, i64, u32, u32, String, usize, usize, usize) = (
                        reference_name.to_string(),
                        ref_start,
                        ref_end,
                        read_start,
                        read_end,
                        read_strand.to_string(),
                        is_sec,
                        is_supp,
                        kmer_cnt,
                    );
                    counted_alignments.push(info);
                    continue;
                }
            };
            let result: Result<(), Error> =
                tbx_reader.fetch(tid as u64, ref_start as u64, ref_end as u64);
            match result {
                Ok(_) => {
                    ();
                }
                Err(_) => {
                    let info: (String, i64, i64, u32, u32, String, usize, usize, usize) = (
                        reference_name.to_string(),
                        ref_start,
                        ref_end,
                        read_start,
                        read_end,
                        read_strand.to_string(),
                        is_sec,
                        is_supp,
                        kmer_cnt,
                    );
                    // println!("{} {}", read_id, kmer_cnt);
                    counted_alignments.push(info);
                    continue;
                }
            }
        }

        let mut tbx_sequences: HashMap<String, (u32, u32)> = HashMap::new();
        for tbx_record in tbx_reader.records() {
            let in_record = tbx_record.unwrap();
            let chunks: Vec<_> = in_record.split(|&x| x == delimiter).collect();
            let start: i64 = convert_u82String(chunks[1]).parse().unwrap();
            let end: i64 = convert_u82String(chunks[2]).parse().unwrap();
            let strand = convert_u82String(chunks[4]);
            let seq = convert_u82String(chunks[5]);
            if start < ref_start {
                continue;
            }
            if end > ref_end {
                continue;
            }
            if read_strand == strand {
                let tmp_start: u32 = start.try_into().unwrap();
                let tmp_end: u32 = end.try_into().unwrap();
                tbx_sequences.insert(seq, (tmp_start, tmp_end));
            }
        }

        if read_length < kmer_size {
            continue;
        }
        let it_start: usize = read_start.try_into().unwrap();
        let it_end: usize = read_end.try_into().unwrap();
        let k: usize = kmer_size.try_into().unwrap();
        for i in it_start..(it_end - k + 1) {
            let slice = if read_strand == "-" {
                reverse_complement(&read_seq[i..(i + k)])
            } else {
                (&read_seq[i..(i + k)]).to_string()
            };
            let slice = (&read_seq[i..(i + k)]).to_string();
            if let Some(value) = tbx_sequences.get(&slice) {
                if get_current_ref_pos(
                    &cigartuples,
                    ref_start,
                    ref_end,
                    i,
                    i + k,
                    read_strand.to_string(),
                ) == *value
                {
                    kmer_cnt += 1;
                }
            }
        }

        let info: (String, i64, i64, u32, u32, String, usize, usize, usize) = (
            reference_name.to_string(),
            ref_start,
            ref_end,
            read_start,
            read_end,
            read_strand.to_string(),
            is_sec,
            is_supp,
            kmer_cnt,
        );
        counted_alignments.push(info);
    }
    counted_alignments
}

fn filter(
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
                        let t_supp_dist = s_diff_start + s_diff_end;
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
                            let id = supp_id as usize;
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
