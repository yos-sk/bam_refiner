use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::tbx;
use rust_htslib::tbx::Read as tbx_read;
use rust_htslib::tpool::Error;
use std::collections::HashMap;
use std::error::Error as stdError;

use bam_refiner::convert_u82String;
use bam_refiner::get_cigartuples;
use bam_refiner::get_current_ref_pos;
use bam_refiner::get_read_position;
use bam_refiner::reverse_complement;
use bam_refiner::filter;

#[path = "./write_bam.rs"]
mod write_bam;

pub fn run(
    input_bam: &str,
    output_bam: &str,
    ref_tabix: &str,
    kmer_size: u32,
) -> Result<(), Box<dyn stdError>> {
    let mut alignments = cal_count_marker(input_bam, ref_tabix, kmer_size);
    let filtered_alignments = filter(&mut alignments);
    write_bam::process_write_bam(input_bam, output_bam, &filtered_alignments);
    Ok(())
}

fn cal_count_marker(
    bamfile: &str,
    ref_tabix: &str,
    kmer_size: u32,
) -> HashMap<String, Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>> {
    let mut ref_tbx_reader =
        tbx::Reader::from_path(ref_tabix).expect(&format!("Could not open {}", ref_tabix));

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
                &mut ref_tbx_reader,
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
        &mut ref_tbx_reader,
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
    ref_tbx_reader: &mut tbx::Reader,
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
        let read_pos: (u32, u32, u32) = get_read_position(&cigartuples);
        let read_start = read_pos.0;
        let read_end = read_pos.1;
        let r_read_length = read_pos.2;
        let r_read_start = if read_strand == "+" {
            read_start
        } else {
            r_read_length - read_end
        };
        let r_read_end = if read_strand == "+" {
            read_end
        } else {
            r_read_length - read_start
        };
        let read_length = read_end - read_start;
        /*
        if reference_name == "h2tg000046l" {
            if ref_start == 18025761 && ref_end == 18033661 {
                eprintln!("{} {} {}", read_start, read_end, r_read_length);
            }
        }
        */
        let mut kmer_cnt: usize = 0;

        let delimiter: u8 = 9; // '\t' for ASCII code
        let tid = match ref_tbx_reader.tid(reference_name) {
            Ok(tid) => tid,
            Err(_) => {
                let info: (String, i64, i64, u32, u32, String, usize, usize, usize) = (
                    reference_name.to_string(),
                    ref_start,
                    ref_end,
                    r_read_start,
                    r_read_end,
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
            ref_tbx_reader.fetch(tid as u64, ref_start as u64, ref_end as u64);
        match result {
            Ok(_) => {
                ();
            }
            Err(_) => {
                let info: (String, i64, i64, u32, u32, String, usize, usize, usize) = (
                    reference_name.to_string(),
                    ref_start,
                    ref_end,
                    r_read_start,
                    r_read_end,
                    read_strand.to_string(),
                    is_sec,
                    is_supp,
                    kmer_cnt,
                );
                counted_alignments.push(info);
                continue;
            }
        }

        let mut tbx_sequences: HashMap<String, (u32, u32)> = HashMap::new();
        for tbx_record in ref_tbx_reader.records() {
            let in_record = tbx_record.unwrap();
            let chunks: Vec<_> = in_record.split(|&x| x == delimiter).collect();
            //eprintln!("{:?}", in_record);
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
            // let slice = (&read_seq[i..(i + k)]).to_string();
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
            r_read_start,
            r_read_end,
            read_strand.to_string(),
            is_sec,
            is_supp,
            kmer_cnt,
        );
        counted_alignments.push(info);
    }
    counted_alignments
}
