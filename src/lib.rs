use flate2::read::MultiGzDecoder;
use rust_htslib::bam;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Data {
    pub reference_name: String,
    pub reference_start: i64,
    pub reference_end: i64,
    pub read_name: String,
    pub read_start: u32,
    pub read_end: u32,
    pub is_reverse: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
    pub cigar_tuples: Vec<(u8, u32)>,
    pub rk_cnt: usize,
    pub pk_sk_cnt: usize,
    pub pk_sk_vec: Vec<usize>,
}

pub fn reverse_complement(sequence: &Vec<u8>) -> Vec<u8> {
    let mut complement_seq: Vec<u8> = Vec::with_capacity(sequence.len());

    for &base in sequence.iter().rev() {
        let complement_base = match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => base,
        };

        complement_seq.push(complement_base);
    }

    complement_seq
}

pub fn get_cigartuples(record: &bam::Record) -> Vec<(u8, u32)> {
    let mut cigartuples: Vec<(u8, u32)> = vec![];

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
            }
            /*
            bam::record::Cigar::Back(len) => {
                cigaråtuples.push((9, len));
                eprintln!("The backward operation exists.");
            },
            _ => {
                eprintln!("Unepected cigar.");
            },
            */
        }
    }
    cigartuples
}

pub fn get_read_position(cigartuples: &Vec<(u8, u32)>) -> (u32, u32, u32) {
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

pub fn get_deletion_ref_pos(cigartuples: &Vec<(u8, u32)>, ref_start: i64) -> Vec<(u32, u32)> {
    let mut ref_length: u32 = 0;
    let mut del_ref_db: Vec<(u32, u32)> = Vec::new();

    for (op, len) in cigartuples.iter() {
        match op {
            0 | 7 | 8 => {
                // read_length += len;
                ref_length += len;
            }
            /*
            1 | 4 | 5 => {
                read_length += len;
            }
            */
            2 | 3 => {
                let start = ref_start as u32 + ref_length;
                let end = start + len - 1;
                del_ref_db.push((start, end));
                ref_length += len;
            }
            _ => (),
        }
    }
    del_ref_db
}

pub fn get_current_ref_pos(
    cigartuples: &Vec<(u8, u32)>,
    ref_start: i64,
    _ref_end: i64,
    it_start: usize,
    it_end: usize,
    _strand: String,
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

pub fn open_file<P: AsRef<Path>>(p: P) -> Result<Box<dyn BufRead>, Box<dyn Error>> {
    let r = File::open(p.as_ref())?;
    let ext = p.as_ref().extension();

    if ext == Some(std::ffi::OsStr::new("gz")) {
        let gz = MultiGzDecoder::new(r);
        let buf_reader = BufReader::new(gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::new(r);
        Ok(Box::new(buf_reader))
    }
}

pub fn get_read_name_list(read_name_list: &str) -> Result<HashSet<String>, Box<dyn Error>> {
    let reader = open_file(read_name_list).expect(&format!("Could not open {}", read_name_list));
    let mut read_set: HashSet<String> = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        read_set.insert(line);
    }
    Ok(read_set)
}
