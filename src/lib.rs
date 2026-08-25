use rust_htslib::bam;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::error::Error;
use flate2::read::MultiGzDecoder;
use std::fs::File;

#[derive(Clone)]
pub struct RefineInfo {
    pub reference_name: String,
    pub ref_start: i64,
    pub ref_end: i64,
    pub read_start: u32,
    pub read_end: u32,
    pub read_strand: String,
    pub is_secondary: usize,
    pub is_supplementary: usize,
    pub kmer_cnt: usize,
    pub ref_kmer_cnt: usize,
}

#[derive(Clone)]
pub struct TempRecord {
    pub reference_name: String,
    pub ref_start: i64,
    pub ref_end: i64,
    pub read_start: u32,
    pub read_end: u32,
    pub read_strand: String,
    pub is_secondary: usize,
    pub is_supplementary: usize,
    pub kmer_cnt: usize,
    pub ref_kmer_cnt: usize,
    pub flag: usize,
}
impl TempRecord {
    pub fn new() -> Self {
        TempRecord {
            reference_name: String::new(),
            ref_start: 0,
            ref_end: 0,
            read_start: 0,
            read_end: 0,
            read_strand: String::new(),
            is_secondary: 0,
            is_supplementary: 0,
            kmer_cnt: 0,
            ref_kmer_cnt: 0,
            flag: 0,
        }
    }
}

#[derive(Clone)]
pub struct NewRecord {
    pub reference_name: String,
    pub ref_start: i64,
    pub ref_end: i64,
    pub read_start: u32,
    pub read_end: u32,
    pub read_strand: String,
    pub is_secondary: usize,
    pub is_supplementary: usize,
    pub kmer_cnt: usize,
    pub ref_kmer_cnt: usize,
    pub flag: usize,
    pub kmers_list: Vec<usize>,
}

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

pub fn get_deletion_ref_pos(cigartuples: &Vec<(usize, u32)>,
    ref_start: i64,
) -> Vec<(u32, u32)> {
    let mut ref_length: u32 = 0;
    let mut del_ref_db: Vec<(u32, u32)> = Vec::new();

    for (op, len) in cigartuples.iter() {
        match op {
            0 | 7 | 8 => {
                // read_length += len;
                ref_length += len;
            },
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
            },
            _ => (),
        }
    }
    del_ref_db
}

pub fn get_current_ref_pos(
    cigartuples: &Vec<(usize, u32)>,
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

/// The haplotype-specific k-mers of one region, grouped into blocks so that one
/// distinguishing base counts once instead of up to k times.
///
/// A block is a run of consecutive reference start positions, cut by a gap > 1
/// or after `kmer_size` k-mers (a longer run must span more than one locus).
/// Blocks come from the reference k-mer set, not from what a read matched, so
/// counts stay comparable between competing placements.
pub struct KmerBlocks {
    ids: HashMap<u32, usize>,
}

impl KmerBlocks {
    /// Build the blocks from the reference start positions of the
    /// haplotype-specific k-mers found in a region (order does not matter).
    pub fn new(starts: &[u32], kmer_size: u32) -> Self {
        let mut sorted: Vec<u32> = starts.to_vec();
        sorted.sort_unstable();
        sorted.dedup();

        let mut ids: HashMap<u32, usize> = HashMap::with_capacity(sorted.len());
        let mut n_blocks: usize = 0;
        let mut block_len: u32 = 0;
        let mut prev: Option<u32> = None;

        for &s in sorted.iter() {
            let extends_block = match prev {
                Some(p) => s == p + 1 && block_len < kmer_size,
                None => false,
            };
            if !extends_block {
                n_blocks += 1;
                block_len = 0;
            }
            block_len += 1;
            ids.insert(s, n_blocks - 1);
            prev = Some(s);
        }
        KmerBlocks { ids }
    }

    /// Number of distinct blocks (~ distinguishing loci) covered by `starts`.
    /// Positions that are not part of the reference set are ignored.
    pub fn count_hits(&self, starts: &[u32]) -> usize {
        let mut hit: HashSet<usize> = HashSet::new();
        for s in starts.iter() {
            if let Some(&id) = self.ids.get(s) {
                hit.insert(id);
            }
        }
        hit.len()
    }
}

/// Whether a winning placement carries enough markers: `min_markers` of them,
/// or every marker its region offered (`kmer_cnt` never exceeds
/// `ref_kmer_cnt`). `min_markers == 1` disables the rule.
pub fn has_enough_markers(kmer_cnt: usize, ref_kmer_cnt: usize, min_markers: usize) -> bool {
    kmer_cnt >= min_markers || kmer_cnt == ref_kmer_cnt
}

/// Whether the best of the placements competing for the same stretch of a read
/// wins clearly enough: `max / (max + second_max) >= threshold`. A tie, or no
/// marker at all, leaves the segment undetermined.
pub fn is_confident_placement(counts: &[usize], threshold: f64) -> bool {
    let mut max: usize = 0;
    let mut second: usize = 0;
    for &c in counts.iter() {
        if c > max {
            second = max;
            max = c;
        } else if c > second {
            second = c;
        }
    }
    if max == 0 || max == second {
        return false;
    }
    max as f64 / (max + second) as f64 >= threshold
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

#[cfg(test)]
mod tests {
    use super::*;

    const K: u32 = 21;

    #[test]
    fn one_run_shorter_than_k_is_one_block() {
        let starts: Vec<u32> = (100..115).collect(); // 15 consecutive k-mers
        let blocks = KmerBlocks::new(&starts, K);
        assert_eq!(blocks.count_hits(&starts), 1);
        // Matching a single k-mer of the run scores the locus once.
        assert_eq!(blocks.count_hits(&[107]), 1);
    }

    #[test]
    fn a_gap_separates_loci() {
        let mut starts: Vec<u32> = (100..110).collect();
        starts.extend(200..205);
        let blocks = KmerBlocks::new(&starts, K);
        assert_eq!(blocks.count_hits(&starts), 2);
        assert_eq!(blocks.count_hits(&[105, 106]), 1);
        assert_eq!(blocks.count_hits(&[105, 201]), 2);
    }

    #[test]
    fn a_run_longer_than_k_is_cut_every_k_kmers() {
        // 45 consecutive k-mers cannot come from one distinguishing base:
        // 21 + 21 + 3 => 3 blocks.
        let starts: Vec<u32> = (0..45).collect();
        let blocks = KmerBlocks::new(&starts, K);
        assert_eq!(blocks.count_hits(&starts), 3);
        assert_eq!(blocks.count_hits(&[0, 20]), 1);
        assert_eq!(blocks.count_hits(&[0, 21]), 2);
        assert_eq!(blocks.count_hits(&[0, 21, 42]), 3);
    }

    #[test]
    fn unknown_and_duplicate_positions_do_not_inflate_the_count() {
        let starts: Vec<u32> = (100..110).collect();
        let blocks = KmerBlocks::new(&starts, K);
        assert_eq!(blocks.count_hits(&[]), 0);
        assert_eq!(blocks.count_hits(&[100, 100, 101]), 1);
        assert_eq!(blocks.count_hits(&[500]), 0);
    }

    #[test]
    fn threshold_of_half_keeps_the_previous_strict_majority_rule() {
        assert!(is_confident_placement(&[5, 4], 0.5));
        assert!(is_confident_placement(&[1], 0.5));
        assert!(!is_confident_placement(&[4, 4], 0.5));
        assert!(!is_confident_placement(&[0, 0], 0.5));
        assert!(!is_confident_placement(&[], 0.5));
    }

    #[test]
    fn the_marker_floor_only_bites_when_markers_were_missed() {
        // Disabled by default.
        assert!(has_enough_markers(1, 50, 1));
        // Matched everything the region offered: kept however few that was.
        assert!(has_enough_markers(1, 1, 3));
        assert!(has_enough_markers(2, 2, 3));
        // Few markers *and* more were available: rejected.
        assert!(!has_enough_markers(1, 2, 3));
        assert!(!has_enough_markers(2, 50, 3));
        // At or above the floor: kept regardless of what was missed.
        assert!(has_enough_markers(3, 50, 3));
    }

    #[test]
    fn a_higher_threshold_rejects_a_thin_margin() {
        // 5 / (5 + 4) = 0.56
        assert!(!is_confident_placement(&[5, 4], 0.8));
        // 9 / (9 + 1) = 0.9
        assert!(is_confident_placement(&[9, 1], 0.8));
        // Only the runner-up matters, not the rest of the field.
        assert!(is_confident_placement(&[9, 1, 1, 1], 0.8));
        assert!(!is_confident_placement(&[9, 5, 1], 0.8));
    }
}
