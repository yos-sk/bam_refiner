use rust_htslib::{tbx, tbx::Read};
use rust_htslib::tpool::Error;
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error as stdError;
use std::io::BufRead;
use std::sync::{Arc, Mutex};
use std::thread;

use bam_refiner::get_deletion_ref_pos;
use bam_refiner::get_current_ref_pos;
use bam_refiner::reverse_complement;
use bam_refiner::open_file;

use bam_refiner::Data;

#[path = "./collect_alignments.rs"]
mod collect_alignments;

pub fn run(
    alignments: Vec<Data>,
    sequences: HashMap<String, Vec<u8>>,
    hap1_tabix: &str,
    hap2_tabix: &str,
    hap1_list: &str,
    hap2_list: &str,
    kmer_size: u32,
    threads: usize,
) -> Result<(), Box<dyn stdError>> {
    let shared_alignments: Arc<Mutex<Vec<Data>>> = Arc::new(Mutex::new(alignments));
    let shared_sequences: Arc<Mutex<HashMap<String, Vec<u8>>>> = Arc::new(Mutex::new(sequences));
    let shared_hap1_tabix = Arc::new(hap1_tabix.to_string());
    let shared_hap2_tabix = Arc::new(hap2_tabix.to_string());
    let shared_hap1_list = Arc::new(hap1_list.to_string());
    let shared_hap2_list = Arc::new(hap2_list.to_string());

    let threads: Vec<_> = (0..threads)
        .map(|i| {
            let shared_map = Arc::clone(&shared_alignments);
            let shared_seq = Arc::clone(&shared_sequences);
            let shared_h1_tbx = Arc::clone(&shared_hap1_tabix);
            let shared_h2_tbx = Arc::clone(&shared_hap2_tabix);
            let shared_h1_ls = Arc::clone(&shared_hap1_list);
            let shared_h2_ls = Arc::clone(&shared_hap2_list);

            thread::spawn(move || {
                // process by thread
                let mut map = shared_map.lock().unwrap();
                let seq = shared_seq.lock().unwrap();
                count_kmers(&mut map, &seq, i, &shared_h1_tbx, &shared_h2_tbx, &shared_h1_ls, &shared_h2_ls, threads, kmer_size);  
            })
        })
        .collect();
    
    for handle in threads {
        handle.join().unwrap();
    }

    let map = shared_alignments.lock().unwrap();
    println!("{:?}", *map);
    Ok(())
}

fn count_kmers(alignments: &mut Vec<Data>, sequences: &HashMap<String, Vec<u8>>, index:usize, hap1_tabix: &str, hap2_tabix: &str, hap1_list: &str, hap2_list: &str, threads: usize, kmer_size: u32) {
    let mut hap1_tbx_reader =
        tbx::Reader::from_path(hap1_tabix).expect(&format!("Could not open {}", hap1_tabix));
    let mut hap2_tbx_reader =
        tbx::Reader::from_path(hap2_tabix).expect(&format!("Could not open {}", hap2_tabix));
    let hap1_set: HashSet<String> = get_read_name_list(hap1_list).expect(&format!("Could not read {}", hap1_list));
    let hap2_set: HashSet<String> = get_read_name_list(hap2_list).expect(&format!("Could not read {}", hap2_list));

    let start = alignments.len() / threads * index;
    let end = if index != threads - 1 {
        alignments.len() / threads * (index + 1)
    } else {
        alignments.len()
    };
    let mut alignments_thread = (&alignments[start..end]).to_vec();

    for read in alignments_thread.iter_mut() {
        if let Some(value) = sequences.get(&read.read_name) {
            let (ref_kmer_cnt, read_kmer_cnt) = count_kmers_tbx(read, value, &mut hap1_tbx_reader, &mut hap2_tbx_reader, &hap1_set, &hap2_set, kmer_size);
            eprintln!("{} {} {}", read.read_name, ref_kmer_cnt, read_kmer_cnt);
            read.rk_cnt = ref_kmer_cnt;
            read.pk_sk_cnt = read_kmer_cnt;
        } else {
            continue;
        }        
    }
}

fn count_kmers_tbx(read: &mut Data, read_seq: &Vec<u8>, hap1_tbx_reader: &mut tbx::Reader, hap2_tbx_reader: &mut tbx::Reader, hap1_set: &HashSet<String>, hap2_set: &HashSet<String>, kmer_size: u32) -> (usize, usize) {
    let tbx_reader = if hap1_set.contains(&read.reference_name) {
        &mut *hap1_tbx_reader
    } else if hap2_set.contains(&read.reference_name) {
        &mut *hap2_tbx_reader
    } else {
        return (0, 0);
    };

    let tid = match tbx_reader.tid(&read.reference_name) {
        Ok(tid) => tid,
        Err(_) => {
            return (0, 0);
        },
    };

    let result: Result<(), Error> = tbx_reader.fetch(tid as u64, read.reference_start as u64, read.reference_end as u64);
    match result {
        Ok(_) => (),
        Err(_) => {
            return (0, 0);
        },
    }

    let mut tbx_sequences: HashMap<String, (u32, u32)> = HashMap::new();
    let del_ref_pos = get_deletion_ref_pos(&read.cigar_tuples, read.reference_start);
    let delimiter: u8 = 9;
    let read_strand = if read.is_reverse {
        "-"
    } else {
        "+"
    };
    let mut ref_kmer_cnt = 0;
    let mut read_kmer_cnt = 0;
    for tbx_record in tbx_reader.records() {
        let in_record = tbx_record.unwrap();
        let chunks: Vec<_> = in_record.split(|&x| x == delimiter).collect();
        let start: i64 = String::from_utf8_lossy(chunks[1]).to_string().parse().unwrap();
        let end: i64 = String::from_utf8_lossy(chunks[2]).to_string().parse().unwrap();
        let strand = String::from_utf8_lossy(chunks[4]).to_string();
        let kmer_seq = String::from_utf8_lossy(chunks[5]).to_string();
        if start < read.reference_start {
            continue;
        }
        if end > read.reference_end {
            continue;
        }
        if read_strand == strand {
            let tmp_start: u32 = start.try_into().unwrap();
            let tmp_end: u32 = end.try_into().unwrap();
            tbx_sequences.insert(kmer_seq, (tmp_start, tmp_end));
            let mut cnt_flag = true;
            for del in del_ref_pos.iter() {
                if del.1 >= tmp_start && del.0 <= tmp_end {
                    cnt_flag = false;
                }
            }
            if cnt_flag {
                ref_kmer_cnt += 1;
            }
        }
    }

    let it_start: usize = read.read_start as usize;
    let it_end: usize = read.read_end as usize;
    let k: usize = kmer_size as usize;
    for i in it_start..(it_end - k + 1) {
        let slice = if read_strand == "-" {
            // eprintln!("{} {} {} {} {:?} {:?}", i, i + k, read_seq.len(), read.read_name, read.is_secondary, read.is_supplementary);
            String::from_utf8_lossy(&reverse_complement(&read_seq[i..(i + k)].to_vec())).to_string()
        } else {
            // eprintln!("{} {} {} {} {:?} {:?}", i, i + k, read_seq.len(), read.read_name, read.is_secondary, read.is_supplementary);
            String::from_utf8_lossy(&read_seq[i..(i + k)].to_vec()).to_string()
        };
        if let Some(value) = tbx_sequences.get(&slice) {
            if get_current_ref_pos(
                &read.cigar_tuples,
                read.reference_start,
                read.reference_end,
                i,
                i + k,
                read_strand.to_string(),
            ) == *value
            {
                read_kmer_cnt += 1;
            }
        }
    }
    (ref_kmer_cnt, read_kmer_cnt)
}

fn get_read_name_list(read_name_list: &str) -> Result<HashSet<String>, Box<dyn stdError>> {
    let reader = open_file(read_name_list).expect(&format!("Could not open {}", read_name_list));
    let mut read_set: HashSet<String> = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        read_set.insert(line);
    }
    Ok(read_set)
}