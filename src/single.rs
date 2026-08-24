use rust_htslib::{
    bam,
    bam::Read,
    tbx,
    tbx::Read as tbx_read,
    tpool::Error,
};
use std::collections::HashMap;
use std::error::Error as stdError;
use std::sync::Arc;
use std::thread;

use crossbeam_channel::{bounded, unbounded};

use bam_refiner::{
    convert_u82String,
    get_cigartuples,
    get_current_ref_pos,
    get_read_position,
    has_enough_markers,
    is_confident_placement,
    reverse_complement,
    KmerBlocks,
    NewRecord,
    RefineInfo,
    TempRecord,
    get_deletion_ref_pos,
};

#[path = "./single/write_bam.rs"]
mod write_bam;

pub fn run(
    input_bam: &str,
    output_bam: &str,
    target_tabix: &str,
    kmer_size: u32,
    threads: usize,
    ratio_threshold: f64,
    min_markers: usize,
) -> Result<(), Box<dyn stdError>> {
    if !(0.0..=1.0).contains(&ratio_threshold) {
        return Err(format!("--ratio-threshold must be within [0.0, 1.0], got {}", ratio_threshold).into());
    }
    let mut alignments = cal_count_marker(input_bam, target_tabix, kmer_size, threads);
    let filtered_alignments = filter(&mut alignments, ratio_threshold, min_markers);
    write_bam::process_write_bam(input_bam, output_bam, &filtered_alignments, threads);
    Ok(())
}

fn cal_count_marker(
    bamfile: &str,
    target_tabix: &str,
    kmer_size: u32,
    threads: usize,
) -> HashMap<String, Vec<RefineInfo>> {
    // Build the tid -> reference name map once (shared, read-only).
    let headers: Arc<HashMap<u32, String>> = {
        let bam = bam::Reader::from_path(bamfile).expect(&format!("Could not open {}", bamfile));
        let mut h: HashMap<u32, String> = HashMap::new();
        for name in bam.header().target_names() {
            let r_tid = bam.header().tid(name).unwrap();
            h.insert(r_tid, convert_u82String(name));
        }
        Arc::new(h)
    };

    let n_workers = threads.max(1);

    // work channel: producer -> workers. Bounded to cap the number of heavy
    // raw-record groups held in memory at once (backpressure).
    let (work_tx, work_rx) = bounded::<Vec<bam::record::Record>>(n_workers * 2);
    // result channel: workers -> collector. Carries only lightweight RefineInfo.
    let (res_tx, res_rx) = unbounded::<(String, Vec<RefineInfo>)>();

    // Spawn the worker pool. Each worker owns its own tabix reader, since
    // tbx::Reader holds mutable fetch state and cannot be shared across threads.
    let mut handles = Vec::with_capacity(n_workers);
    for _ in 0..n_workers {
        let work_rx = work_rx.clone();
        let res_tx = res_tx.clone();
        let headers = Arc::clone(&headers);
        let target_tabix = target_tabix.to_string();
        let handle = thread::spawn(move || {
            let mut target_tbx_reader = tbx::Reader::from_path(&target_tabix)
                .expect(&format!("Could not open {}", target_tabix));

            for read_alignments in work_rx.iter() {
                let read_id = convert_u82String(read_alignments[0].qname());
                let t_alignments = process_read_alignments(
                    &read_alignments,
                    &headers,
                    &mut target_tbx_reader,
                    kmer_size,
                );
                res_tx.send((read_id, t_alignments)).expect("result channel closed");
            }
        });
        handles.push(handle);
    }
    // Drop the main thread's extra handles so the channels close once the
    // producer/workers are done.
    drop(work_rx);
    drop(res_tx);

    // Producer: stream the (name-sorted) BAM in this thread and dispatch one
    // group of records per read. Only this thread touches the BAM reader.
    let mut bam = bam::Reader::from_path(bamfile).expect(&format!("Could not open {}", bamfile));
    bam.set_threads(n_workers).expect(&format!("Failure set {} threads", n_workers));

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
            work_tx
                .send(std::mem::take(&mut read_alignments))
                .expect("work channel closed");
            read_alignments.push(record);
            prev_read_id = read_id;
        } else {
            read_alignments.push(record);
        }
    }
    if !read_alignments.is_empty() {
        work_tx.send(read_alignments).expect("work channel closed");
    }
    // Closing the work channel lets the workers finish and drop their result
    // senders, which in turn ends the collection loop below.
    drop(work_tx);

    // Collector: drain results into the map. The result channel is unbounded,
    // so workers never block on send and cannot deadlock against the producer.
    let mut alignments: HashMap<String, Vec<RefineInfo>> = HashMap::new();
    for (read_id, t_alignments) in res_rx.iter() {
        alignments.insert(read_id, t_alignments);
    }
    for handle in handles {
        handle.join().expect("worker thread panicked");
    }

    for (key, value) in alignments.iter_mut() {
        value.sort_by(|t1, t2| {
            let key_a = (t1.read_start as isize, -(t1.read_end as isize));
            let key_b = (t2.read_start as isize, -(t2.read_end as isize));
            key_a.cmp(&key_b)
        });

        eprint!("{}\t", key);
        for (i, info) in value.iter().enumerate() {
            if i != value.len() - 1 {
                eprint!(
                    "{},{},{},{},{},{},{},{},{},{}\t",
                    info.reference_name, info.ref_start, info.ref_end, info.read_start, info.read_end, info.read_strand, info.is_secondary, info.is_supplementary, info.kmer_cnt, info.ref_kmer_cnt,
                );
            } else {
                eprintln!(
                    "{},{},{},{},{},{},{},{},{},{}",
                    info.reference_name, info.ref_start, info.ref_end, info.read_start, info.read_end, info.read_strand, info.is_secondary, info.is_supplementary, info.kmer_cnt, info.ref_kmer_cnt,
                );
            }
        }
    }
    alignments
}


fn process_read_alignments(
    read_alignments: &Vec<bam::record::Record>,
    headers: &HashMap<u32, String>,
    target_tbx_reader: &mut tbx::Reader,
    kmer_size: u32,
) -> Vec<RefineInfo> {
    let mut records = read_alignments.clone();
    records.sort_by_key(|t| t.flags());

    let mut counted_alignments: Vec<RefineInfo> = Vec::new();
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
        // Reference start positions of the haplotype-specific k-mers matched by
        // the read; resolved into blocks at the end so one distinguishing base
        // is counted once.
        let mut matched_starts: Vec<u32> = Vec::new();

        let delimiter: u8 = 9; // '\t' for ASCII code
        let tid = match target_tbx_reader.tid(reference_name) {
            Ok(tid) => tid,
            Err(_) => {
                let info = RefineInfo {
                    reference_name: reference_name.to_string(),
                    ref_start: ref_start,
                    ref_end: ref_end,
                    read_start: r_read_start,
                    read_end: r_read_end,
                    read_strand: read_strand.to_string(),
                    is_secondary: is_sec,
                    is_supplementary: is_supp,
                    kmer_cnt: 0,
                    ref_kmer_cnt: 0,
                };
                counted_alignments.push(info);
                continue;
            }
        };
        let result: Result<(), Error> =
            target_tbx_reader.fetch(tid as u64, ref_start as u64, ref_end as u64);
        match result {
            Ok(_) => {
                ();
            }
            Err(_) => {
                let info = RefineInfo {
                    reference_name: reference_name.to_string(),
                    ref_start: ref_start,
                    ref_end: ref_end,
                    read_start: r_read_start,
                    read_end: r_read_end,
                    read_strand: read_strand.to_string(),
                    is_secondary: is_sec,
                    is_supplementary: is_supp,
                    kmer_cnt: 0,
                    ref_kmer_cnt: 0,
                };
                counted_alignments.push(info);
                continue;
            }
        }

        let mut tbx_sequences: HashMap<String, (u32, u32)> = HashMap::new();
        // Every haplotype-specific k-mer of this region, used to define the
        // blocks, plus the subset the read can actually observe (k-mers spanned
        // by a deletion are unobservable and excluded from ref_kmer_cnt).
        let mut region_starts: Vec<u32> = Vec::new();
        let mut ref_starts: Vec<u32> = Vec::new();
        let del_ref_pos = get_deletion_ref_pos(&cigartuples, ref_start);
        for tbx_record in target_tbx_reader.records() {
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
                region_starts.push(tmp_start);
                let mut cnt_flag = true;
                for del in del_ref_pos.iter() {
                    if del.1 >= tmp_start && del.0 < tmp_end {
                        cnt_flag = false;
                    }
                }
                 if cnt_flag {
                    ref_starts.push(tmp_start);
                }
            }
        }
        let kmer_blocks = KmerBlocks::new(&region_starts, kmer_size);
        let ref_kmer_cnt = kmer_blocks.count_hits(&ref_starts);

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
                    matched_starts.push(value.0);
                }
            }
        }
        let kmer_cnt = kmer_blocks.count_hits(&matched_starts);

        let info = RefineInfo {
            reference_name: reference_name.to_string(),
            ref_start: ref_start,
            ref_end: ref_end,
            read_start: r_read_start,
            read_end: r_read_end,
            read_strand: read_strand.to_string(),
            is_secondary: is_sec,
            is_supplementary: is_supp,
            kmer_cnt: kmer_cnt,
            ref_kmer_cnt: ref_kmer_cnt,
        };
        counted_alignments.push(info);
    }
    counted_alignments
}

fn filter(
    alignments: &mut HashMap<String, Vec<RefineInfo>>,
    ratio_threshold: f64,
    min_markers: usize,
) -> HashMap<String, Vec<NewRecord>>
{
    let mut new_results: HashMap<String, Vec<NewRecord>> = HashMap::new();

    for (key, value) in alignments.iter_mut() {
        value.sort_by_key(|info| info.is_secondary);
        let mut prim_info = TempRecord::new();
        let mut supp_info: Vec<TempRecord> = Vec::new();
        let mut prim_kmer_cnts: Vec<usize> = Vec::new();
        let mut supp_kmer_cnts: Vec<Vec<usize>> = Vec::new();
        for info in value {
            // Secondary alignments
            if info.is_secondary == 1 {
                let mut supp_id: isize = -1;
                let mut supp_dist: isize = -1;
                let mut supp_cnt: isize = -1;
                let mut prim_dist: isize = -1;
                let mut prim_cnt: isize = -1;
                if !prim_info.reference_name.is_empty() {
                    let diff_start = prim_info.read_start as isize - info.read_start as isize;
                    let diff_end = prim_info.read_end as isize - info.read_end as isize;
                    prim_dist = diff_start.abs() + diff_end.abs();
                    prim_cnt = prim_info.kmer_cnt as isize;
                }
                if !supp_info.is_empty() {
                    for (i, t_supp_info) in supp_info.iter().enumerate() {
                        let s_diff_start = t_supp_info.read_start as isize - info.read_start as isize;
                        let s_diff_end = t_supp_info.read_end as isize - info.read_end as isize;
                        let t_supp_dist = s_diff_start.abs() + s_diff_end.abs();
                        if supp_dist == -1 {
                            supp_dist = t_supp_dist;
                            supp_id = i as isize;
                            supp_cnt = t_supp_info.kmer_cnt as isize;
                        } else if supp_dist > t_supp_dist {
                            supp_dist = t_supp_dist;
                            supp_id = i as isize;
                            supp_cnt = t_supp_info.kmer_cnt as isize;
                        }
                    }
                    if prim_dist <= supp_dist {
                        if prim_cnt < info.kmer_cnt.try_into().unwrap() {
                            prim_info =  TempRecord {
                                reference_name: info.reference_name.clone(),
                                ref_start: info.ref_start,
                                ref_end: info.ref_end,
                                read_start: info.read_start,
                                read_end: info.read_end,
                                read_strand: info.read_strand.clone(),
                                is_secondary: 0,
                                is_supplementary: 0,
                                kmer_cnt: info.kmer_cnt,
                                ref_kmer_cnt: info.ref_kmer_cnt,
                                flag: 1,
                            };
                        } else if prim_cnt == info.kmer_cnt.try_into().unwrap() {
                            prim_info = TempRecord {
                                reference_name: prim_info.reference_name.clone(),
                                ref_start: prim_info.ref_start,
                                ref_end: prim_info.ref_end,
                                read_start: prim_info.read_start,
                                read_end: prim_info.read_end,
                                read_strand: prim_info.read_strand.clone(),
                                is_secondary: 0,
                                is_supplementary: 0,
                                kmer_cnt: prim_info.kmer_cnt,
                                ref_kmer_cnt: prim_info.ref_kmer_cnt,
                                flag: 0,
                            };
                        }
                        prim_kmer_cnts.push(info.kmer_cnt);
                    } else {
                        if supp_cnt < info.kmer_cnt.try_into().unwrap() {
                            let id = supp_id as usize;
                            supp_info[id] = TempRecord {
                                reference_name: info.reference_name.clone(),
                                ref_start: info.ref_start,
                                ref_end: info.ref_end,
                                read_start: info.read_start,
                                read_end: info.read_end,
                                read_strand: info.read_strand.clone(),
                                is_secondary: 0,
                                is_supplementary: 1,
                                kmer_cnt: info.kmer_cnt,
                                ref_kmer_cnt: info.ref_kmer_cnt,
                                flag: 1,
                            };
                        } else if supp_cnt == info.kmer_cnt.try_into().unwrap() {
                            let id = supp_id as usize;
                            supp_info[id] = TempRecord {
                                reference_name: supp_info[id].reference_name.clone(),
                                ref_start: supp_info[id].ref_start,
                                ref_end: supp_info[id].ref_end,
                                read_start: supp_info[id].read_start,
                                read_end: supp_info[id].read_end,
                                read_strand: supp_info[id].read_strand.clone(),
                                is_secondary: 0,
                                is_supplementary: 1,
                                kmer_cnt: supp_info[id].kmer_cnt,
                                ref_kmer_cnt: supp_info[id].ref_kmer_cnt,
                                flag: 0,
                            };
                        }
                        let id = supp_id as usize;
                        supp_kmer_cnts[id].push(info.kmer_cnt)
                    }
                } else {
                    if prim_cnt < info.kmer_cnt.try_into().unwrap() {
                        prim_info = TempRecord {
                            reference_name: info.reference_name.clone(),
                            ref_start: info.ref_start,
                            ref_end: info.ref_end,
                            read_start: info.read_start,
                            read_end: info.read_end,
                            read_strand: info.read_strand.clone(),
                            is_secondary: 0,
                            is_supplementary: 0,
                            kmer_cnt: info.kmer_cnt,
                            ref_kmer_cnt: info.ref_kmer_cnt,
                            flag: 1,
                        };
                    }
                    prim_kmer_cnts.push(info.kmer_cnt);
                }
            // supplementary alignments
            } else if info.is_supplementary == 1 {
                if info.kmer_cnt == 0 {
                    let t_info = TempRecord {
                        reference_name: info.reference_name.clone(),
                        ref_start: info.ref_start,
                        ref_end: info.ref_end,
                        read_start: info.read_start,
                        read_end: info.read_end,
                        read_strand: info.read_strand.clone(),
                        is_secondary: 0,
                        is_supplementary: 1,
                        kmer_cnt: info.kmer_cnt,
                        ref_kmer_cnt: info.ref_kmer_cnt,
                        flag: 0,
                    };
                    supp_info.push(t_info);
                } else {
                    let t_info =  TempRecord {
                        reference_name: info.reference_name.clone(),
                        ref_start: info.ref_start,
                        ref_end: info.ref_end,
                        read_start: info.read_start,
                        read_end: info.read_end,
                        read_strand: info.read_strand.clone(),
                        is_secondary: 0,
                        is_supplementary: 1,
                        kmer_cnt: info.kmer_cnt,
                        ref_kmer_cnt: info.ref_kmer_cnt,
                        flag: 1,
                    };
                    supp_info.push(t_info);
                }
                supp_kmer_cnts.push(vec![info.kmer_cnt]);
            // primary alignment
            } else {
                if info.kmer_cnt  == 0 {
                    prim_info = TempRecord {
                        reference_name: info.reference_name.clone(),
                        ref_start: info.ref_start,
                        ref_end: info.ref_end,
                        read_start: info.read_start,
                        read_end: info.read_end,
                        read_strand: info.read_strand.clone(),
                        is_secondary: 0,
                        is_supplementary: 0,
                        kmer_cnt: info.kmer_cnt,
                        ref_kmer_cnt: info.ref_kmer_cnt,
                        flag: 0,
                    };
                } else {
                    prim_info = TempRecord {
                        reference_name: info.reference_name.clone(),
                        ref_start: info.ref_start,
                        ref_end: info.ref_end,
                        read_start: info.read_start,
                        read_end: info.read_end,
                        read_strand: info.read_strand.clone(),
                        is_secondary: 0,
                        is_supplementary: 0,
                        kmer_cnt: info.kmer_cnt,
                        ref_kmer_cnt: info.ref_kmer_cnt,
                        flag: 1,
                    };
                }
                prim_kmer_cnts.push(info.kmer_cnt);
            }
        }
        eprintln!("Primary kmer cnts\t{}: {:?}", key, prim_kmer_cnts);
        eprintln!("Supplemntary kmer cnts\t{}: {:?}", key, supp_kmer_cnts);

        // The winner of each segment is the placement with the highest count,
        // already picked above. Whether that win is decisive can only be judged
        // once every competitor is known, so the flag is (re)computed here from
        // the collected counts instead of incrementally in the loop.
        prim_info.flag = if is_confident_placement(&prim_kmer_cnts, ratio_threshold)
            && has_enough_markers(prim_info.kmer_cnt, prim_info.ref_kmer_cnt, min_markers)
        { 1 } else { 0 };
        for (i, t_supp_info) in supp_info.iter_mut().enumerate() {
            t_supp_info.flag = if is_confident_placement(&supp_kmer_cnts[i], ratio_threshold)
                && has_enough_markers(t_supp_info.kmer_cnt, t_supp_info.ref_kmer_cnt, min_markers)
            { 1 } else { 0 };
        }

        let f_prim_info = NewRecord {
            reference_name: prim_info.reference_name.clone(),
            ref_start: prim_info.ref_start,
            ref_end: prim_info.ref_end,
            read_start: prim_info.read_start,
            read_end: prim_info.read_end,
            read_strand: prim_info.read_strand.clone(),
            is_secondary: prim_info.is_secondary,
            is_supplementary: prim_info.is_supplementary,
            kmer_cnt: prim_info.kmer_cnt,
            ref_kmer_cnt: prim_info.ref_kmer_cnt,
            flag: prim_info.flag,
            kmers_list: prim_kmer_cnts,
        };

        let mut new_result: Vec<NewRecord> = vec![f_prim_info];

        for (i, t_supp_info) in supp_info.iter().enumerate() {
            let f_supp_info = NewRecord {
                reference_name: t_supp_info.reference_name.clone(),
                ref_start: t_supp_info.ref_start,
                ref_end: t_supp_info.ref_end,
                read_start: t_supp_info.read_start,
                read_end: t_supp_info.read_end,
                read_strand: t_supp_info.read_strand.clone(),
                is_secondary: t_supp_info.is_secondary,
                is_supplementary: t_supp_info.is_supplementary,
                kmer_cnt: t_supp_info.kmer_cnt,
                ref_kmer_cnt: t_supp_info.ref_kmer_cnt,
                flag: t_supp_info.flag,
                kmers_list: supp_kmer_cnts[i].clone(),
            };
            new_result.push(f_supp_info);
        }
        new_results.insert(key.to_string(), new_result);
    }

    for (key, value) in new_results.iter_mut() {
        value.sort_by(|t1, t2| {
            let key_a = (t1.read_start as isize, -(t1.read_end as isize));
            let key_b = (t2.read_start as isize, -(t2.read_end as isize));
            key_a.cmp(&key_b)
        });

        print!("{}\t", key);
        for (i, info) in value.iter().enumerate() {
            if i != value.len() - 1 {
                print!(
                    "{},{},{},{},{},{},{},{},{},{},{}\t",
                    &info.reference_name,
                    info.ref_start,
                    info.ref_end,
                    info.read_start,
                    info.read_end,
                    &info.read_strand,
                    info.is_secondary,
                    info.is_supplementary,
                    info.kmer_cnt,
                    info.ref_kmer_cnt,
                    info.flag,
                );
            } else {
                println!(
                    "{},{},{},{},{},{},{},{},{},{},{}",
                    &info.reference_name,
                    info.ref_start,
                    info.ref_end,
                    info.read_start,
                    info.read_end,
                    &info.read_strand,
                    info.is_secondary,
                    info.is_supplementary,
                    info.kmer_cnt,
                    info.ref_kmer_cnt,
                    info.flag,
                );
            }
        }
    }
    new_results
}
