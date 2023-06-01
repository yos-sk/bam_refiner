#[allow(unused_imports)]
use std::path::Path;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::tbx;
use rust_htslib::tbx::Read as tbx_read;
use rust_htslib::tpool::Error;
use std::env;
use std::collections::HashMap;


#[allow(unused)]
fn main() {
    let args: Vec<String> = env::args().collect();
    let bam_file = &args[1];
    let hap1_tabix = &args[2];
    let hap2_tabix = &args[3];
    let out_bam = &args[4];
    let kmer_size: u32 = 21;
    //cal_count_marker(bam_file);
    let mut alignments = cal_count_marker(bam_file, hap1_tabix, hap2_tabix, kmer_size);
    let mut filtered_alignments = filter(&mut alignments);
    write_bam(bam_file, out_bam, &mut filtered_alignments);
}

#[allow(non_snake_case)]
fn convert_u82String(query: &[u8]) -> String {
    let mut converted_query = String::new();

    for item in query {
        let ch = *item as char; 
        converted_query.push(ch)
    }
    converted_query
}

fn get_cigartuples(record: &bam::Record) -> Vec<(usize, u32)> {
    let mut cigartuples: Vec<(usize, u32)> = vec![]; 
    if record.is_reverse() {
        for op in record.cigar().iter().rev() {
            match op {
                bam::record::Cigar::Match(len) => {
                    cigartuples.push((0, *len));
                },
                bam::record::Cigar::Ins(len) => {
                    cigartuples.push((1, *len));
                },
                bam::record::Cigar::Del(len) => {
                    cigartuples.push((2, *len));
                },
                bam::record::Cigar::RefSkip(len) => {
                    cigartuples.push((3, *len));
                },
                bam::record::Cigar::SoftClip(len) => {
                    cigartuples.push((4, *len));
                },
                bam::record::Cigar::HardClip(len) => {
                    cigartuples.push((5, *len));
                },
                bam::record::Cigar::Pad(len) => {
                    cigartuples.push((6, *len));
                },
                bam::record::Cigar::Equal(len) => {
                    cigartuples.push((7, *len));
                },
                bam::record::Cigar::Diff(len) => {
                    cigartuples.push((8, *len));
                },
                /*
                bam::record::Cigar::Back(len) => {
                    cigartuples.push((9, len));
                    eprintln!("The backward operation exists.");
                }
                
                _ => {
                    eprintln!("Unepected cigar.");
                },
                */
            }
        }
    } else {
        for op in record.cigar().iter() {
            match op {
                bam::record::Cigar::Match(len) => {
                    cigartuples.push((0, *len));
                },
                bam::record::Cigar::Ins(len) => {
                    cigartuples.push((1, *len));
                },
                bam::record::Cigar::Del(len) => {
                    cigartuples.push((2, *len));
                },
                bam::record::Cigar::RefSkip(len) => {
                    cigartuples.push((3, *len));
                },
                bam::record::Cigar::SoftClip(len) => {
                    cigartuples.push((4, *len));
                },
                bam::record::Cigar::HardClip(len) => {
                    cigartuples.push((5, *len));
                },
                bam::record::Cigar::Pad(len) => {
                    cigartuples.push((6, *len));
                },
                bam::record::Cigar::Equal(len) => {
                    cigartuples.push((7, *len));
                },
                bam::record::Cigar::Diff(len) => {
                    cigartuples.push((8, *len));
                },
                /*
                bam::record::Cigar::Back(len) => {
                    cigaråtuples.push((9, len));
                    eprintln!("The backward operation exists.");
                }
                
                _ => {
                    eprintln!("Unepected cigar.");
                },
                */
            }
        }
    }
    cigartuples
}


fn get_read_position(cigartuples: &Vec<(usize, u32)>) -> (u32, u32) {
    let mut read_start: u32 = 0;
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
    }
    // eprintln!("Read length: {}", read_length);
    let read_end = read_start + read_length;
    (read_start, read_end)
}

fn get_current_ref_pos(cigartuples: &Vec<(usize, u32)>, ref_start: i64, ref_end: i64, it_start: usize, it_end: usize, strand: String) -> (u32, u32) {
    let mut out_start: u32 = if strand == "+" {
        ref_start as u32
    } else {
        ref_end as u32
    } ;
    let mut out_end: u32 = out_start;
    let r_start: u32 = it_start as u32;
    let r_end: u32 = it_end as u32;

    let mut read_length: u32 = 0;
    let mut ref_length: u32 = 0;

    let mut f_start = false;
    let mut f_end = false;

    for (op, len) in cigartuples.iter() {
        match op {
            0 | 7 | 8  => {
                read_length += len;
                ref_length += len;
            },
            1 | 4 | 5 => {
                read_length += len;
            },
            2 | 3 => {
                ref_length += len;
            }
            _ => (),
        }

        if read_length >= r_start {
            if !f_start {
                if *op == 0 || *op == 7 || *op == 8 {
                    if strand == "-" {
                        out_end -= ref_length - (read_length - r_start);
                    } else {
                        out_start += ref_length - (read_length - r_start);
                    }
                    f_start = true;
                } else {
                    return (0, 0);
                }
            }
        }
        if read_length >= r_end {
            if !f_end {
                if *op == 0 || *op == 7 || *op == 8 {
                    if strand == "-" {
                        out_start -= ref_length - (read_length - r_end);
                    } else {
                        out_end += ref_length - (read_length - r_end);
                    }
                } else {
                    return (0, 0);
                }
                f_end = true;
            }
        }
    }
    (out_start, out_end)
}

#[allow(unused)]
fn cal_count_marker(bamfile: &str, hap1_tabix: &str, hap2_tabix: &str, kmer_size: u32) ->  HashMap<String, Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>> {
    
    let mut hap1_tbx_reader = tbx::Reader::from_path(hap1_tabix)
        .expect(&format!("Could not open {}", hap1_tabix));
    
    let mut hap2_tbx_reader = tbx::Reader::from_path(hap2_tabix)
        .expect(&format!("Could not open {}", hap2_tabix));
    
    let mut bam = bam::Reader::from_path(bamfile)
        .expect(&format!("Could not open {}", bamfile));
    
    let header = bam.header().target_names();
    
    let mut headers: HashMap<u32, String> = HashMap::new();
    for name in header {
        let r_tid = bam.header().tid(name).unwrap();
        let r_string = convert_u82String(name);
        headers.insert(r_tid, r_string);
    }


    let mut sequences: HashMap<String, String> = HashMap::new();
    let mut alignments: HashMap<String, Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>> = HashMap::new();

    let mut read_alignments: Vec<bam::record::Record> = Vec::new();
    let mut prev_read_id = String::new();

    let mut line_num = 0;
    for r_record in bam.records() {
        line_num += 1;
        eprintln!("Processing line {}", line_num);
        let record = r_record.unwrap();
        if record.is_unmapped() {
            continue;
        }
        let read_id = convert_u82String(record.qname());

        if prev_read_id.len() == 0 {
            prev_read_id = read_id;
            read_alignments.push(record);
            continue;
        }

        if read_id != prev_read_id {
            let t_alignments = process_read_alignments(&read_alignments, &headers, &mut hap1_tbx_reader, &mut hap2_tbx_reader, kmer_size);
            alignments.insert(prev_read_id, t_alignments);
            read_alignments = Vec::new();
            read_alignments.push(record);
            prev_read_id = read_id;
        } else {
            read_alignments.push(record);
        }
    }
    let t_alignments = process_read_alignments(&read_alignments, &headers, &mut hap1_tbx_reader, &mut hap2_tbx_reader, kmer_size);
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
                eprint!("{},{},{},{},{},{},{},{},{}\t", tuple.0, tuple.1, tuple.2, tuple.3, tuple.4, tuple.5, tuple.6, tuple.7, tuple.8);
            } else {
                eprintln!("{},{},{},{},{},{},{},{},{}", tuple.0, tuple.1, tuple.2, tuple.3, tuple.4, tuple.5, tuple.6, tuple.7, tuple.8);
            }
        }
    }

    alignments
}

fn process_read_alignments(read_alignments: &Vec<bam::record::Record>, headers:&HashMap<u32, String>, hap1_tbx_reader: &mut tbx::Reader, hap2_tbx_reader: &mut tbx::Reader, kmer_size: u32) 
-> Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>
{
    let mut records = read_alignments.clone();
    records.sort_by_key(|t| {
        t.flags()
    });

    let mut counted_alignments: Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)> = Vec::new();
    let mut sequence = String::new();
    for (i, record) in records.iter().enumerate() {
        if i == 0 {
            sequence = convert_u82String(&record.seq().as_bytes());
            // eprintln!("{}, {}, {}", read_id, record.flags(), sequence.len());
        }
        let reference_id = record.tid() as u32;
        let reference_name = headers.get(&reference_id).unwrap(); 
        let ref_start = record.pos();
        let ref_end = record.cigar().end_pos();
        let read_strand = if record.is_reverse() {
            "-"
        } else {
            "+"
        };
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
        
        let is_supp: usize = if record.is_supplementary() {
            1
        } else {
            0
        };

        let is_sec: usize = if record.is_secondary() {
            1
        } else {
            0
        };
        

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
                    let info: (String, i64, i64, u32, u32, String, usize, usize, usize)
                        = (reference_name.to_string(), ref_start, ref_end, read_start, read_end, read_strand.to_string(), is_sec, is_supp, kmer_cnt);
                    counted_alignments.push(info);
                    continue;
                },
            };
            let result: Result<(), Error> = tbx_reader
                                .fetch(tid as u64, ref_start as u64, ref_end as u64);
            match result {
                Ok(_) => {
                    ();
                },
                Err(_) => {
                    let info: (String, i64, i64, u32, u32, String, usize, usize, usize) 
                        = (reference_name.to_string(), ref_start, ref_end, read_start, read_end, read_strand.to_string(), is_sec, is_supp, kmer_cnt);
                    counted_alignments.push(info);
                    continue;
                }
            }
        } else {
            tbx_reader = hap2_tbx_reader;
            let tid = match tbx_reader.tid(reference_name) {
                Ok(tid) => tid,
                Err(_) => {
                    let info: (String, i64, i64, u32, u32, String, usize, usize, usize) 
                        = (reference_name.to_string(), ref_start, ref_end, read_start, read_end, read_strand.to_string(), is_sec, is_supp, kmer_cnt);
                    counted_alignments.push(info);
                    continue;
                },
            };
            let result: Result<(), Error> = tbx_reader
                                .fetch(tid as u64, ref_start as u64, ref_end as u64);
            match result {
                Ok(_) => {
                    ();
                },
                Err(_) => {
                    let info: (String, i64, i64, u32, u32, String, usize, usize, usize)
                        = (reference_name.to_string(), ref_start, ref_end, read_start, read_end, read_strand.to_string(), is_sec, is_supp, kmer_cnt);
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
            let hap_name = convert_u82String(chunks[3]);
            let strand = convert_u82String(chunks[4]);
            let seq = convert_u82String(chunks[5]);
            // println!("{}\t{}", read_strand, strand);
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
        for i in it_start..(it_end-k+1) {
            let slice = (&read_seq[i..(i+k)]).to_string();
            // eprintln!("{}\t{}\t{:?}", slice, seq, slice==seq);
            if let Some(value) = tbx_sequences.get(&slice) {
                /*
                if read_strand == "-" {
                    println!("{} {} {} {} {} {}", read_id, reference_name, ref_start, ref_end, i, i + k);
                    println!("{:?} {:?}", get_current_ref_pos(&cigartuples, ref_start, ref_end, i, i + k, read_strand.to_string()), value);
                }*/
                if get_current_ref_pos(&cigartuples, ref_start, ref_end, i, i + k, read_strand.to_string()) == *value {
                    kmer_cnt += 1;
                }
            }
        }
        // println!("{} {}", read_id, kmer_cnt);
        /*
        if read_id == "m64288_220429_181717/88016460/ccs" {
            println!("ref position: {}:{}-{}", reference_name, ref_start, ref_end);
            println!("Originel: {:?}", convert_u82String(&record.seq().as_bytes()));
            println!("Sequence: {:?}", read_seq);
            println!("kmer_db: {:?}", tbx_sequences);
        }
        */

        let info: (String, i64, i64, u32, u32, String, usize, usize, usize) 
            = (reference_name.to_string(), ref_start, ref_end, read_start, read_end, read_strand.to_string(), is_sec, is_supp, kmer_cnt);
        counted_alignments.push(info);
    }
    counted_alignments
}

fn reverse_complement(sequence: &str) -> String {
    // complement
    let complement = sequence.chars().map(|c| match c {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        _ => c,
    }).collect::<String>();

    // reverse
    let rev_comp = complement.chars().rev().collect::<String>();

    rev_comp
}

fn filter(alignments: &mut HashMap<String, Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>>) 
-> HashMap<String, Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>>  {
    let mut new_results: HashMap<String, Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>> = HashMap::new();
    for (key, value) in alignments.iter_mut() {
        value.sort_by_key(|tuple| {
            tuple.6
        });
        let mut prim_info: (String, i64, i64, u32, u32, String, usize, usize, usize) 
        = (Default::default(), Default::default(), Default::default(), Default::default(), Default::default(), Default::default(), Default::default(), Default::default(), Default::default());
        let mut supp_info: Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)> = vec![]; 
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
                    if prim_dist < supp_dist {
                        if prim_cnt < tuple.8.try_into().unwrap() {
                            let info: (String, i64, i64, u32, u32, String, usize, usize, usize) = (tuple.0.clone(), tuple.1, tuple.2, tuple.3, tuple.4, tuple.5.clone(), 0, tuple.7, tuple.8);
                            prim_info = info;
                        }
                    } else {
                        if supp_cnt < tuple.8.try_into().unwrap() {
                            let info: (String, i64, i64, u32, u32, String, usize, usize, usize) = (tuple.0.clone(), tuple.1, tuple.2, tuple.3, tuple.4, tuple.5.clone(), 0, 1, tuple.8);
                            let id = supp_id as usize;
                            supp_info[id] = info;
                        }
                    }
                } else {
                    if prim_cnt < tuple.8.try_into().unwrap() {
                        let info: (String, i64, i64, u32, u32, String, usize, usize, usize) = (tuple.0.clone(), tuple.1, tuple.2, tuple.3, tuple.4, tuple.5.clone(), 0, tuple.7, tuple.8);
                        prim_info = info;
                    }
                }
            } else if tuple.7 == 1 {
                supp_info.push(tuple.clone());
            } else {
                prim_info = tuple.clone();
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
                print!("{},{},{},{},{},{},{},{},{}\t", tuple.0, tuple.1, tuple.2, tuple.3, tuple.4, tuple.5, tuple.6, tuple.7, tuple.8);
            } else {
                println!("{},{},{},{},{},{},{},{},{}", tuple.0, tuple.1, tuple.2, tuple.3, tuple.4, tuple.5, tuple.6, tuple.7, tuple.8);
            }
        }
    }
    new_results
}


fn write_bam(bamfile: &str, output_bam: &str, filtered_alignments: &mut HashMap<String, Vec<(String, i64, i64, u32, u32, String, usize, usize, usize)>>) {
    let mut bam = bam::Reader::from_path(bamfile)
        .expect(&format!("Could not open {}", bamfile));
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(output_bam, &header, bam::Format::Bam)
        .expect(&format!("Could not open {}", output_bam));
    
    let mut headers: HashMap<u32, String> = HashMap::new();
    for name in bam.header().target_names() {
        let r_tid = bam.header().tid(name).unwrap();
        let r_string = convert_u82String(name);
        headers.insert(r_tid, r_string);
    }
    
    for rd in bam.records() {
        let r = rd.unwrap();
        if r.is_unmapped() {
            continue;
        }
        let read_id = convert_u82String(r.qname());
        let reference_id = r.tid() as u32;
        let reference_name = headers.get(&reference_id).unwrap(); 
        let ref_start = r.pos();
        let ref_end = r.cigar().end_pos();
        let read_strand = if r.is_reverse() {
            "-"
        } else {
            "+"
        };


        let info = if let Some(value) = filtered_alignments.get(&read_id) {
            value
        } else {
            continue;
        };
        
        for i in info.iter() {
            if *reference_name != i.0 {
                continue;
            }
            if ref_start != i.1 {
                continue;
            }
            if ref_end != i.2 {
                continue;
            }
            if read_strand != i.5 {
                continue
            }

            let mut record = bam::record::Record::new();
            // reference_id
            record.set_tid(r.tid());
            // reference_start
            record.set_pos(r.pos());
            // flag
            if r.is_reverse() {
                if r.is_supplementary() {
                    let f: u16 = 2064;
                    record.set_flags(f);
                } else {
                    let f: u16 = 16;
                    record.set_flags(f);
                }
            } else {
                if r.is_supplementary() {
                    let f: u16 = 2048;
                    record.set_flags(f);
                } else {
                    let f: u16 = 0;
                    record.set_flags(f);
                }
            }
            // qname, cigar, query_sequence, quality
            let cigar_string: bam::record::CigarString = bam::record::CigarString::from(r.cigar().iter().cloned().collect::<Vec<_>>());
            record.set(r.qname(), Some(&cigar_string), &r.seq().as_bytes(), r.qual());
            // mapping quality
            record.set_mapq(r.mapq());
            // next_reference_id
            // next_reference_start
            // template length
            // tags
            for aux in r.aux_iter() {
                let (tag, value) = aux.unwrap();
                record.push_aux(tag, value).unwrap();
            }
            out.write(&record).unwrap();
            break;
        }
    }
}
