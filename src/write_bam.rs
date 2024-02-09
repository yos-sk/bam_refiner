use rust_htslib::{bam, bam::Read, htslib};
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;

use bam_refiner::Data;
use bam_refiner::get_read_name_list;
use bam_refiner::reverse_complement;
use bam_refiner::get_read_position;
use bam_refiner::get_cigartuples;


pub fn run(
    bamfile: &str,
    output_bam: &str,
    refined_alignments: &HashMap<String, Vec<Data>>,
    sequences: &HashMap<String, Vec<u8>>,
    qualities: &HashMap<String, Vec<u8>>,
    hap1_list: &str,
    hap2_list: &str,
) -> Result<(), Box<dyn Error>> {

    let hap1_set: HashSet<String> = get_read_name_list(hap1_list).expect(&format!("Could not read {}", hap1_list));
    let hap2_set: HashSet<String> = get_read_name_list(hap2_list).expect(&format!("Could not read {}", hap2_list));

    let mut bam = bam::Reader::from_path(bamfile).expect(&format!("Could not open {}", bamfile));
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(output_bam, &header, bam::Format::Bam)
        .expect(&format!("Could not open {}", output_bam));

    let header = bam.header().clone();

    let filter_closure: Box<dyn Fn(&bam::Record) -> bool> = Box::new(|record: &bam::Record| {
        record.flags() & htslib::BAM_FUNMAP as u16 == 0
    });

    for record in bam.records()
        .map(|r| r.expect("Failure parsing Bam file"))
        .filter(|read| filter_closure(read)) {
        
        // read name from bam
        let read_id: String = String::from_utf8_lossy(record.qname()).to_string();

        // reference name
        let ref_name = String::from_utf8_lossy(header.tid2name(record.tid() as u32)).to_string();

        // reference start
        let ref_start = record.pos();

        // reference end
        let ref_end = record.cigar().end_pos();

        // cigar_tuples
        let cigartuples = get_cigartuples(&record);

        // read position
        let read_pos: (u32, u32, u32) = get_read_position(&cigartuples);
        let r_start = read_pos.0;
        let r_end = read_pos.1;
        let r_read_length = read_pos.2;
        let r_read_start = if !record.is_reverse() {
            r_start
        } else {
            r_read_length - r_end
        };
        let r_read_end = if !record.is_reverse() {
            r_end
        } else {
            r_read_length - r_start
        };
        // sequence
        let sequence = if let Some(value) = sequences.get(&read_id) {
            value
        } else {
            panic!("Read name: {} doesn't have sequence", read_id);
        };

        // quality
        let quality = if let Some(value) = qualities.get(&read_id) {
            value
        } else {
            panic!("Read name: {} doesn't have sequence", read_id);
        };

        if let Some(value) = refined_alignments.get(&read_id) {
            for alignment in value.iter() {
                if ref_name != alignment.reference_name {
                    continue;
                }
                if ref_start != alignment.reference_start {
                    continue;
                }
                if ref_end != alignment.reference_end {
                    continue;
                }
                if record.is_reverse() != alignment.is_reverse {
                    continue;
                }
                if r_read_start != alignment.read_start {
                    continue;
                }
                if r_read_end != alignment.read_end {
                    continue;
                }
                let mut out_record = bam::record::Record::new();
                // reference_id: i32
                out_record.set_tid(record.tid());
                // reference_start: i64
                out_record.set_pos(record.pos());
                // flag: u16
                if record.is_reverse() {
                    if alignment.is_supplementary {
                        out_record.set_flags(2064);
                    } else {
                        out_record.set_flags(16);
                    }
                } else {
                    if alignment.is_supplementary {
                        out_record.set_flags(2048);
                    } else {
                        out_record.set_flags(0);
                    }
                }
                // qname: &[u8], cigar: Option<&CigarString>, query_sequence: &[u8], quality: &[u8]
                let new_cigar = make_new_cigar(record.cigar(), record.is_secondary(), alignment.is_supplementary);
                let new_seq = make_new_seq(&sequence, alignment.is_reverse, record.is_secondary(), alignment.is_supplementary, &new_cigar);
                let new_qual = make_new_qual(&quality, alignment.is_reverse, record.is_secondary(), alignment.is_supplementary, &new_cigar);

                
                out_record.set(
                    record.qname(),
                    Some(&CigarString(new_cigar)),
                    &new_seq,
                    &new_qual,
                );
    
                // mapping quality
                let max_kmer_value = alignment.pk_sk_vec.iter().max().unwrap();
                let count_of_max = alignment.pk_sk_vec.iter().filter(|&x| *x == *max_kmer_value).count();
                let mapq: u8 = 60 / count_of_max as u8;
                out_record.set_mapq(mapq);

                // next_reference_id
                // next_reference_start
                // template length
                // tags
                for aux in record.aux_iter() {
                    let (tag, value) = aux.unwrap();
                    out_record.push_aux(tag, value).unwrap();
                }

                // New tag: HP
                if *max_kmer_value == 0  {
                    let aux_hp_tag = bam::record::Aux::U8(0);
                    out_record.push_aux(b"HP", aux_hp_tag).unwrap();
                } else {
                    if hap1_set.contains(&ref_name) {
                        let aux_hp_tag = bam::record::Aux::U8(1);
                        out_record.push_aux(b"HP", aux_hp_tag).unwrap();
                    } else if hap2_set.contains(&ref_name) {
                        let aux_hp_tag = bam::record::Aux::U8(2);
                        out_record.push_aux(b"HP", aux_hp_tag).unwrap();
                    } else {
                        eprintln!("{} does not contain in hap1 and hap2 list", ref_name);
                    }
                }

                // New tag: PK and SK
                if alignment.is_supplementary {
                    let supp_kmers: String = alignment.pk_sk_vec.iter().map(|&x| x.to_string()).collect::<Vec<String>>().join(",");
                    let supp_kmers_tag = bam::record::Aux::String(&supp_kmers);
                    out_record.push_aux(b"SK", supp_kmers_tag).unwrap();
                } else {
                    let prim_kmers: String = alignment.pk_sk_vec.iter().map(|&x| x.to_string()).collect::<Vec<String>>().join(",");
                    let prim_kmers_tag = bam::record::Aux::String(&prim_kmers);
                    out_record.push_aux(b"PK", prim_kmers_tag).unwrap();
                }

                // New tag: reference kmer
                let aux_rk_tag = bam::record::Aux::U32(alignment.rk_cnt as u32);
                out_record.push_aux(b"RK", aux_rk_tag).unwrap();
                out.write(&out_record).unwrap();
                break;
            }
        } else {
            continue;
        }
    }
    Ok(())
}

fn make_new_cigar(cigar: CigarStringView, input_is_secondary: bool, output_is_supplementary: bool) -> Vec<Cigar> {
    let mut new_cigar: Vec<Cigar> = Vec::new();

    for op in cigar.iter() {
        match op {
            Cigar::SoftClip(len) => {
                if input_is_secondary && output_is_supplementary {
                    new_cigar.push(Cigar::HardClip(*len));
                } else {
                    new_cigar.push(*op);
                }
            },
            _ => new_cigar.push(*op),
        }
    }
    new_cigar
}

fn make_new_seq(sequence: &Vec<u8>, is_reverse: bool, input_is_secondary: bool, output_is_supplementary: bool, cigarstring: &Vec<Cigar>) -> Vec<u8> {
    let mut start: usize = 0;
    let mut end: usize = sequence.len();
    let seq = if is_reverse {
        reverse_complement(sequence)
    } else {
        sequence.to_vec()
    };
    
    for (i, op) in cigarstring.iter().enumerate() {
        match op {
            Cigar::SoftClip(len) => {
                if input_is_secondary && output_is_supplementary {
                    if i == 0 {
                        start = *len as usize;
                    } else if i == cigarstring.len() - 1 {
                        end = end - *len as usize;
                    }
                }
            },
            Cigar::HardClip(len) => {
                if i == 0 {
                    start = *len as usize;
                } else if i == cigarstring.len() - 1 {
                    end = end - *len as usize;
                }
            }
            _ => (),
        }
    }

    seq[start..end].to_vec()
}

fn make_new_qual(quality: &Vec<u8>, is_reverse: bool, input_is_secondary: bool, output_is_supplementary: bool, cigarstring: &Vec<Cigar>) -> Vec<u8> {
    let mut start: usize = 0;
    let mut end: usize = quality.len();

    let mut qual = quality.to_vec();
    if is_reverse {
        qual.reverse();
    }

    for (i, op) in cigarstring.iter().enumerate() {
        match op {
            Cigar::SoftClip(len) => {
                if input_is_secondary && output_is_supplementary {
                    if i == 0 {
                        start = *len as usize;
                    } else if i == cigarstring.len() - 1 {
                        end = end - *len as usize;
                    }
                }
            },
            Cigar::HardClip(len) => {
                if i == 0 {
                    start = *len as usize;
                } else if i == cigarstring.len() - 1 {
                    end = end - *len as usize;
                }
            }
            _ => (),
        }
    }
    qual[start..end].to_vec()
}