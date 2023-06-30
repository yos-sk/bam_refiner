use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::HashMap;

use marker_filter::convert_u82String;
use marker_filter::get_cigartuples;
use marker_filter::get_read_position;
use marker_filter::reverse_complement;

pub fn process_write_bam(
    bamfile: &str,
    output_bam: &str,
    filtered_alignments: &HashMap<
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
    >,
) {
    let mut bam = bam::Reader::from_path(bamfile).expect(&format!("Could not open {}", bamfile));
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(output_bam, &header, bam::Format::Bam)
        .expect(&format!("Could not open {}", output_bam));

    let mut headers: HashMap<u32, String> = HashMap::new();
    for name in bam.header().target_names() {
        let r_tid = bam.header().tid(name).unwrap();
        let r_string = convert_u82String(name);
        headers.insert(r_tid, r_string);
    }
    let mut prev_read_id = String::new();
    let mut read_alignments: Vec<bam::record::Record> = Vec::new();
    let mut line_num = 0;
    for rd in bam.records() {
        let r = rd.unwrap();
        line_num += 1;
        let read_id = convert_u82String(r.qname());
        eprintln!("Processing line {}, {}", line_num, &read_id);
        if r.is_unmapped() {
            continue;
        }

        if prev_read_id.len() == 0 {
            prev_read_id = read_id;
            read_alignments.push(r);
            continue;
        }

        if read_id != prev_read_id {
            write_bam(
                &prev_read_id,
                &mut out,
                &read_alignments,
                &headers,
                filtered_alignments,
            );
            read_alignments = Vec::new();
            read_alignments.push(r);
            prev_read_id = read_id;
        } else {
            read_alignments.push(r);
        }
    }
    write_bam(
        &prev_read_id,
        &mut out,
        &read_alignments,
        &headers,
        filtered_alignments,
    );
}

fn write_bam(
    read_id: &str,
    out: &mut bam::Writer,
    read_alignments: &Vec<bam::record::Record>,
    headers: &HashMap<u32, String>,
    filtered_alignments: &HashMap<
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
    >,
) {
    let mut records = read_alignments.clone();
    records.sort_by_key(|t| t.flags());

    let mut sequence = String::new();
    let mut qual: Option<&[u8]> = None;
    for (i, r) in records.iter().enumerate() {
        if i == 0 {
            if r.is_reverse() {
                sequence = reverse_complement(&convert_u82String(&r.seq().as_bytes()));
            } else {
                sequence = convert_u82String(&r.seq().as_bytes());
            }
            qual = Some(r.qual());
        }

        let mut read_seq = String::new();
        if r.is_secondary() {
            if r.is_reverse() {
                if sequence.len() != 0 {
                    read_seq = reverse_complement(&sequence);
                } else {
                    eprintln!("Bad alignment: BamFile needs to be sorted");
                }
            } else {
                if sequence.len() != 0 {
                    read_seq = sequence.clone();
                } else {
                    eprintln!("Bad alignment: BamFile needs to be sorted");
                }
            }
        }

        let reference_id = r.tid() as u32;
        let reference_name = headers.get(&reference_id).unwrap();
        let ref_start = r.pos();
        let ref_end = r.cigar().end_pos();
        let read_strand = if r.is_reverse() { "-" } else { "+" };

        let cigartuples = get_cigartuples(&r);
        let read_pos: (u32, u32) = get_read_position(&cigartuples);
        let read_start = read_pos.0;
        let read_end = read_pos.1;

        let info = if let Some(value) = filtered_alignments.get(read_id) {
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
                continue;
            }
            if read_start != i.3 {
                continue;
            }
            if read_end != i.4 {
                continue;
            }

            let mut record = bam::record::Record::new();
            // reference_id
            record.set_tid(r.tid());
            // reference_start
            record.set_pos(r.pos());
            // flag
            if r.is_reverse() {
                if i.7 == 1 {
                    let f: u16 = 2064;
                    record.set_flags(f);
                } else {
                    let f: u16 = 16;
                    record.set_flags(f);
                }
            } else {
                if i.7 == 1 {
                    let f: u16 = 2048;
                    record.set_flags(f);
                } else {
                    let f: u16 = 0;
                    record.set_flags(f);
                }
            }
            // qname, cigar, query_sequence, quality
            let cigar_string: bam::record::CigarString =
                bam::record::CigarString::from(r.cigar().iter().cloned().collect::<Vec<_>>());
            if r.is_secondary() {
                let bytes: Vec<u8> = read_seq.into_bytes();
                if let Some(quality) = qual {
                    record.set(r.qname(), Some(&cigar_string), &bytes, &quality);
                }
            } else {
                record.set(
                    r.qname(),
                    Some(&cigar_string),
                    &r.seq().as_bytes(),
                    r.qual(),
                );
            }
            // mapping quality
            if i.9 == 0 {
                let mapq: u8 = 30;
                record.set_mapq(mapq);
            } else {
                let mapq: u8 = 60;
                record.set_mapq(mapq);
            }

            // next_reference_id
            // next_reference_start
            // template length
            // tags
            for aux in r.aux_iter() {
                let (tag, value) = aux.unwrap();
                record.push_aux(tag, value).unwrap();
            }

            // New tag: HP
            if i.9 == 0 {
                let aux_hp_tag = bam::record::Aux::String("Amb");
                record.push_aux(b"HP", aux_hp_tag).unwrap();
            } else {
                if &reference_name[0..2] == "h1" {
                    let aux_hp_tag = bam::record::Aux::String("HP1");
                    record.push_aux(b"HP", aux_hp_tag).unwrap();
                } else {
                    let aux_hp_tag = bam::record::Aux::String("HP2");
                    record.push_aux(b"HP", aux_hp_tag).unwrap();
                }
            }
            out.write(&record).unwrap();
            break;
        }
    }
}
