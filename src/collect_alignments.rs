use rust_htslib::{bam, bam::Read, htslib};
use std::error::Error;
use std::collections::HashMap;
use bam_refiner::get_cigartuples;
use bam_refiner::get_read_position;
use bam_refiner::reverse_complement;
use bam_refiner::Data;

pub fn run(input_bam: &str, threads: usize) -> Result<(Vec<Data>, HashMap<String, Vec<u8>>), Box<dyn Error>> {
    let mut records: Vec<Data> = Vec::new(); 
    let mut sequences: HashMap<String, Vec<u8>> = HashMap::new();
    let mut bam = bam::Reader::from_path(input_bam).expect(&format!("Could not open {}", input_bam));

    let header = bam.header().clone();
    bam.set_threads(threads).expect(&format!("Failure set {} threads", threads));

    let filter_closure: Box<dyn Fn(&bam::Record) -> bool> = Box::new(|record: &bam::Record| {
        record.flags() & htslib::BAM_FUNMAP as u16 == 0
    });

    let mut all_counts = 0;
    for record in bam
        .rc_records()
        .map(|r| r.expect("Failure parsing Bam file"))
        .inspect(|_| all_counts += 1)
        .filter(|read| filter_closure(read))
    {
        // reference name
        let ref_name = String::from_utf8_lossy(header.tid2name(record.tid() as u32)).to_string();

        // reference start
        let ref_start = record.pos();

        // reference end
        let ref_end = record.cigar().end_pos();

        // cigar_tuples
        let cigartuples = get_cigartuples(&record);

         // read name
        let read_id: String = String::from_utf8_lossy(record.qname()).to_string();

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

        // sequence of read
        if !record.is_supplementary() && !record.is_secondary() {
            let seq: Vec<u8> = if record.is_reverse() {
                reverse_complement(&record.seq().as_bytes())
            } else {
                record.seq().as_bytes()
            };
            sequences.insert(read_id.clone(), seq);
        }

        eprintln!("{},{},{},{},{},{},{:?},{:?},{:?}", read_id, ref_name, ref_start, ref_end, r_read_start, r_read_end, record.is_reverse(), record.is_supplementary(), record.is_secondary());

        // Save the above infromation as Data structure
        let save = Data {
            reference_name: ref_name,
            reference_start: ref_start,
            reference_end: ref_end,
            read_name: read_id,
            read_start: r_read_start,
            read_end: r_read_end,
            is_reverse: record.is_reverse(),
            is_supplementary: record.is_supplementary(),
            is_secondary: record.is_secondary(),
            cigar_tuples: cigartuples, 
            rk_cnt: 0,
            pk_sk_cnt: 0,
        };
        
        records.push(save);
        
    }
    Ok((records, sequences))
}