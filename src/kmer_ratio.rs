use rust_htslib::{bam, bam::Read, htslib};
use rust_htslib::bam::record::Aux;
use std::error::Error;
use std::collections::HashMap;


pub fn run(
    input_bam: &str,
    threads: usize,
    prior_mean: f64,
    prior_weight: f64,
) -> Result<(), Box<dyn Error>> {
    if !(0.0..=1.0).contains(&prior_mean) {
        return Err(format!("--prior-mean must be within [0.0, 1.0], got {}", prior_mean).into());
    }
    if prior_weight < 0.0 {
        return Err(format!("--prior-weight must not be negative, got {}", prior_weight).into());
    }
    let mut kmers: HashMap<String, Vec<u32>> = HashMap::new();
    let mut bam = bam::Reader::from_path(input_bam).expect(&format!("Could not open {}", input_bam));

    bam.set_threads(threads).expect(&format!("Failure set {} threads", threads));

    let filter_closure: Box<dyn Fn(&bam::Record) -> bool> = Box::new(|record: &bam::Record| {
        record.flags() & (htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY) as u16 == 0
    });

    let mut all_counts = 0;
    for read in bam
        .rc_records()
        .map(|r| r.expect("Failure parsing Bam file"))
        .inspect(|_| all_counts += 1)
        .filter(|read| filter_closure(read))
    {
        let read_id: String = String::from_utf8_lossy(read.qname()).to_string();
        //let seq: &[u8] = read.seq().as_bytes();

        let rk_tag = get_rk_tag(&read);
        let pk_sk_tag = get_pk_sk_tag(&read);

        if let Some(value) = kmers.get_mut(&read_id) {
            value[0] += rk_tag;
            value[1] += pk_sk_tag;
        } else {
            kmers.insert(read_id, vec![rk_tag, pk_sk_tag]);
        }
    
    }

    // Fraction of the loci a read could observe that it matched, smoothed by the
    // prior so that available == 0 falls back to prior_mean instead of 1.0 and
    // 1/1 claims less than 1000/1000.
    for (key, value) in kmers.iter() {
        let freq: f64 =
            (value[1] as f64 + prior_mean * prior_weight) / (value[0] as f64 + prior_weight);
        println!("{}\t{}\t{}\t{}", key, value[0], value[1], freq);
    }
    Ok(())
}


fn get_rk_tag(record: &bam::Record) -> u32 {
    match record.aux(b"RK") {
        Ok(value) => match value {
            Aux::U32(v) => v,
            _ => panic!("Unexpected type of  RK tag: {:?}", value),
        },
        Err(_e) => panic!("Colud not access the RK tag"),
    }
}

fn get_pk_sk_tag(record: &bam::Record) -> u32 {
    if record.is_supplementary() {
        match record.aux(b"SK") {
            Ok(value) => match value {
                Aux::String(v) => {
                    let kmers: Vec<u32> = v.split(',')
                                           .map(|s| s.parse::<u32>().unwrap())
                                           .collect();
                    
                    if let Some(max_value) = kmers.iter().cloned().max() {
                        max_value
                    } else {
                        panic!("No values for the SK tags");
                    }
                },
                _ => panic!("Unexpected type of  SK tag: {:?}", value),
            },
            Err(_e) => panic!("Colud not access the SK tag"),
        }
    } else {
        match record.aux(b"PK") {
            Ok(value) => match value {
                Aux::String(v) => {
                    let kmers: Vec<u32> = v.split(',')
                                           .map(|s| s.parse::<u32>().unwrap())
                                           .collect();
                    
                    if let Some(max_value) = kmers.iter().cloned().max() {
                        max_value
                    } else {
                        panic!("No values for the PK tags");
                    }
                },
                _ => panic!("Unexpected type of  PK tag: {:?}", value),
            },
            Err(_e) => panic!("Colud not access the PK tag"),
        }
    }
}