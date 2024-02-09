use std::collections::HashMap;
use std::error::Error;

use bam_refiner::Data;

pub fn run(
    alignments: Vec<Data>,
    indices: HashMap<String, Vec<usize>>,
) -> Result<HashMap<String, Vec<Data>>, Box<dyn Error>> {
    let mut new_records: HashMap<String, Vec<Data>> = HashMap::new();

    for (read_id, index) in indices.iter() {
        let mut primary: Vec<Data> = Vec::new();
        let mut secondary: Vec<Data> = Vec::new();
        let mut supplementary: Vec<Data> = Vec::new();

        let mut prim_kmer_cnts: Vec<usize> = Vec::new();
        let mut supp_kmer_cnts: Vec<Vec<usize>> = Vec::new();

        for i in index.iter() {
            assert_eq!(
                alignments[*i].read_name,
                read_id.to_string(),
                "Read names are inconsistent: {} {}",
                alignments[*i].read_name,
                read_id
            );
            // eprintln!("Refine: {} {} {}", alignments[*i].read_name, alignments[*i].pk_sk_cnt, alignments[*i].rk_cnt);
            if alignments[*i].is_secondary {
                secondary.push(alignments[*i].clone());
            } else if alignments[*i].is_supplementary {
                supplementary.push(alignments[*i].clone());
                supp_kmer_cnts.push(vec![alignments[*i].pk_sk_cnt]);
            } else {
                primary.push(alignments[*i].clone());
                prim_kmer_cnts.push(alignments[*i].pk_sk_cnt);
            }
        }

        assert!(
            !primary.is_empty(),
            "Read: {} doesn't have any primary alignment",
            read_id
        );

        // Compare secondary alignments with primary/supplementary alignments
        for t_secondary in secondary.iter_mut() {
            let mut supp_id: isize = -1;
            let mut supp_dist: isize = -1;
            let mut supp_cnt: usize = 0;
            let prim_cnt: usize = primary[0].pk_sk_cnt;

            // secondary vs primary
            let diff_start = primary[0].read_start as isize - t_secondary.read_start as isize;
            let diff_end = primary[0].read_end as isize - t_secondary.read_end as isize;
            let prim_dist = diff_start.abs() + diff_end.abs();

            // secondary vs supplementary
            if !supplementary.is_empty() {
                for (i, t_supplementary) in supplementary.iter().enumerate() {
                    let s_diff_start =
                        t_supplementary.read_start as isize - t_secondary.read_start as isize;
                    let s_diff_end =
                        t_supplementary.read_end as isize - t_secondary.read_end as isize;
                    let t_supp_dist = s_diff_start.abs() + s_diff_end.abs();
                    if supp_dist == -1 || supp_dist < t_supp_dist {
                        supp_dist = t_supp_dist;
                        supp_id = i as isize;
                        supp_cnt = t_supplementary.pk_sk_cnt;
                    }
                }
            }

            // define what the secondary alignment is for
            if supp_dist == -1 || prim_dist <= supp_dist {
                // secondary alignment for primary alignment
                if prim_cnt < t_secondary.pk_sk_cnt {
                    t_secondary.is_secondary = false;
                    primary[0] = t_secondary.clone();
                }
                prim_kmer_cnts.push(t_secondary.pk_sk_cnt);
            } else {
                // secondary alignment for supplementary alignment
                if supp_cnt < t_secondary.pk_sk_cnt {
                    t_secondary.is_secondary = false;
                    t_secondary.is_supplementary = true;
                    supplementary[supp_id as usize] = t_secondary.clone();
                }
                supp_kmer_cnts[supp_id as usize].push(t_secondary.pk_sk_cnt);
            }
        }
        primary[0].pk_sk_vec = prim_kmer_cnts;
        if !supplementary.is_empty() {
            for (i, t_supplementary) in supplementary.iter_mut().enumerate() {
                t_supplementary.pk_sk_vec = supp_kmer_cnts[i].clone();
            }
        }
        new_records.insert(read_id.to_string(), vec![primary, supplementary].concat());
    }

    // Print new alignment results to standard output
    for (read_id, records) in new_records.iter() {
        print!("{}", read_id);
        for record in records.iter() {
            print!(
                "\t{},{},{},{},{},{:?},{:?},{:?},{},{},{:?}",
                record.reference_name,
                record.reference_start,
                record.reference_end,
                record.read_start,
                record.read_end,
                record.is_reverse,
                record.is_secondary,
                record.is_supplementary,
                record.pk_sk_cnt,
                record.rk_cnt,
                record.pk_sk_vec,
            )
        }
        print!("\n");
    }
    Ok(new_records)
}
