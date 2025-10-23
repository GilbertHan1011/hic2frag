use clap::Parser;
use std::path::PathBuf;
use noodles_bam as bam;
use noodles_sam::alignment::Record;
use std::error::Error;
use bed_utils::bed::{io::Reader, BEDLike, BED};
use bed_utils::intervaltree::{Interval, Lapper};
use std::collections::HashMap;
use std::fs::File;
use log::warn;

#[derive(Parser, Debug)]
#[clap(author = "GilbertHan", version, about = "Bam to HiC fragments")]
struct Cli {

    #[clap(short, long)]
    fragment_file: Vec<String>,

    #[clap(short, long)]
    bam: Vec<String>,

    #[clap(short, long)]
    out_dir: Option<PathBuf>,
}

fn get_read_pos(read:&bam::Record, st: &str) -> Option<usize>{
    match st {
        "middle" => {
            let start =  read.alignment_start().transpose().ok().flatten().unwrap().get();
            let start0based = start -1;
            let span = read.alignment_span().transpose().ok().flatten().unwrap();
            let pos = start0based + span / 2;
            Some(pos)
        }
        "start" => {Some(get_read_start(read))}
        "left" => {Some(read.alignment_start().transpose().ok().flatten().unwrap().get())}
        _ => None,
    }
}

/// return 5' start of read
fn get_read_start(read: &bam::Record) -> usize{
    let start = read.alignment_start().transpose().ok().flatten().unwrap().get();
    let pos = if read.flags().is_reverse_complemented(){
        let span = read.alignment_span().transpose().ok().flatten().unwrap();
        start + span - 1
    } else{
        start
    };
    pos
}

fn get_read_strand(read: &bam::Record) -> &str{
    let mut strand = "+";
    if read.flags().is_reverse_complemented(){
        strand = "-"
    }
    strand
}

fn is_intra_chrom(read1:  &bam::Record, read2 :  &bam::Record) -> Option<bool>{
    let tid1_opt = read1.reference_sequence_id().transpose().ok().flatten().unwrap();
    let tid2_opt = read2.reference_sequence_id().transpose().ok().flatten().unwrap();
    
    if tid2_opt == tid2_opt {
        return Some(false)
    } else{
        return Some(true)
    }
}

/*
Calculte the contact distance between two intrachromosomal reads

    read1 : [Record]
    read2 : [Record]

*/    
fn get_cis_distance(read1:  &bam::Record, read2 :  &bam::Record) -> Option<usize>{
    let mut dist = None;
    let unmap1 = read1.flags().is_unmapped();
    let unmap2 = read2.flags().is_unmapped();
    if !unmap1 && !unmap2{
        let r1pos = get_read_pos(read1,"middle").unwrap();
        let r2pos = get_read_pos(read2, "middle").unwrap();
        if r1pos > r2pos {
            dist = Some(r1pos - r2pos);
        } else {
            dist = Some(r2pos - r1pos);
        }
    }
    dist
}

fn get_ordered_reads<'a>(read1: &'a bam::Record, read2: &'a bam::Record) 
    -> Option<(&'a bam::Record, &'a bam::Record)>{
    let tid1_opt = read1.reference_sequence_id().transpose().ok().flatten();
    let tid2_opt = read2.reference_sequence_id().transpose().ok().flatten();
    tid1_opt.zip(tid2_opt).and_then(|(tid1, tid2)| {
        if tid1 < tid2 {
            Some((read1,read2))
        } else if tid1 > tid2{
            Some((read2, read1))
        } else {
            let r1pos_opt = get_read_pos(read1, "middle");
            let r2pos_opt = get_read_pos(read2, "middle");
            r1pos_opt.zip(r2pos_opt).map(|(r1pos, r2pos)| {
                // We have valid positions. Sort by position.
                if r1pos <= r2pos {
                    (read1, read2)
                } else {
                    (read2, read1)
                }
            })
        }
    })
}

fn are_contiguous_fragments(frag1:BED<3>, frag2:BED<3>, chr1:usize, chr2:usize) -> bool{
    if chr1 != chr2{
        return false
    } 
    let frag1_touch_frag2 = frag1.end() == frag2.start();
    let frag2_touch_frag1 = frag1.start() == frag2.end();
    frag2_touch_frag1 || frag1_touch_frag2
}

fn is_religation(read1: &bam::Record, read2: &bam::Record,frag1:BED<3>, frag2:BED<3>)-> bool{
    let tid1_opt = read1.reference_sequence_id().transpose().ok().flatten().unwrap();
    let tid2_opt = read2.reference_sequence_id().transpose().ok().flatten().unwrap();
    let mut ret = false;
    if are_contiguous_fragments(frag1, frag2, tid1_opt,tid2_opt){
        ret = true
    }
    ret
}

fn is_self_circle(read1: &bam::Record, read2: &bam::Record) -> bool{
    let mut ret = false;
    if let Some((r1,r2)) = get_ordered_reads(read1, read2) {
        (get_read_strand(r1) == "-") && (get_read_strand(r2) == "+")
    } else{
        false
    }
}


fn get_overlapping_restriction_fragment(res_frag : &HashMap<String, Lapper<u64,BED<3>>>, 
    chrom: &str, read:  &bam::Record) -> Option<BED<3>> {
    let pos = get_read_pos(read, "middle").unwrap();
    if let Some(lapper) = res_frag.get(chrom) {
        let overlapping_frag : Vec<_> = lapper.find(pos as u64, pos as u64 + 1).collect();
        if overlapping_frag.len() > 1{
            warn!("Warning: {} restriction fragments found for {} - skipped", 
                     overlapping_frag.len(), read.name().unwrap().to_string());
            return(None)
        } else if overlapping_frag.len() == 0 {
            warn!("Warning: {} restriction fragments found for {} - skipped", 
                     overlapping_frag.len(), read.name().unwrap().to_string());
            return(None)
        } else{
            //let test = &overlapping_frag[0].val;
            return(Some(overlapping_frag[0].val.clone()))
        }
    } else{
        warn!("Warning: No restriction fragments found for {} - skipped", 
                     read.name().unwrap().to_string());
        return(None)
    }
}

fn convert_vec_to_lapper<B: BEDLike + Clone>(
    bed_records: &Vec<B>
 ) -> HashMap<String, Lapper<u64,B>> {
    let mut chrom_to_interval : HashMap<String, Vec<Interval<u64, B>>> = HashMap::new();

    for bed in bed_records {
        let chrom = bed.chrom().to_string();
        let interval = Interval {
            start : bed.start(),
            stop : bed.end(),
            val : bed.clone()
        };

        chrom_to_interval
            .entry(chrom)
            .or_insert_with(Vec::new)
            .push(interval);
    }

    chrom_to_interval
        .into_iter()
        .map(|(chrom,intervals)| (chrom, Lapper::new(intervals)))
        .collect()
}


fn main() ->  Result<(), Box<dyn Error>>{
    let mut reader = bam::io::reader::Builder::default().build_from_path("/home/zhanglab/example/Rosalind_rust/hic2frag/data/test.bam")?;
    let bed_file = "/home/zhanglab/example/Rosalind_rust//tmp/test.bed";
    let bed_file_open = File::open(bed_file).unwrap();
    let mut bed_reader = Reader::new(bed_file_open,None);
    let bed_rec:Vec<BED<3>> = bed_reader.into_records::<BED<3>>()
        .map(|r| r.unwrap())
        .collect();
    let bed_ladder = convert_vec_to_lapper(&bed_rec);
    let test = &bed_rec[0];
    println!("Bed file created: {},{},{}", &test.chrom(),&test.start(),&test.end());

    let headers =  reader.read_header()?;
    for result in reader.records() {
        let record = result?;
        
        if (record.flags().is_first_segment()){
            let r1 = record;
            if (!r1.flags().is_unmapped()){
                if let Some(result) = r1.reference_sequence(&headers){
                    let (name_bytes,_) = result?;
                    let name_str = std::str::from_utf8(name_bytes)?;
                }
            }
        }
    }

        Ok(())
}
