use clap::Parser;
use std::hash::Hash;
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
