use clap::Parser;
use std::path::PathBuf;
use noodles_bam as bam;
use noodles_sam::alignment::Record;
use std::error::Error;
use bed_utils::bed::{io::Reader, BEDLike, BED};
use bed_utils::intervaltree::{Interval, Lapper};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use log::{warn, info};
use std::ptr;

#[derive(Debug, Default)]
struct Statistics {
    reads_counter: u64,
    de_counter: u64,
    re_counter: u64,
    sc_counter: u64,
    valid_counter: u64,
    valid_counter_ff: u64,
    valid_counter_rr: u64,
    valid_counter_fr: u64,
    valid_counter_rf: u64,
    single_counter: u64,
    dump_counter: u64,
    filt_counter: u64,
    // Allele specific counters
    g1g1_ascounter: u64,
    g2g2_ascounter: u64,
    g1u_ascounter: u64,
    ug1_ascounter: u64,
    g2u_ascounter: u64,
    ug2_ascounter: u64,
    g1g2_ascounter: u64,
    g2g1_ascounter: u64,
    uu_ascounter: u64,
    cf_ascounter: u64,
}

struct OutputHandlers {
    valid: BufWriter<File>,
    de: Option<BufWriter<File>>,
    re: Option<BufWriter<File>>,
    sc: Option<BufWriter<File>>,
    dump: Option<BufWriter<File>>,
    single: Option<BufWriter<File>>,
    filt: Option<BufWriter<File>>,
    sam: Option<File>,
}

#[derive(Parser, Debug, Clone)]
#[clap(author = "GilbertHan", version, about = "Bam to HiC fragments")]
struct Cli {
    #[clap(short = 'f', long, help = "Restriction fragment file (BED format)")]
    fragment_file: String,

    #[clap(short = 'r', long, help = "BAM/SAM file of mapped reads")]
    bam: String,

    #[clap(short, long, help = "Output directory. Default is current directory")]
    out_dir: Option<PathBuf>,

    #[clap(short = 's', long, help = "Shortest insert size of mapped reads to consider")]
    min_insert_size: Option<u64>,

    #[clap(short = 'l', long, help = "Longest insert size of mapped reads to consider")]
    max_insert_size: Option<u64>,

    #[clap(short = 't', long, help = "Shortest restriction fragment length to consider")]
    min_frag_size: Option<u64>,

    #[clap(short = 'm', long, help = "Longest restriction fragment length to consider")]
    max_frag_size: Option<u64>,

    #[clap(short = 'd', long, help = "Minimum distance between intrachromosomal contact to consider")]
    min_cis_dist: Option<u64>,

    #[clap(short = 'g', long, help = "Genotype tag for allele specific classification")]
    gtag: Option<String>,

    #[clap(short = 'a', long, help = "Write all additional output files with discarded reads info")]
    all: bool,

    #[clap(short = 'S', long, help = "Output an additional SAM file with flag 'CT' for pairs classification")]
    sam: bool,

    #[clap(short, long, help = "Verbose output")]
    verbose: bool,
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
    let tid1_opt = read1.reference_sequence_id().transpose().ok().flatten();
    let tid2_opt = read2.reference_sequence_id().transpose().ok().flatten();
    
    if let (Some(tid1), Some(tid2)) = (tid1_opt, tid2_opt) {
        Some(tid1 == tid2)
    } else {
        None
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
        let r1pos = get_read_pos(read1,"start").unwrap();
        let r2pos = get_read_pos(read2, "start").unwrap();
        if r1pos > r2pos {
            dist = Some(r1pos - r2pos); //maintain same as hic-pro statistics
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
            let r1pos_opt = get_read_pos(read1, "start");
            let r2pos_opt = get_read_pos(read2, "start");
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

fn are_contiguous_fragments(frag1:BED<6>, frag2:BED<6>, chr1:usize, chr2:usize) -> bool{
    if chr1 != chr2{
        return false
    } 
    let frag1_touch_frag2 = frag1.end() == frag2.start();
    let frag2_touch_frag1 = frag1.start() == frag2.end();
    frag2_touch_frag1 || frag1_touch_frag2
}

fn is_religation(read1: &bam::Record, read2: &bam::Record,frag1:BED<6>, frag2:BED<6>)-> bool{
    let tid1_opt = read1.reference_sequence_id().transpose().ok().flatten().unwrap();
    let tid2_opt = read2.reference_sequence_id().transpose().ok().flatten().unwrap();
    let mut ret = false;
    if are_contiguous_fragments(frag1, frag2, tid1_opt,tid2_opt){
        ret = true
    }
    ret
}

fn is_self_circle(read1: &bam::Record, read2: &bam::Record) -> bool{
    if let Some((r1,r2)) = get_ordered_reads(read1, read2) {
        (get_read_strand(r1) == "-") && (get_read_strand(r2) == "+")
    } else{
        false
    }
}

/*
    Both reads are expected to be on the same restriction fragments
    Check the orientation of reads -><-

    read1 : [AlignedRead]
    read2 : [AlignedRead]
 */
fn is_dangling_end(read1: &bam::Record, read2: &bam::Record) -> bool{
    if let Some((r1, r2)) = get_ordered_reads(read1, read2) {
        get_read_strand(r1) == "+" && get_read_strand(r2) == "-"
    } else {
        false
    }
}

/*
    Both reads are expected to be on the different restriction fragments
    Check the orientation of reads ->-> / <-<- / -><- / <-->

    read1 : [AlignedRead]
    read2 : [AlignedRead]
 */

fn get_valid_orientation(read1: &bam::Record, read2: &bam::Record) -> &'static str{
    let (r1, r2) = get_ordered_reads(read1, read2).unwrap();
    let r1_strand = get_read_strand(r1);
    let r2_strand = get_read_strand(r2);
    let direction: &'static str = match (r1_strand, r2_strand){
        ("+", "+") => "FF",
        ("-", "-") => "RR",
        ("+", "-") => "FR",
        ("-", "+") => "RF",
        _ => "Unknown",
    };
    direction
}

fn get_pe_fragment_size(read1: &bam::Record, read2: &bam::Record, 
    res_frag1: Option<BED<6>>, res_frag2: Option<BED<6>>,
    interaction_type: &str) -> Option<u64>{
 // 1. Get ordered reads. If this is None, the chain stops and returns None.
    get_ordered_reads(read1, read2).and_then(|(r1, r2)| {
        
        // 2. Pair up fragments with the *ordered* reads.
        //    (This now assigns references, which is cheap and safe)
        let (rfrag1, rfrag2) = if ptr::eq(r1, read2) { // Check if r1 is the original read2
            (res_frag2, res_frag1) 
        } else {
            (res_frag1, res_frag2)
        };
        
        // Check if we have valid fragments
        let (rfrag1, rfrag2) = match (rfrag1, rfrag2) {
            (Some(f1), Some(f2)) => (f1, f2),
            _ => return None,
        };
        let r1pos_opt = get_read_pos(r1, "start");
        let r2pos_opt = get_read_pos(r2, "start");

        r1pos_opt.zip(r2pos_opt).map(|(r1pos_usize, r2pos_usize)| {
            
            // These variables are defined *inside* this closure
            let r1pos = r1pos_usize as u64;
            let r2pos = r2pos_usize as u64;

            // 4. Calculate size. Use `saturating_sub` to PREVENT panics.
            if interaction_type == "DE" || interaction_type == "RE" {
                // r2 is the rightmost read, so r2pos >= r1pos
                r2pos.saturating_sub(r1pos)
            } else if interaction_type == "SC" {
                let d1 = r1pos.saturating_sub(rfrag1.start() + 1);
                let d2 = rfrag2.end().saturating_sub(r2pos + 1);
                d1 + d2
            } else if interaction_type == "VI" {
                let dr1 = if get_read_strand(r1) == "+" {
                    rfrag1.end().saturating_sub(r1pos + 1)
                } else {
                    r1pos.saturating_sub(rfrag1.start()+1)
                };
                
                let dr2 = if get_read_strand(r2) == "+" {
                    rfrag2.end().saturating_sub(r2pos+1)
                } else {
                    r2pos.saturating_sub(rfrag2.start()+1)
                };
                dr1 + dr2
            } else {
                0 // Safe default
            }
        })
    })
}

fn get_interaction_type(read1: &bam::Record, _read1_chrom: &str, res_frag1 : Option<BED<6>>,
    read2: &bam::Record, _read2_chrom: &str, res_frag2 : Option<BED<6>>, _verbose: bool)-> Option<&'static str>{
        match (read1.flags().is_unmapped(), read2.flags().is_unmapped(),&res_frag1, &res_frag2) {
            (false, false, Some(rf1), Some(rf2)) =>{
                if ptr::eq(rf1, rf2) {
                    if is_self_circle(read1, read2){
                        Some("SC")
                    } else if is_dangling_end(read1, read2){
                        Some("DE")
                    } else{
                        None
                    }
                } else {
                    if is_religation(read1, read2, rf1.clone(), rf2.clone()) {
                        Some("RE")
                    } else {
                        // This is the valid interaction case
                        Some("VI")
                    }
                }
            },
            (true,_,_,_) | (_,true,_,_) =>{
                Some("SI")
            },
            _ => None,
        }
    }



fn get_overlapping_restriction_fragment(res_frag : &HashMap<String, Lapper<u64,BED<6>>>, 
    chrom: &str, read:  &bam::Record) -> Option<BED<6>> {
    let pos = get_read_pos(read, "middle").unwrap();
    if let Some(lapper) = res_frag.get(chrom) {
        let overlapping_frag : Vec<_> = lapper.find(pos as u64, pos as u64 + 1).collect();
        if overlapping_frag.len() > 1{
            warn!("Warning: {} restriction fragments found for {} - skipped", 
                     overlapping_frag.len(), read.name().unwrap().to_string());
            return None
        } else if overlapping_frag.len() == 0 {
            warn!("Warning: {} restriction fragments found for {} - skipped", 
                     overlapping_frag.len(), read.name().unwrap().to_string());
            return None
        } else{
            //let test = &overlapping_frag[0].val;
            return Some(overlapping_frag[0].val.clone())
        }
    } else{
        warn!("Warning: No restriction fragments found for {} - skipped", 
                     read.name().unwrap().to_string());
        return None
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

fn create_output_handlers(output_dir: &PathBuf, base_name: &str, all_output: bool, sam_output: bool) -> Result<OutputHandlers, Box<dyn Error>> {
    let valid_file = output_dir.join(format!("{}.validPairs", base_name));
    let valid = BufWriter::new(OpenOptions::new().create(true).write(true).truncate(true).open(valid_file)?);

    let de = if all_output {
        let de_file = output_dir.join(format!("{}.DEPairs", base_name));
        Some(BufWriter::new(OpenOptions::new().create(true).write(true).truncate(true).open(de_file)?))
    } else {
        None
    };

    let re = if all_output {
        let re_file = output_dir.join(format!("{}.REPairs", base_name));
        Some(BufWriter::new(OpenOptions::new().create(true).write(true).truncate(true).open(re_file)?))
    } else {
        None
    };

    let sc = if all_output {
        let sc_file = output_dir.join(format!("{}.SCPairs", base_name));
        Some(BufWriter::new(OpenOptions::new().create(true).write(true).truncate(true).open(sc_file)?))
    } else {
        None
    };

    let dump = if all_output {
        let dump_file = output_dir.join(format!("{}.DumpPairs", base_name));
        Some(BufWriter::new(OpenOptions::new().create(true).write(true).truncate(true).open(dump_file)?))
    } else {
        None
    };

    let single = if all_output {
        let single_file = output_dir.join(format!("{}.SinglePairs", base_name));
        Some(BufWriter::new(OpenOptions::new().create(true).write(true).truncate(true).open(single_file)?))
    } else {
        None
    };

    let filt = if all_output {
        let filt_file = output_dir.join(format!("{}.FiltPairs", base_name));
        Some(BufWriter::new(OpenOptions::new().create(true).write(true).truncate(true).open(filt_file)?))
    } else {
        None
    };

    let sam = if sam_output {
        let sam_file = output_dir.join(format!("{}_interaction.bam", base_name));
        let file = OpenOptions::new().create(true).write(true).truncate(true).open(sam_file)?;
        Some(file)
    } else {
        None
    };

    Ok(OutputHandlers {
        valid,
        de,
        re,
        sc,
        dump,
        single,
        filt,
        sam,
    })
}

fn write_valid_pair(
    handler: &mut BufWriter<File>,
    read1: &bam::Record,
    _read2: &bam::Record,
    r1_chrom: &str,
    r2_chrom: &str,
    r1_pos: usize,
    r2_pos: usize,
    r1_strand: &str,
    r2_strand: &str,
    dist: Option<u64>,
    r1_fragname: &str,
    r2_fragname: &str,
    r1_mapq: u8,
    r2_mapq: u8,
    htag: &str,
) -> Result<(), Box<dyn Error>> {
    writeln!(
        handler,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        read1.name().map(|n| n.to_string()).unwrap_or_else(|| "Unknown".to_string()),
        r1_chrom,
        r1_pos, 
        r1_strand,
        r2_chrom,
        r2_pos, 
        r2_strand,
        dist.map_or_else(|| "*".to_string(), |d| d.to_string()),
        r1_fragname,
        r2_fragname,
        r1_mapq,
        r2_mapq,
        htag
    )?;
    Ok(())
}

fn write_single_pair(
    handler: &mut BufWriter<File>,
    read: &bam::Record,
    chrom: &str,
    pos: usize,
    strand: &str,
    fragname: &str,
    mapq: u8,
) -> Result<(), Box<dyn Error>> {
    writeln!(
        handler,
        "{}\t{}\t{}\t{}\t*\t*\t*\t*\t{}\t*\t{}\t*",
        read.name().map(|n| n.to_string()).unwrap_or_else(|| "Unknown".to_string()),
        chrom,
        pos + 1, // Convert to 1-based
        strand,
        fragname,
        mapq
    )?;
    Ok(())
}


fn write_statistics(stats: &Statistics, output_dir: &PathBuf, base_name: &str, gtag: Option<&String>) -> Result<(), Box<dyn Error>> {
    let stat_file = output_dir.join(format!("{}.RSstat", base_name));
    let mut stat_writer = BufWriter::new(OpenOptions::new().create(true).write(true).truncate(true).open(stat_file)?);
    
    writeln!(stat_writer, "## Hi-C processing")?;
    writeln!(stat_writer, "Valid_interaction_pairs\t{}", stats.valid_counter)?;
    writeln!(stat_writer, "Valid_interaction_pairs_FF\t{}", stats.valid_counter_ff)?;
    writeln!(stat_writer, "Valid_interaction_pairs_RR\t{}", stats.valid_counter_rr)?;
    writeln!(stat_writer, "Valid_interaction_pairs_RF\t{}", stats.valid_counter_rf)?;
    writeln!(stat_writer, "Valid_interaction_pairs_FR\t{}", stats.valid_counter_fr)?;
    writeln!(stat_writer, "Dangling_end_pairs\t{}", stats.de_counter)?;
    writeln!(stat_writer, "Religation_pairs\t{}", stats.re_counter)?;
    writeln!(stat_writer, "Self_Cycle_pairs\t{}", stats.sc_counter)?;
    writeln!(stat_writer, "Single-end_pairs\t{}", stats.single_counter)?;
    writeln!(stat_writer, "Filtered_pairs\t{}", stats.filt_counter)?;
    writeln!(stat_writer, "Dumped_pairs\t{}", stats.dump_counter)?;

    if let Some(_gtag) = gtag {
        writeln!(stat_writer, "## ======================================")?;
        writeln!(stat_writer, "## Allele specific information")?;
        writeln!(stat_writer, "Valid_pairs_from_ref_genome_(1-1)\t{}", stats.g1g1_ascounter)?;
        writeln!(stat_writer, "Valid_pairs_from_ref_genome_with_one_unassigned_mate_(0-1/1-0)\t{}", 
                 stats.ug1_ascounter + stats.g1u_ascounter)?;
        writeln!(stat_writer, "Valid_pairs_from_alt_genome_(2-2)\t{}", stats.g2g2_ascounter)?;
        writeln!(stat_writer, "Valid_pairs_from_alt_genome_with_one_unassigned_mate_(0-2/2-0)\t{}", 
                 stats.ug2_ascounter + stats.g2u_ascounter)?;
        writeln!(stat_writer, "Valid_pairs_from_alt_and_ref_genome_(1-2/2-1)\t{}", 
                 stats.g1g2_ascounter + stats.g2g1_ascounter)?;
        writeln!(stat_writer, "Valid_pairs_with_both_unassigned_mated_(0-0)\t{}", stats.uu_ascounter)?;
        writeln!(stat_writer, "Valid_pairs_with_at_least_one_conflicting_mate_(3-)\t{}", stats.cf_ascounter)?;
    }

    Ok(())
}


fn process_read_pair(
    r1: &bam::Record,
    r1_chrom: Option<&str>,
    r1_resfrag: Option<&BED<6>>,
    r2: &bam::Record,
    r2_chrom: Option<&str>,
    r2_resfrag: Option<&BED<6>>,
    handlers: &mut OutputHandlers,
    stats: &mut Statistics,
    cli: Cli,
) -> Result<(), Box<dyn Error>> {
    let interaction_type = get_interaction_type(
        r1, 
        r1_chrom.unwrap_or("*"), 
        r1_resfrag.cloned(),
        r2, 
        r2_chrom.unwrap_or("*"), 
        r2_resfrag.cloned(),
        cli.verbose
    );
    
    let dist = get_pe_fragment_size(r1, r2, 
        r1_resfrag.cloned(),
        r2_resfrag.cloned(),
        interaction_type.unwrap_or("UNKNOWN")
    );
    
    let cdist = get_cis_distance(r1, r2);
    
    // Apply filters
    let mut final_interaction_type = interaction_type;
    
    // Check insert size criteria
    if let Some(distance) = dist {
        if let Some(min_size) = cli.min_insert_size {
            if distance < min_size {
                final_interaction_type = Some("FILT");
            }
        }
        if let Some(max_size) = cli.max_insert_size {
            if distance > max_size {
                final_interaction_type = Some("FILT");
            }
        }
    }
    
    // Check distance criteria for valid interactions
    if final_interaction_type == Some("VI") {
        if let Some(min_dist) = cli.min_cis_dist {
            if let Some(cis_dist) = cdist {
                if (cis_dist as u64) < min_dist {
                    final_interaction_type = Some("FILT");
                }
            }
        }
    }
    
    // Update statistics and write output
    match final_interaction_type {
        Some("VI") => {
            stats.valid_counter += 1;
            let valid_type = get_valid_orientation(r1, r2);
            match valid_type {
                "FF" => stats.valid_counter_ff += 1,
                "RR" => stats.valid_counter_rr += 1,
                "FR" => stats.valid_counter_fr += 1,
                "RF" => stats.valid_counter_rf += 1,
                _ => {}
            }
            
            // Handle allele specific counting if gtag is provided
            if let Some(ref _gtag) = cli.gtag {
                // This would require implementing get_read_tag function
                // For now, we'll skip allele specific counting
            }
            
            write_output_pair(r1, r2, r1_chrom, r2_chrom, r1_resfrag, r2_resfrag, 
                            dist, &mut handlers.valid, cli.gtag.as_deref())?;
        }
        Some("DE") => {
            stats.de_counter += 1;
            if let Some(ref mut handler) = handlers.de {
                write_output_pair(r1, r2, r1_chrom, r2_chrom, r1_resfrag, r2_resfrag, 
                                dist, handler, cli.gtag.as_deref())?;
            }
        }
        Some("RE") => {
            stats.re_counter += 1;
            if let Some(ref mut handler) = handlers.re {
                write_output_pair(r1, r2, r1_chrom, r2_chrom, r1_resfrag, r2_resfrag, 
                                dist, handler, cli.gtag.as_deref())?;
            }
        }
        Some("SC") => {
            stats.sc_counter += 1;
            if let Some(ref mut handler) = handlers.sc {
                write_output_pair(r1, r2, r1_chrom, r2_chrom, r1_resfrag, r2_resfrag, 
                                dist, handler, cli.gtag.as_deref())?;
            }
        }
        Some("SI") => {
            stats.single_counter += 1;
            if let Some(ref mut handler) = handlers.single {
                write_single_output(r1, r2, r1_chrom, r2_chrom, r1_resfrag, r2_resfrag, handler)?;
            }
        }
        Some("FILT") => {
            stats.filt_counter += 1;
            if let Some(ref mut handler) = handlers.filt {
                write_output_pair(r1, r2, r1_chrom, r2_chrom, r1_resfrag, r2_resfrag, 
                                dist, handler, cli.gtag.as_deref())?;
            }
        }
        _ => {
            stats.dump_counter += 1;
            if let Some(ref mut handler) = handlers.dump {
                write_output_pair(r1, r2, r1_chrom, r2_chrom, r1_resfrag, r2_resfrag, 
                                dist, handler, cli.gtag.as_deref())?;
            }
        }
    }

        Ok(())
}


fn write_output_pair(
    r1: &bam::Record,
    r2: &bam::Record,
    r1_chrom: Option<&str>,
    r2_chrom: Option<&str>,
    r1_resfrag: Option<&BED<6>>,
    r2_resfrag: Option<&BED<6>>,
    dist: Option<u64>,
    handler: &mut BufWriter<File>,
    gtag: Option<&str>,
) -> Result<(), Box<dyn Error>> {
    if !r1.flags().is_unmapped() && !r2.flags().is_unmapped() {
        // Get ordered reads
        if let Some((or1, or2)) = get_ordered_reads(r1, r2) {
            let or1_chrom = if ptr::eq(or1, r1) { r1_chrom } else { r2_chrom };
            let or2_chrom = if ptr::eq(or1, r1) { r2_chrom } else { r1_chrom };
            
            let or1_resfrag = if ptr::eq(or1, r1) { r1_resfrag } else { r2_resfrag };
            let or2_resfrag = if ptr::eq(or1, r1) { r2_resfrag } else { r1_resfrag };
            
            let or1_pos = get_read_pos(or1, "start").unwrap_or(0);
            let or2_pos = get_read_pos(or2, "start").unwrap_or(0);
            let or1_strand = get_read_strand(or1);
            let or2_strand = get_read_strand(or2);
            
            let or1_fragname = or1_resfrag.map(|f| f.name().unwrap().to_string()).unwrap_or_else(|| "None".to_string());
            let or2_fragname = or2_resfrag.map(|f| f.name().unwrap().to_string()).unwrap_or_else(|| "None".to_string());
            
            let htag = if gtag.is_some() { "0-0" } else { "" };
            
            write_valid_pair(
                handler,
                or1, or2,
                or1_chrom.unwrap_or("*"),
                or2_chrom.unwrap_or("*"),
                or1_pos, or2_pos,
                or1_strand, or2_strand,
                dist,
                &or1_fragname, &or2_fragname,
                or1.mapping_quality().map(|q| q.get()).unwrap_or(0),
                or2.mapping_quality().map(|q| q.get()).unwrap_or(0),
                htag,
            )?;
        }
    } else if r2.flags().is_unmapped() && !r1.flags().is_unmapped() {
        let r1_pos = get_read_pos(r1, "start").unwrap_or(0);
        let r1_strand = get_read_strand(r1);
        let r1_fragname = r1_resfrag.map(|f| f.chrom().to_string()).unwrap_or_else(|| "None".to_string());
        
        write_single_pair(
            handler,
            r1,
            r1_chrom.unwrap_or("*"),
            r1_pos,
            r1_strand,
            &r1_fragname,
            r1.mapping_quality().map(|q| q.get()).unwrap_or(0),
        )?;
    } else if r1.flags().is_unmapped() && !r2.flags().is_unmapped() {
        let r2_pos = get_read_pos(r2, "start").unwrap_or(0);
        let r2_strand = get_read_strand(r2);
        let r2_fragname = r2_resfrag.map(|f| f.chrom().to_string()).unwrap_or_else(|| "None".to_string());
        
        write_single_pair(
            handler,
            r2,
            r2_chrom.unwrap_or("*"),
            r2_pos,
            r2_strand,
            &r2_fragname,
            r2.mapping_quality().map(|q| q.get()).unwrap_or(0),
        )?;
    }
    
    Ok(())
}

fn write_single_output(
    r1: &bam::Record,
    r2: &bam::Record,
    r1_chrom: Option<&str>,
    r2_chrom: Option<&str>,
    r1_resfrag: Option<&BED<6>>,
    r2_resfrag: Option<&BED<6>>,
    handler: &mut BufWriter<File>,
) -> Result<(), Box<dyn Error>> {
    if !r1.flags().is_unmapped() {
        let r1_pos = get_read_pos(r1, "start").unwrap_or(0);
        let r1_strand = get_read_strand(r1);
        let r1_fragname = r1_resfrag.map(|f| f.chrom().to_string()).unwrap_or_else(|| "None".to_string());
        
        write_single_pair(
            handler,
            r1,
            r1_chrom.unwrap_or("*"),
            r1_pos,
            r1_strand,
            &r1_fragname,
            r1.mapping_quality().map(|q| q.get()).unwrap_or(0),
        )?;
    }
    
    if !r2.flags().is_unmapped() {
        let r2_pos = get_read_pos(r2, "start").unwrap_or(0);
        let r2_strand = get_read_strand(r2);
        let r2_fragname = r2_resfrag.map(|f| f.chrom().to_string()).unwrap_or_else(|| "None".to_string());
        
        write_single_pair(
            handler,
            r2,
            r2_chrom.unwrap_or("*"),
            r2_pos,
            r2_strand,
            &r2_fragname,
            r2.mapping_quality().map(|q| q.get()).unwrap_or(0),
        )?;
    }

        Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    
    let cli = Cli::parse();
    
    // Set up output directory
    let output_dir = cli.out_dir.clone().unwrap_or_else(|| PathBuf::from("."));
    std::fs::create_dir_all(&output_dir)?;
    
    // Get base name for output files
    let bam_path = PathBuf::from(&cli.bam);
    let base_name = bam_path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output");
    
    if cli.verbose {
        info!("## HiC-Pro Rust Implementation");
        info!("## mappedReadsFile= {}", cli.bam);
        info!("## fragmentFile= {}", cli.fragment_file);
        info!("## minInsertSize= {:?}", cli.min_insert_size);
        info!("## maxInsertSize= {:?}", cli.max_insert_size);
        info!("## minFragSize= {:?}", cli.min_frag_size);
        info!("## maxFragSize= {:?}", cli.max_frag_size);
        info!("## allOutput= {}", cli.all);
        info!("## SAM output= {}", cli.sam);
        info!("## verbose= {}", cli.verbose);
    }
    
    // Load restriction fragments
    if cli.verbose {
        info!("## Loading Restriction File Intervals {} ...", cli.fragment_file);
    }
    let bed_file_open = File::open(&cli.fragment_file)?;
    let bed_reader = Reader::new(bed_file_open, None);
    let bed_rec: Vec<BED<6>> = bed_reader.into_records::<BED<6>>()
        .map(|r| r.unwrap())
        .collect();
    
    // Filter fragments by size if specified
    let filtered_bed_rec: Vec<BED<6>> = bed_rec.into_iter()
        .filter(|bed| {
            let frag_len = bed.end() - bed.start();
            if let Some(min_size) = cli.min_frag_size {
                if frag_len < min_size { return false; }
            }
            if let Some(max_size) = cli.max_frag_size {
                if frag_len > max_size { return false; }
            }
            true
        })
        .collect();
    
    let bed_ladder = convert_vec_to_lapper(&filtered_bed_rec);
    
    // Create output handlers
    let mut handlers = create_output_handlers(&output_dir, base_name, cli.all, cli.sam)?;
    
    // Initialize statistics
    let mut stats = Statistics::default();
    
    // Open BAM file
    if cli.verbose {
        info!("## Opening BAM file {} ...", cli.bam);
    }
    let mut reader = bam::io::reader::Builder::default().build_from_path(&cli.bam)?;
    let headers = reader.read_header()?;
    
    if cli.verbose {
        info!("## Classifying Interactions ...");
    }
    
    // Process reads
    let mut r1: Option<bam::Record> = None;
    let mut r1_chrom: Option<String> = None;
    let mut r1_resfrag: Option<BED<6>> = None;
    
    for result in reader.records() {
        let record = result?;
        stats.reads_counter += 1;
        
        if record.flags().is_first_segment() {
            r1 = Some(record);
            if let Some(ref r1_ref) = r1 {
                if !r1_ref.flags().is_unmapped() {
                    if let Some(result) = r1_ref.reference_sequence(&headers) {
                        let (name_bytes, _) = result?;
                        r1_chrom = Some(std::str::from_utf8(name_bytes)?.to_string());
                        r1_resfrag = get_overlapping_restriction_fragment(&bed_ladder, r1_chrom.as_ref().unwrap(), r1_ref);
                    } else {
                        r1_chrom = None;
                        r1_resfrag = None;
                    }
                } else {
                    r1_chrom = None;
                    r1_resfrag = None;
                }
            }
        } else if record.flags().is_last_segment() {
            let r2 = record;
            let (r2_chrom, r2_resfrag) = if !r2.flags().is_unmapped() {
                if let Some(result) = r2.reference_sequence(&headers) {
                    let (name_bytes, _) = result?;
                    let r2_chrom = std::str::from_utf8(name_bytes)?.to_string();
                    let res_frag = get_overlapping_restriction_fragment(&bed_ladder, &r2_chrom, &r2);
                    (Some(r2_chrom), res_frag)
                } else {
                    (None, None)
                }
            } else {
                (None, None)
            };
            
            // Process the read pair
            if let Some(r1_ref) = r1.take() {
                process_read_pair(
                    &r1_ref,
                    r1_chrom.as_deref(),
                    r1_resfrag.as_ref(),
                    &r2,
                    r2_chrom.as_deref(),
                    r2_resfrag.as_ref(),
                    &mut handlers,
                    &mut stats,
                    cli.clone(),
                )?;
            }
            
            // Reset for next pair
            r1 = None;
            r1_chrom = None;
            r1_resfrag = None;
        }
        
        if stats.reads_counter % 100000 == 0 && cli.verbose {
            info!("## {}", stats.reads_counter);
        }
    }
    
    // Write statistics
    write_statistics(&stats, &output_dir, base_name, cli.gtag.as_ref())?;
    
    if cli.verbose {
        info!("## Processing complete!");
        info!("## Total reads processed: {}", stats.reads_counter);
        info!("## Valid interactions: {}", stats.valid_counter);
    }
    
    Ok(())
}

