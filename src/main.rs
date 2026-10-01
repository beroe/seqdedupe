use clap::Parser;
use anyhow::Result;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use chrono::Utc;
use rayon::prelude::*;

#[derive(Parser)]
#[command(name = "seqdedupe", version)]
#[command(about = "Remove duplicate and substring sequences from FASTA files")]
struct Args {
    #[arg(help = "Input FASTA file")]
    input: String,
    
    #[arg(short, long, help = "Output file (stdout if not specified)")]
    output: Option<String>,
    
    #[arg(short, long, help = "Treat as DNA sequences (check reverse complements)")]
    dna: bool,
    
    #[arg(short, long, help = "Remove substring sequences (slower for large files)")]
    substring: bool,
    
    #[arg(long, help = "Batch size for processing (default 10000)")]
    batch_size: Option<usize>,
    
    #[arg(long, help = "Number of CPU cores to use (default: half of available cores)")]
    cores: Option<usize>,
}

#[derive(Debug, Clone)]
struct FastaRecord {
    header: String,
    sequence: String,
}

fn timestamp() -> String {
    Utc::now().format("%Y-%m-%d %H:%M:%S UTC").to_string()
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'N' => 'N',
            '-' => '-',
            other => other,
        })
        .collect()
}

fn get_memory_usage() -> String {
    #[cfg(target_os = "linux")]
    {
        if let Ok(contents) = std::fs::read_to_string("/proc/self/status") {
            for line in contents.lines() {
                if line.starts_with("VmRSS:") {
                    return line.trim().to_string();
                }
            }
        }
    }
    "Memory usage: N/A".to_string()
}

// Streaming approach - process file in batches to reduce memory
fn process_streaming_duplicates(filename: &str, output_file: Option<&str>, is_dna: bool, batch_size: usize) -> Result<()> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    
    let mut seen_sequences = HashSet::new();
    let mut current_header = String::new();
    let mut current_sequence = String::new();
    let mut batch_records = Vec::new();
    let mut total_processed = 0;
    let mut total_kept = 0;
    
    // Setup output writer
    let mut output_writer: Box<dyn Write> = if let Some(output_path) = output_file {
        Box::new(File::create(output_path)?)
    } else {
        Box::new(std::io::stdout())
    };
    
    let mut line_count = 0;
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        line_count += 1;
        
        if line_count % 100000 == 0 {
            eprintln!("[{}] Processed {} lines, {} sequences kept, {}", 
                     timestamp(), line_count, total_kept, get_memory_usage());
        }
        
        if line.starts_with('>') {
            // Process previous record
            if !current_header.is_empty() {
                batch_records.push(FastaRecord {
                    header: current_header.clone(),
                    sequence: current_sequence.clone(),
                });
                
                // Process batch when full
                if batch_records.len() >= batch_size {
                    let (kept, processed) = process_batch(&mut batch_records, &mut seen_sequences, is_dna, &mut output_writer)?;
                    total_kept += kept;
                    total_processed += processed;
                    batch_records.clear();
                }
            }
            current_header = line.to_string();
            current_sequence.clear();
        } else if !line.is_empty() {
            current_sequence.push_str(line);
        }
    }
    
    // Process final record and batch
    if !current_header.is_empty() {
        batch_records.push(FastaRecord {
            header: current_header,
            sequence: current_sequence,
        });
    }
    
    if !batch_records.is_empty() {
        let (kept, processed) = process_batch(&mut batch_records, &mut seen_sequences, is_dna, &mut output_writer)?;
        total_kept += kept;
        total_processed += processed;
    }
    
    eprintln!("[{}] Final: {} sequences processed, {} unique kept", 
             timestamp(), total_processed, total_kept);
    
    Ok(())
}

fn process_batch(
    batch: &mut Vec<FastaRecord>, 
    seen_sequences: &mut HashSet<String>, 
    is_dna: bool,
    output_writer: &mut Box<dyn Write>
) -> Result<(usize, usize)> {
    let mut kept_count = 0;
    let processed_count = batch.len();
    
    for record in batch.drain(..) {
        let sequence = &record.sequence;
        let mut is_duplicate = false;
        
        // Check if sequence already seen
        if seen_sequences.contains(sequence) {
            is_duplicate = true;
        }
        
        // For DNA, also check reverse complement
        if is_dna && !is_duplicate {
            let rev_comp = reverse_complement(sequence);
            if seen_sequences.contains(&rev_comp) {
                is_duplicate = true;
            }
        }
        
        if !is_duplicate {
            seen_sequences.insert(sequence.clone());
            if is_dna {
                seen_sequences.insert(reverse_complement(sequence));
            }
            
            // Write immediately to reduce memory usage
            writeln!(output_writer, "{}", record.header)?;
            writeln!(output_writer, "{}", record.sequence)?;
            kept_count += 1;
        }
    }
    
    Ok((kept_count, processed_count))
}

fn remove_exact_duplicates(records: Vec<FastaRecord>, is_dna: bool) -> Vec<FastaRecord> {
    let mut seen_sequences = HashSet::new();
    let mut unique_records = Vec::new();
    
    for record in records {
        let is_duplicate = seen_sequences.contains(&record.sequence)
            || (is_dna && seen_sequences.contains(&reverse_complement(&record.sequence)));
        
        if !is_duplicate {
            seen_sequences.insert(record.sequence.clone());
            if is_dna {
                seen_sequences.insert(reverse_complement(&record.sequence));
            }
            unique_records.push(record);
        }
    }
    
    unique_records
}

// Parallel substring removal - processes batches of sequences in parallel
fn remove_substrings_parallel(input_file: &str, output_file: Option<&str>, is_dna: bool, num_cores: usize) -> Result<()> {
    eprintln!("[{}] Starting parallel substring removal using {} cores", timestamp(), num_cores);
    eprintln!("[{}] Warning: Substring removal on large files requires significant memory and time", timestamp());
    
    // Set rayon thread pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_cores)
        .build_global()
        .expect("Failed to set thread pool size");
    
    // Read all sequences (unavoidable for substring checking)
    let mut records = Vec::new();
    let file = File::open(input_file)?;
    let reader = BufReader::new(file);
    
    let mut current_header = String::new();
    let mut current_sequence = String::new();
    
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        
        if line.starts_with('>') {
            if !current_header.is_empty() {
                records.push(FastaRecord {
                    header: current_header.clone(),
                    sequence: current_sequence.clone(),
                });
            }
            current_header = line.to_string();
            current_sequence.clear();
        } else if !line.is_empty() {
            current_sequence.push_str(line);
        }
    }
    
    if !current_header.is_empty() {
        records.push(FastaRecord {
            header: current_header,
            sequence: current_sequence,
        });
    }
    
    eprintln!("[{}] Loaded {} sequences, {}", 
             timestamp(), records.len(), get_memory_usage());
    
    // Remove exact duplicates first so identical sequences are collapsed too
    let mut records = remove_exact_duplicates(records, is_dna);
    eprintln!("[{}] After removing exact duplicates: {} sequences", 
             timestamp(), records.len());
    
    // Sort by length (longest first)
    records.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));
    
    // Containment is transitive, so a sequence is redundant iff ANY strictly
    // longer sequence contains it (kept or not). Each check is independent,
    // so there is no shared state and no ordering race between threads.
    let total = records.len();
    let progress = AtomicUsize::new(0);
    let report_every = (total / 20).max(1);
    
    eprintln!("[{}] Checking {} sequences for substrings", timestamp(), total);
    
    let keep: Vec<bool> = records
        .par_iter()
        .map(|record| {
            let current_seq = &record.sequence;
            let rev_comp = if is_dna { Some(reverse_complement(current_seq)) } else { None };
            let n_longer = records.partition_point(|r| r.sequence.len() > current_seq.len());
            
            let is_substring = records[..n_longer].iter().any(|longer| {
                longer.sequence.contains(current_seq.as_str())
                    || rev_comp.as_ref().map_or(false, |rc| longer.sequence.contains(rc.as_str()))
            });
            
            let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
            if done % report_every == 0 {
                eprintln!("[{}] Checked {}/{} sequences, {}", 
                         timestamp(), done, total, get_memory_usage());
            }
            
            !is_substring
        })
        .collect();
    
    let final_vec: Vec<&FastaRecord> = records
        .iter()
        .zip(keep)
        .filter_map(|(record, k)| if k { Some(record) } else { None })
        .collect();
    
    // Write results
    eprintln!("[{}] Writing {} final sequences to output", timestamp(), final_vec.len());
    
    let mut output_writer: Box<dyn Write> = if let Some(output_path) = output_file {
        Box::new(File::create(output_path)?)
    } else {
        Box::new(std::io::stdout())
    };
    
    for record in final_vec.iter() {
        writeln!(output_writer, "{}", record.header)?;
        writeln!(output_writer, "{}", record.sequence)?;
    }
    
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    let batch_size = args.batch_size.unwrap_or(10000);
    
    // Determine number of cores to use
    let available_cores = num_cpus::get();
    let num_cores = args.cores.unwrap_or(available_cores / 2).max(1);
    
    eprintln!("[{}] Starting seqdedupe v{} with batch size {}", timestamp(), env!("CARGO_PKG_VERSION"), batch_size);
    eprintln!("[{}] Available cores: {}, using: {}", timestamp(), available_cores, num_cores);
    eprintln!("[{}] Initial {}", timestamp(), get_memory_usage());
    
    if args.substring {
        // For substring removal, use parallel processing
        remove_substrings_parallel(&args.input, args.output.as_deref(), args.dna, num_cores)?;
    } else {
        // Use streaming approach for exact duplicates only
        process_streaming_duplicates(&args.input, args.output.as_deref(), args.dna, batch_size)?;
        eprintln!("[{}] Note: only exact duplicates were removed. Use -s to also remove identical substrings.", timestamp());
    }
    
    eprintln!("[{}] Complete. Final {}", timestamp(), get_memory_usage());
    
    Ok(())
}