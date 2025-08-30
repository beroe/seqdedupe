use clap::Parser;
use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::{Arc, Mutex};
use chrono::Utc;
use rayon::prelude::*;

#[derive(Parser)]
#[command(name = "seqdedupe")]
#[command(about = "Remove duplicate and substring sequences from FASTA files (memory optimized)")]
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
    
    eprintln!("[{}] Loaded {} sequences for substring checking, {}", 
             timestamp(), records.len(), get_memory_usage());
    
    // Sort by length (longest first)
    records.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));
    
    // Process sequences in parallel batches
    let batch_size = (records.len() / num_cores).max(1000); // At least 1000 per batch
    let kept_sequences = Arc::new(Mutex::new(HashSet::<String>::new()));
    let final_records = Arc::new(Mutex::new(Vec::new()));
    
    eprintln!("[{}] Processing in batches of {} sequences", timestamp(), batch_size);
    
    // Process records in chunks using parallel iteration
    records
        .par_chunks(batch_size)
        .enumerate()
        .for_each(|(chunk_idx, chunk)| {
            let mut local_kept: Vec<String> = Vec::new();
            
            for record in chunk {
                let current_seq = &record.sequence;
                let mut is_substring = false;
                
                // Check against globally kept sequences
                {
                    let kept_set = kept_sequences.lock().unwrap();
                    for kept_seq in kept_set.iter() {
                        if kept_seq != current_seq && kept_seq.contains(current_seq) {
                            is_substring = true;
                            break;
                        }
                        
                        if is_dna {
                            let rev_comp = reverse_complement(current_seq);
                            if kept_seq.contains(&rev_comp) {
                                is_substring = true;
                                break;
                            }
                        }
                    }
                }
                
                if !is_substring {
                    // Check against locally kept sequences in this batch
                    for local_seq in &local_kept {
                        if local_seq != current_seq && local_seq.contains(current_seq) {
                            is_substring = true;
                            break;
                        }
                        
                        if is_dna {
                            let rev_comp = reverse_complement(current_seq);
                            if local_seq.contains(&rev_comp) {
                                is_substring = true;
                                break;
                            }
                        }
                    }
                }
                
                if !is_substring {
                    local_kept.push(current_seq.clone());
                    
                    // Add to final results
                    {
                        let mut final_vec = final_records.lock().unwrap();
                        final_vec.push(record.clone());
                    }
                }
            }
            
            // Add local kept sequences to global set
            {
                let mut kept_set = kept_sequences.lock().unwrap();
                for seq in local_kept {
                    kept_set.insert(seq);
                }
            }
            
            if chunk_idx % 10 == 0 {
                let final_count = final_records.lock().unwrap().len();
                eprintln!("[{}] Processed chunk {}, {} sequences kept so far, {}", 
                         timestamp(), chunk_idx, final_count, get_memory_usage());
            }
        });
    
    // Write results
    let final_vec = final_records.lock().unwrap();
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
    
    eprintln!("[{}] Starting seqdedupe (optimized) with batch size {}", timestamp(), batch_size);
    eprintln!("[{}] Available cores: {}, using: {}", timestamp(), available_cores, num_cores);
    eprintln!("[{}] Initial {}", timestamp(), get_memory_usage());
    
    if args.substring {
        // For substring removal, use parallel processing
        remove_substrings_parallel(&args.input, args.output.as_deref(), args.dna, num_cores)?;
    } else {
        // Use streaming approach for exact duplicates only
        process_streaming_duplicates(&args.input, args.output.as_deref(), args.dna, batch_size)?;
    }
    
    eprintln!("[{}] Complete. Final {}", timestamp(), get_memory_usage());
    
    Ok(())
}