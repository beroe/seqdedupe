use clap::Parser;
use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use chrono::Utc;

#[derive(Parser)]
#[command(name = "seqdedupe")]
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
}

#[derive(Debug, Clone)]
struct FastaRecord {
    header: String,
    sequence: String,
}

fn timestamp() -> String {
    Utc::now().format("%Y-%m-%d %H:%M:%S UTC").to_string()
}

fn parse_fasta(filename: &str) -> Result<Vec<FastaRecord>> {
    let file = File::open(filename)
        .with_context(|| format!("Failed to open file: {}", filename))?;
    let reader = BufReader::new(file);
    
    let mut records = Vec::new();
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
    
    Ok(records)
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

fn remove_exact_duplicates(records: Vec<FastaRecord>, is_dna: bool) -> Vec<FastaRecord> {
    let mut seen_sequences = HashSet::new();
    let mut unique_records = Vec::new();
    
    for record in records {
        let mut is_duplicate = false;
        
        if seen_sequences.contains(&record.sequence) {
            is_duplicate = true;
        }
        
        if is_dna && !is_duplicate {
            let rev_comp = reverse_complement(&record.sequence);
            if seen_sequences.contains(&rev_comp) {
                is_duplicate = true;
            }
        }
        
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

fn remove_substring_sequences(mut records: Vec<FastaRecord>, is_dna: bool) -> Vec<FastaRecord> {
    records.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));
    
    let mut kept_records: Vec<FastaRecord> = Vec::new();
    
    for i in 0..records.len() {
        let mut is_substring = false;
        let current_seq = &records[i].sequence;
        
        for j in 0..kept_records.len() {
            let longer_seq = &kept_records[j].sequence;
            
            if longer_seq.contains(current_seq) {
                is_substring = true;
                break;
            }
            
            if is_dna {
                let rev_comp = reverse_complement(current_seq);
                if longer_seq.contains(&rev_comp) {
                    is_substring = true;
                    break;
                }
            }
        }
        
        if !is_substring {
            kept_records.push(records[i].clone());
        }
    }
    
    kept_records
}

fn write_fasta(records: &[FastaRecord], output: Option<&str>) -> Result<()> {
    if let Some(filename) = output {
        let mut file = File::create(filename)
            .with_context(|| format!("Failed to create output file: {}", filename))?;
        
        for record in records {
            writeln!(file, "{}", record.header)?;
            writeln!(file, "{}", record.sequence)?;
        }
    } else {
        for record in records {
            println!("{}", record.header);
            println!("{}", record.sequence);
        }
    }
    
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    eprintln!("[{}] Reading FASTA file: {}", timestamp(), args.input);
    let records = parse_fasta(&args.input)?;
    eprintln!("[{}] Found {} sequences", timestamp(), records.len());
    
    eprintln!("[{}] Removing exact duplicates...", timestamp());
    let mut final_records = remove_exact_duplicates(records, args.dna);
    eprintln!("[{}] After removing exact duplicates: {} sequences", timestamp(), final_records.len());
    
    if args.substring {
        eprintln!("[{}] Removing substring sequences...", timestamp());
        final_records = remove_substring_sequences(final_records, args.dna);
        eprintln!("[{}] After removing substring sequences: {} sequences", timestamp(), final_records.len());
    }
    
    eprintln!("[{}] Final result: {} unique sequences", timestamp(), final_records.len());
    write_fasta(&final_records, args.output.as_deref())?;
    
    Ok(())
}
