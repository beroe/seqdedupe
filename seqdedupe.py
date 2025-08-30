#!/usr/bin/env python3

import argparse
import sys
from datetime import datetime
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Remove duplicate and substring sequences from FASTA files"
    )
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("-o", "--output", help="Output file (stdout if not specified)")
    parser.add_argument("-d", "--dna", action="store_true", 
                       help="Treat as DNA sequences (check reverse complements)")
    parser.add_argument("-s", "--substring", action="store_true",
                       help="Remove substring sequences (slower for large files)")
    return parser.parse_args()

def remove_exact_duplicates(records, is_dna=False):
    """Remove exact duplicate sequences, including reverse complements for DNA."""
    seen_sequences = set()
    unique_records = []
    
    for record in records:
        sequence = str(record.seq)
        is_duplicate = False
        
        # Check if sequence already seen
        if sequence in seen_sequences:
            is_duplicate = True
        
        # For DNA, also check reverse complement
        if is_dna and not is_duplicate:
            rev_comp = str(Seq(sequence).reverse_complement())
            if rev_comp in seen_sequences:
                is_duplicate = True
        
        if not is_duplicate:
            seen_sequences.add(sequence)
            if is_dna:
                seen_sequences.add(str(Seq(sequence).reverse_complement()))
            unique_records.append(record)
    
    return unique_records

def remove_substring_sequences(records, is_dna=False):
    """Remove sequences that are substrings of other sequences."""
    # Sort by sequence length (longest first) for efficient substring checking
    sorted_records = sorted(records, key=lambda r: len(r.seq), reverse=True)
    kept_records = []
    
    for current_record in sorted_records:
        current_seq = str(current_record.seq)
        is_substring = False
        
        # Check if current sequence is a substring of any kept sequence
        for kept_record in kept_records:
            kept_seq = str(kept_record.seq)
            
            if current_seq in kept_seq:
                is_substring = True
                break
            
            # For DNA, also check if reverse complement is a substring
            if is_dna:
                rev_comp = str(Seq(current_seq).reverse_complement())
                if rev_comp in kept_seq:
                    is_substring = True
                    break
        
        if not is_substring:
            kept_records.append(current_record)
    
    return kept_records

def timestamp():
    """Return current timestamp as string."""
    return datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")

def write_fasta(records, output_file=None):
    """Write FASTA records to file or stdout."""
    if output_file:
        SeqIO.write(records, output_file, "fasta")
    else:
        SeqIO.write(records, sys.stdout, "fasta")

def main():
    args = parse_arguments()
    
    # Read FASTA file
    print(f"[{timestamp()}] Reading FASTA file: {args.input}", file=sys.stderr)
    records = list(SeqIO.parse(args.input, "fasta"))
    print(f"[{timestamp()}] Found {len(records)} sequences", file=sys.stderr)
    
    # Remove exact duplicates
    print(f"[{timestamp()}] Removing exact duplicates...", file=sys.stderr)
    final_records = remove_exact_duplicates(records, args.dna)
    print(f"[{timestamp()}] After removing exact duplicates: {len(final_records)} sequences", file=sys.stderr)
    
    # Remove substring sequences if requested
    if args.substring:
        print(f"[{timestamp()}] Removing substring sequences...", file=sys.stderr)
        final_records = remove_substring_sequences(final_records, args.dna)
        print(f"[{timestamp()}] After removing substring sequences: {len(final_records)} sequences", file=sys.stderr)
    
    print(f"[{timestamp()}] Final result: {len(final_records)} unique sequences", file=sys.stderr)
    
    # Write output
    write_fasta(final_records, args.output)

if __name__ == "__main__":
    main()