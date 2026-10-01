# Sequence deduplication. 

By default it only removes identical sequences (not substrings), streaming
the file with low memory use. With `--dna` it also checks for identical
forward and reverse-complement sequences in the file.

Use `-s` / `--substring` to also remove sequences that are identical
substrings of longer ones (exact duplicates are always removed first).
Suggested workflow for large files is to first remove exact identicals and save
that as a separate file, then run the --substring flag on that file in parallel mode.

## Building

Precompiled binaries for Linux, macOS and Windows are available under
**Releases** on the right side of the GitHub repository page
(https://github.com/beroe/seqdedupe/releases). Download the archive for
your platform, unpack it, and run `seqdedupe`. No Rust needed.

To build from source instead, you need Rust (https://rustup.rs). Build once after cloning (and again after pulling changes):

   ```cargo build --release```

The binary is written to `./target/release/seqdedupe`. Alternatively, install it onto your PATH:

   ```cargo install --git https://github.com/beroe/seqdedupe```

## Usage

  - It will detect available cores and use half of them. 
  - Can be overridden with --cores flag
  - `seqdedupe -v` (or `--version`) prints the version and help

  For exact duplicates only (streaming):
   ```./target/release/seqdedupe --dna large_file.fna -o deduped.fna```
  
  For substring removal (multithreaded):
  > Use 8 cores

  ```./target/release/seqdedupe --dna --substring --cores 8 deduped.fna -o final.fna``` 

  > Or let it use half (default)

  ```./target/release/seqdedupe --dna --substring deduped.fna -o final.fna```

  Recommended Workflow:

  1. First pass: Exact duplicates only (fast, low memory)
  2. Second pass: Substring removal on smaller deduplicated file (parallel, faster)


