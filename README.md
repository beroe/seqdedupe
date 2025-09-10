# Sequence deduplication. 

For DNA files, it will only remove identical sequences (not substrings). 
It checks for identical forward and reverse-complement sequences in the file.

For Amino Acid files there is the option to remove identical substrings. 
Suggested workflow is to first remove exact identicals and save that as a
a separate file, then run the --substring flag on that file in parallel mode.

  - It will detect available cores and use half of them. 
  - Can be overridden with --cores flag

  For exact duplicates only (streaming):
   ```./target/release/seqdedupe-opt --dna large_file.fna -o deduped.fna```
  
  For substring removal (multithreaded):
  > Use all available cores

  ```release/seqdedupe-opt --dna --substring --cores 8 deduped.fna -o final.fna``` 

  > Or let it use half (default)

  ```release/seqdedupe-opt --dna --substring deduped.fna -o final.fna```

  Recommended Workflow:

  1. First pass: Exact duplicates only (fast, low memory)
  2. Second pass: Substring removal on smaller deduplicated file (parallel, faster)


