# 🌀 circure

**circure** is a script for detecting or validating the circularity of long-read assembled sequences.

Circularity is confirmed by mapping (primary) long-reads back to the assembly with **≥95% identity** and **≥99% of their length** aligned, and requiring, at a minimum, that the reads:
- fully cover the assembled sequence  
- map continuously across the contig start and end  
  (i.e., the artificial breakpoint introduced by the assembler)
