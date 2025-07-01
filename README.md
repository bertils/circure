# 🌀 circure

**circure** is a tool for detecting or validating the circularity of long-read assembled sequences.

Circularity is confirmed by mapping long reads back to the assembly and requiring, at a minimum, that reads:
- fully cover the assembled sequence
- map continuously across the contig start and end  
  (i.e., the artificial breakpoint introduced by the assembler)
