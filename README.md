# Introduction

**circure** is a script for validating the circularity of long-read assembled sequences.

Circularity is inferred by mapping (primary) long-reads back to the assembly with **≥95% identity** and **≥99% of their length** aligned, and requiring, at a minimum, that the reads:
- fully cover the assembled sequence  
- map continuously across the contig start and end  
  (i.e., the artificial breakpoint introduced by the assembler)

### [More explanation in  section ***to be added***]

# Usage
**circure** accepts a subset of the [PAF format](https://github.com/lh3/minimap2/blob/master/PAF.md) as input, which must be generated using [minimap2](https://github.com/lh3/minimap2). It is recommended to use accurate (Q20+) long reads, in which case the `-x lr:hq` preset is advised to use. 



### Example command:
```bash
minimap2 -cx lr:hq --secondary=no -t {threads} -o {output.paf} {input.fasta} {input.reads.fastq}
```
Using the `--secondary=no` flag ensures that minimap2 only reports **primary alignments**, which improves downstream performance (speed) when running the circure script.


Next, generate the input file for the **circure** script by extracting a subset of the PAF output from minimap2 with ```awk```:
```bash
awk 'OFS="\\t" {{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}}' {input} > {output}
```
The tab-separated PAF input file must retain the `.paf` extension.
\
Finally, execute the **circure** script:
```bash
Rscript running_circure.R {input} {output}
```
The `{output}` argument must be given as a file. All dependencies required to run the circure script are listed in the circure.yml Conda environment file.


The script will write the results as a tab-separated file.

### [Output explanation]

<!-- 
- **contig**: contig/seqeunce
- **file**: filename.
- **predcition**: TRUE if inferred/predicted circular
- **reads_mapping_over_ab**: # reads that map continuously across the contig the artificial breakpoint (ab)
- **reads_longer_than_contig_no_ab_split**: circularity has been passed but this contig has a read that is larger than the contig mapping to fully from end to end to the contig with no breaks.
- **reads_overhanging**: reads overhanging from the first or last base.
-->


<!-- 
# Annotated description of circure steps
-->

