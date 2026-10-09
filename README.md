# circure
**circure** is a script for validating the circularity of long-read assembled sequences.

Circularity is inferred by mapping (primary) long-reads back to the assembly with **≥95% identity** and **≥99% of their length** aligned, and requiring, at a minimum, that the reads:
- fully cover the assembled sequence  
- map continuously across the contig start and end  
  (i.e., the artificial breakpoint introduced by the assembler)

### [More explanation in  section ***to be added***]

## Usage
**circure** accepts a subset of the [PAF format](https://github.com/lh3/minimap2/blob/master/PAF.md) as input, which must be generated using [minimap2](https://github.com/lh3/minimap2). **circure** accepts a subset of the PAF format as input, which must be generated using minimap2. If using accurate long reads (Q20+), it is advised to use the `-x lr:hq` preset.


### Example dataset
The example dataset (`PAS01578.dorado2.0.0.bmdna_r10.4.1_e8.2_400bps_6.0.0_hac.sim-200000.fastq.gz`) is a subsample of 200,000 Oxford Nanopore reads from the [ZymoBIOMICS HMW DNA Standard](https://zymoresearch.eu/products/zymobiomics-hmw-dna-standard), obtained from [MicroBench](https://github.com/Kirk3gaard/MicroBench) (ENA project [PRJEB85558](https://www.ebi.ac.uk/ena/browser/view/PRJEB85558)). The .fastq can be downloaded from https://zenodo.org/uploads/22944851. 
The sample contains four plasmids, one *Escherichia coli* plasmid (110,007 bp) and three *Staphylococcus aureus* plasmids (6,337 bp, 2,993 bp and 2,216 bp), which are provided in `zymohmw_plasmids.fasta`. 

### Example command:
```bash
minimap2 -cx lr:hq --secondary=no -t {threads} -o {output.paf} {input.fasta} {input.reads.fastq}
```
Using the `--secondary=no` flag ensures that minimap2 only reports **primary alignments**, which improves downstream performance (speed) when running the **circure** script.


Next, generate the input file for the **circure** script by extracting a subset of the PAF output from minimap2 with ```awk```:
```bash
awk 'BEGIN { FS="\t"; OFS="\t" } { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11 }' {input} > {output}
```
The tab-separated PAF input file must retain the `.paf` extension.
\
Finally, execute the **circure** script:
```bash
Rscript running_circure.R {input} {output}
```
The `{output}` argument must be given as a file. All dependencies required to run the **circure** script are listed in the circure.yml Conda environment file.


Results are output as a tab-separated file (examplified underneath from the zymoHMW run):

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

