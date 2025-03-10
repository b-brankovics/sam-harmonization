# Development notes

Metadata aspects for future considerations:
- Comment lines could be used to store input file information
- Reference sequence can hold assembly ID information

## SAM format

> Format description: [https://samtools.github.io/hts-specs/SAMv1.pdf](https://samtools.github.io/hts-specs/SAMv1.pdf)

| Col | Field | Type   | Brief description |
| :-- | :--   | :-- | :-- |
| 1   | QNAME |	String | Query template NAME |
| 2   | FLAG  | Int    | bitwise FLAG |
| 3   | RNAME | String | References sequence NAME | 
| 4   | POS   | Int    | 1-based leftmost mapping POSition |
| 5   | MAPQ  | Int    | MAPping Quality |
| 6   | CIGAR | String | CIGAR String |
| 7   | RNEXT | String | Ref. name of the mate/next read |
| 8   | PNEXT | Int    | Position of the mate/next read |
| 9   | TLEN  | Int    | observed Template LENgth |
| 10  | SEQ   | String | segment SEQuence |
| 11  | QUAL  | String | ASCII of Phred-scaled base QUALity+33 |

Bitwise Flags

| Integer | Binary | Description (Paired Read Interpretation) |
| --:  | :--   | :-- |
| 1    | 000000000001 | template having multiple templates in sequencing (read is paired) |
| 2    | 000000000010 | each segment properly aligned according to the aligner (read mapped in proper pair) |
| 4    | 000000000100 | segment unmapped (read1 unmapped) |
| 8    | 000000001000 | next segment in the template unmapped (read2 unmapped) |
| 16   | 000000010000 | SEQ being reverse complemented (read1 reverse complemented) |
| 32   | 000000100000 | SEQ of the next segment in the template being reverse complemented (read2 reverse complemented) |
| 64   | 000001000000 | the first segment in the template (is read1) |
| 128  | 000010000000 | the last segment in the template (is read2) |
| 256  | 000100000000 | not primary alignment |
| 512  | 001000000000 | alignment fails quality checks |
| 1024 | 010000000000 | PCR or optical duplicate |
| 2048 | 100000000000 | supplementary alignment (e.g. aligner specific, could be a portion of a split read or a tied region) |

- unmapped: 4
- reverse: 16

The CIGAR operations are given in the following table (set ‘*’ if unavailable):

| Op | BAM | Consumes query | Consumes reference | Description |
| --- | --- | ---           | ---                | :--         |
| M  | 0   | yes            | yes                | alignment match (can be a sequence match or mismatch) |
| I  | 1   | yes            | no                 | insertion to the reference |
| D  | 2   | no             | yes                | deletion from the reference |
| N  | 3   | no             | yes                | skipped region from the reference |
| S  | 4   | yes            | no                 | soft clipping (clipped sequences present in SEQ) |
| H  | 5   | no*            | no                 | hard clipping (clipped sequences NOT present in SEQ) |
| P  | 6   | no             | no                 | padding (silent deletion from padded reference) |
| =  | 7   | yes            | yes                | sequence match |
| X  | 8   | yes            | yes                | sequence mismatch |

> __*__ `H` Does consume the query, but it does not consume the SEQ sequences. Sum of `H` operations plus the length of SEQ equeals the length of the original query sequence.

- “Consumes query” and “consumes reference” indicate whether the CIGAR operation causes the alignment to step along the query sequence and the reference sequence respectively.
- `H` can only be present as the first and/or last operation.
- `S` may only have `H` operations between them and the ends of the CIGAR string.
- For mRNA-to-genome alignment, an `N` operation represents an intron.  For other types of alignments, the interpretation of `N` is not defined.
- Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ. Sum of `H` operations plus the length of SEQ equeals the length of the original query sequence.

A CIGAR string is a set of length operation pairs: an integer marking the length of the opration. E.g. 150M means that 150 bp are aligned as match and/or mismatch.

## BLAST to SAM

BLAST can produce a SAM output, but it is **not a valid SAM file**. The main problem that `NM:i:<int>` is not correct it should be the edit distance gaps (indel) plus mismatch, but it is only mismatch.

```
# Build blast DB with taxIDs
makeblastdb -taxid_map test_map.txt  -in set1.fas -dbtype nucl -parse_seqids
# Run blast and print SAM format
blastn -query set2.fas -db set1.fas -outfmt 17 -parse_deflines >blast.sam
```

BLAST edit distance is not correct (NM:i:\d+). It only contains the [IDN] operations, not the mismatches!!

- `AS:i:1808` - integer alignment score, here 1808 (`score`)
- `NM:i:0` - integer edit distance, here 0
Plus non-standard tags:
- `EV:f:0` - float e-value, here 0 (`evalue`)
- `PI:f:96.55` - float percent identify, here 96.55
- `BS:f:701.049` - float bit-score, here 701.049 (`bitscore`)

`blast-sam.pl` is a wrapper to produce SAM format blast hits

`asn2sam.pl` is wrapper for `blast_formatter` to produce a SAM format. By using the ASN format, it is possible to run BLAST once and generate multiple formats from the results.

```
makeblastdb -taxid_map test_map.txt  -in set1.fas -dbtype nucl -parse_seqids
blastn -query set2.fas -db set1.fas -outfmt 11 -parse_deflines >blast.asn
asn2sam.pl blast.asn >blast.sam
```

### Portability of asn files

ASN.1 is not a stand alone file, it requires the reference index files.
The location of the reference file used for creating the blast DB is stored in the ASN.1 as a relative path to the working directory where the blast command was executed. For example:

```
subject database "dir1/ref.fas",
```

In this case it is in `dir1` folder and the original input fasta is called `ref.fas`.
So when ASN.1 is read by `blast_formatter` it will look for the indices in `dir1` folder relative to where `blast_formatter` is called, not relative to the ASN.1 file.
The actual `ref.fas` files is not needed only the index files generated by `makeblastdb`.

It specifically relies on finding `$subject-db\.ndb` and `$subject-db\.nin`. The rest of the files don't seem to be needed for this step. (Probably needed to run a BLAST on the db.)

## Alignment score

Alignment score is stored in the optional field in SAM with `AS:i:` tag as an integer.
Different programs calculate alignment scores in different ways.
To combine the results from different methods, scores have to be calculated in the same way.
Nucmer does not calculate alignment scores, so during conversion from delta a score is calculated using the default use by these tools (see below).

### Defalut scoring in exonerate

Based on test data, the following rewards and penalties can be estimated and the following equation.

| Type | weigth |
| --- | --- |
| identical | 5 |
| mismatch | -4 |
| gap extend | -4 |
| gap open | -12 |

```math
Score = 5 * (l_{alignment} -  d_{edit}) - 4 * d_{edit} - 8 * n_{gaps} = 5 * l_{alignment} - 9 * d_{edit} - 8 * n_{gaps}
```
Where $`l_{alignment}`$ is the alignment length, $`d_{edit}`$ is the edit distance (`NM:i:<int>`) and $`n_{gaps}`$ is the number of gaps (number of [IDN] in the CIGAR)

### Defalut scoring in blast

Based on test data, the following rewards and penalties can be estimated and the following equation.

| Type | weight |
| --- | --- |
| identical | 1 |
| mismatch | -2 |
| gap extend | -2 |
| gap open | -3 |

```math
Score = (l_{alignment} -  d_{edit}) - 2 * d_{edit} - 1 * n_{gaps} = l_{alignment} - 3 * d_{edit} - 1 * n_{gaps}
```

Where $`l_{alignment}`$ is the alignment length, $`d_{edit}`$ is the edit distance (`NM:i:<int>`) and $`n_{gaps}`$ is the number of gaps (number of [IDN] in the CIGAR)

This does not match the scores for inserts for some reason. Could not identify what goes wrong there. It does not seem to use a linear model for introns. Maybe there is a different scoring system that gets rounded to integers.
