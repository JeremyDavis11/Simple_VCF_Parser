# VCF Parser and Filter

A basic command-line tool for parsing and filtering variant call format (VCF) files by quality score, read depth, missing genotype rate, and minor allele frequency.

## Background

VCF files are the standard format for storing genetic variant calls (SNPs, indels, etc.) from sequencing data. Raw VCF files often contain low-confidence variants that can introduce noise into downstream analyses. Filtering by quality metrics is a critical preprocessing step to retain only reliable variant calls.

This tool supports filtering on four common criteria:

- **Quality score (QUAL)** — the Phred-scaled confidence that a variant exists at a site. Higher is better; low quality scores indicate uncertain calls.
- **Read depth (DP)** — the number of reads covering a site. Low depth reduces confidence in the genotype call.
- **Missing genotype rate** — the proportion of samples with no genotype call (`./.`) at a site. High proportions of missing data can bias population-level analyses.
- **Minor allele frequency (MAF)** — the frequency of the less common allele across samples. Very rare variants are often sequencing errors rather than true variants. This assumes that the input VCF file is biallelic. 

## Dependencies

- Python 3.x
- pandas

Install dependencies with:

```
pip install -r requirements.txt
```

## Usage

```
python VCF_parser.py input.vcf [options]
```

### Arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `input` | positional | — | Path to input VCF file |
| `--min_qual` | int | 30 | Minimum quality score |
| `--min_depth` | int | 10 | Minimum read depth |
| `--max_missing` | float | 0.2 | Maximum missing genotype rate (0–1) |
| `--min_maf` | float | 0.05 | Minimum minor allele frequency (0–1) |
| `--output` | str | filtered.vcf | Path to output VCF file |

### Examples

Run with default filters:
```
python VCF_parser.py input.vcf
```

Run with custom filters:
```
python VCF_parser.py input.vcf --min_qual 40 --min_depth 15 --max_missing 0.1 --min_maf 0.1 --output output.vcf
```

### Output

The filtered VCF is written to the specified output file, preserving the original metadata header. A summary is printed to the terminal:

```
Variants before filtering: 10
Variants after filtering: 6
Variants removed: 4
```

## Test Data

A small test VCF (`test.vcf`) with 10 variants across 3 samples is included. Using default filter settings, 6 variants are expected to pass.

```
python VCF_parser.py test.vcf
```

## Limitations
### biallelic VCF
This parser assumes that the input VCF file is a standard biallelic, with a reference and alternative allele and not multiallelic. 

### Depth parsed from INFO field
some VCFs store depth information in the FORMAT/sample fields rather than the INFO field, which this parser does not catch.

### Genotype assumes diploid
The genotype parsing logic splits on ```/``` and expects exactly two alleles. polyploid calls will not be properly parsed.
=======
# Simple_VCF_Parser
Parse a biallelic VCF file and filter variants based on various criteria defined at the command line.
