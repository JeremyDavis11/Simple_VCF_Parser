'''filter a VCF file by minimum quality score, minimum read depth, maximum missing genotypes, and minimum minor allele frequency'''

import pandas as pd
import argparse

def parse_vcf(filepath):
    '''
    parse a VCF file variant call information into a DataFrame
    put header metadata information into a list
    input filepath: path to VCF file to be parsed
    output meta: header lines collected into a list
    output df: dataframe with variant information
    '''
    meta = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('##'):
                meta.append(line)
            elif line.startswith('#CHROM'):
                cols = line.strip().lstrip('#').split('\t')
                break
        df = pd.read_csv(f, sep='\t', names=cols)
    return df, meta


def filter_vcf(df, min_qual, min_depth, max_missing, min_maf):
    '''
    filter variant calls based on defined parameters
    input df: dataframe with variant call information
    input min_qual: minimum quality score to filter for
    input min_depth: minimum read depth to filter for
    input max_missing: maximum proportion of missing genotypes before filtering varaint
    input min_maf: minimum minor allele frequency before filtering variant
    return df: new dataframe containing variants that passed filtering
    '''
    # Filter by QUAL column
    df = df[df['QUAL'] >= min_qual]
    
    # Extract DP from INFO field and filter by depth
    df = df.copy()
    df['DP'] = df['INFO'].str.extract(r'DP=(\d+)').astype(int)
    df = df[df['DP'] >= min_depth]

    # extract sample columns and filter by maximum missing genotypes
    sample_cols = df.columns[df.columns.get_loc('FORMAT') + 1:]
    df[sample_cols] = df[sample_cols].astype(str)
    genotypes = df[sample_cols].apply(lambda col: col.str.split(':').str[0])
    missing_rate = (genotypes == './.').sum(axis=1) / len(sample_cols)
    df = df[missing_rate <= max_missing]

    # calculate minor allele frequency nested function
    def calc_maf(row):
        '''
        calculate the minor allele frequency for a given variant
        input row: a row of data for a given variant
        return maf: the minor allele frequency for the variant
        '''
        allele_counts = {'0': 0, '1': 0}
        for sample in sample_cols:
            gt = row[sample].split(':')[0]
            if gt == './.':
                continue
            for allele in gt.replace('|', '/').split('/'):
                if allele in allele_counts:
                    allele_counts[allele] += 1
        total = sum(allele_counts.values())
        if total == 0:
            return 0
        maf = min(allele_counts.values()) / total
        return maf
    
    # filter by minimum maf frequency
    maf = df.apply(calc_maf, axis=1)
    df = df[maf >= min_maf]
    
    return df


def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='parse a VCF file and filter variants based on defined criteria')
    parser.add_argument('input', type=str, default='test.vcf', help='path to input VCF file')
    parser.add_argument('--min_qual', type=int, default=30, help='define the minimum quality score to filter for')
    parser.add_argument('--min_depth', type=int, default=10, help='define the minimum read depth to filter VCF for')
    parser.add_argument('--max_missing', type=float, default=0.2, help='define maximum missing genotypes per sample before filtering')
    parser.add_argument('--min_maf', type=float, default=0.05, help='define the minimum minor allele frequency for filtering' )
    parser.add_argument('--output', type=str, default='filtered.vcf', help='name the output filtered vcf file')
    args = parser.parse_args()

    # call functions
    df, meta = parse_vcf(args.input)
    filtered = filter_vcf(df, args.min_qual, args.min_depth, args.max_missing, args.min_maf)

    # print summary to terminal
    print(f"Variants before filtering: {len(df)}")
    print(f"Variants after filtering: {len(filtered)}")
    print(f"Variants removed: {len(df) - len(filtered)}")

    # write to output
    with open(args.output, 'w') as o:
        # write the metadata at the top of the file
        for entry in meta:
            o.write(entry)
        
        # write the column header line with # prefix
        o.write('#' + '\t'.join(filtered.columns.tolist()) + '\n')
        
        # write filtered lines from df
        filtered.to_csv(o, sep='\t', index=False, header=False)
    
        
if __name__ == '__main__':
    main()
