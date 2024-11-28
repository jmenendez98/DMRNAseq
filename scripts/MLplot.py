import argparse
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from typing import List, Tuple, Generator
from collections import defaultdict
import random

def sample_reads(bam_file: str, num_reads: int, min_mapq: int) -> List[pysam.AlignedSegment]:
    """
    Sample aligned reads from a BAM file with MAPQ > 10.
    
    Args:
        bam_file (str): Path to the input BAM file
        num_reads (int): Number of reads to sample
    
    Returns:
        List of sampled reads
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Filter for aligned reads with MAPQ > min_mapq
        def is_high_quality_read(read):
            return (not read.is_unmapped) and (read.mapping_quality > min_mapq)
        
        # Get all high-quality reads
        high_quality_reads = [read for read in bam.fetch() if is_high_quality_read(read)]
        
        # If no sampling specified, return all high-quality reads
        if num_reads is None or num_reads >= len(high_quality_reads):
            return high_quality_reads
        
        # Reservoir sampling for high-quality reads
        sampled_reads = []
        for i, read in enumerate(high_quality_reads):
            if i < num_reads:
                sampled_reads.append(read)
            else:
                # Randomly replace reads with decreasing probability
                j = random.randint(0, i)
                if j < num_reads:
                    sampled_reads[j] = read
        
        return sampled_reads

def create_dataframes(reads: List[pysam.AlignedSegment], mod_code: str) -> Tuple[pd.DataFrame, str, int]:
    """
    Create dataframe for a specific modification count.
    
    Args:
        reads (List[pysam.AlignedSegment]): List of sampled reads
        mod_code (str): Modification code to filter (e.g., 'a' or 'm')
    
    Returns:
        Tuple of DataFrame, modification type, and total read length
    """
    # Use more efficient defaultdict with default factory
    mod_dict = defaultdict(int)
    sequenced_bp = 0

    for read in reads:
        if read.modified_bases_forward:
            sequenced_bp += len(read.query_sequence)
            for mod, mod_sites in read.modified_bases_forward.items():
                mod_type = mod[2]
                if mod_type == mod_code:
                    for mod_site in mod_sites:
                        MLscore = mod_site[1]
                        mod_dict[MLscore] += 1

    # Create dataframe
    df = pd.DataFrame.from_dict(mod_dict, orient='index', columns=['count']).sort_index()
    
    return df, mod_code, sequenced_bp

def plot_histogram(dataframe: pd.DataFrame, mod: str, sequenced_bp: int, prefix: str, 
                   normalize_by: int = 1000, ml_threshold: float = None) -> None:
    """
    Plot normalized histogram of modification ML scores.
    
    Args:
        dataframe (pd.DataFrame): Dataframe with modification counts
        mod (str): Modification type
        sequenced_bp (int): Total number of sequenced bases
        prefix (str): Output file prefix
        normalize_by (int): Base pair length for normalization (default: 1000 for kb)
        ml_threshold (float, optional): Custom ml_threshold for vertical line
    """
    # More precise normalization
    scaling_factor = sequenced_bp / normalize_by
    dataframe['normalized_count'] = dataframe['count'] / scaling_factor

    plt.figure(figsize=(10, 6))
    plt.hist(dataframe.index, 
             weights=dataframe['normalized_count'], 
             bins=128, 
             color='skyblue', 
             range=(0, 256),
             edgecolor='black')
    
    plt.xlabel('ML Score')
    plt.ylabel(f'Count per {normalize_by} bp')
    plt.title(f"ML Scores for mod code: {mod}")

    # Add ml_threshold line if specified
    if ml_threshold is not None:
        plt.axvline(x=ml_threshold, color='red', linestyle='-', linewidth=2, 
                    label=f'Threshold (ML: {ml_threshold})')

        # Calculate percentage of data below threshold
        below_threshold = dataframe.loc[dataframe.index <= ml_threshold, 'normalized_count'].sum()
        total = dataframe['normalized_count'].sum()
        percentage = (below_threshold / total) * 100

        plt.text(ml_threshold, plt.gca().get_ylim()[1], 
                 f'{percentage:.2f}% below threshold', 
                 rotation=90, va='top', ha='right', color='red')

    # Improved x-axis ticks
    plt.xticks(np.arange(0, 257, 32))
    plt.xlim(0, 256)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Save high-resolution figure
    output_name = f'{prefix}_ML-{mod}.png'
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    plt.close()  # Use close() instead of clf() to free up memory

def main():
    """
    Main function to parse arguments and process BAM file
    """
    parser = argparse.ArgumentParser(description="Analyze modifications in BAM files")
    
    parser.add_argument('input_bam', type=str, help='Input BAM file path')
    parser.add_argument('output_prefix', type=str, help='Output PNG filename prefix')
    parser.add_argument('mod_code', type=str, help='Modification code to analyze (e.g., "a" for A+a, "m" for C+m)')
    parser.add_argument('-s', '--sample_size', type=int, default=None, 
                        help='Number of reads to sample (default: None, use all reads)')
    parser.add_argument('-n', '--normalize_by', type=int, default=1000, 
                        help='Base pair length for normalization (default: 1000 for kb)')
    parser.add_argument('-t', '--threshold', type=float, default=None, 
                        help='Threshold value between 0 and 1 for vertical line (default: None)')
    parser.add_argument('--mapq', type=int, default=10, 
                        help='Minimum mapping quality for read filtering (default: 10)')

    args = parser.parse_args()

    # Validate inputs
    if args.mod_code not in ['a', 'm']:
        raise ValueError("Modification code must be 'a' or 'm'")

    if args.threshold is not None and (args.threshold < 0 or args.threshold > 1):
        raise ValueError("Threshold must be between 0 and 1")

    # Sample reads 
    samples = sample_reads(args.input_bam, num_reads=args.sample_size, min_mapq=args.mapq)

    # Process and plot modifications for specific mod_code
    df, mod, length = create_dataframes(samples, args.mod_code)
    plot_histogram(df, mod, length, args.output_prefix, 
                   args.normalize_by, args.threshold)

if __name__ == "__main__":
    main()