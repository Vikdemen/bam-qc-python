import logging
import time
from argparse import ArgumentParser
from concurrent.futures import as_completed, ProcessPoolExecutor

from pysam import AlignedSegment, AlignmentFile


def main():
    args = parse_arguments()
    bam_file = args.alignments
    logging.basicConfig(level=logging.INFO)
    start = time.perf_counter()
    gc = calculate_mean_gc(bam_file)
    end = time.perf_counter()
    total = end - start
    logging.info(f"Total: {total // 60} min {total % 60} seconds")
    print(gc)


def calculate_mean_gc(filepath: str) -> float:
    """
    :param filepath: Path to a BAM file
    :return: Mean GC content of uniquely mapped reads
    """
    bam = AlignmentFile(filepath, "rb")
    try:
        references = bam.references
    finally:
        bam.close()

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_region, filepath, refseq) for refseq in references]
        gc_sum = 0
        total_unique = 0
        for future in as_completed(futures):
            region_gc, region_count = future.result()
            gc_sum += region_gc
            total_unique += region_count
        average_gc = gc_sum / total_unique

    return average_gc


def process_region(filepath, region: str) -> tuple[float, int]:
    """
    :param filepath: A path to indexed BAM file
    :param region: Region that would be processed
    :return: Sum of GC content of unique aligned reads and total number of said reads.
    """
    gc_sum = 0
    total_unique = 0
    bam = AlignmentFile(filepath, "rb")
    try:
        for read in bam.fetch(region):
            # skipping non-unique reads
            if not is_valid(read):
                continue
            gc_sum += calculate_gc(read)
            total_unique += 1
    finally:
        bam.close()
    return gc_sum, total_unique


def is_valid(read: AlignedSegment) -> bool:
    """
    :param read: Alignment record
    :return: If it's a first element of paired-end read, properly mapped, and is primary alignment with non-zero
    mapping quality
    """
    if read.is_unmapped:
        return False
    if not read.is_proper_pair:
        return False
    if read.is_secondary:
        return False
    if read.mapping_quality == 0:
        return False
    if read.is_read2:
        return False
    return True


def calculate_gc(read: AlignedSegment) -> float:
    # for now, doesn't account for matching bases from CIGAR
    seq = read.seq
    # should be a little bit faster than sets
    gc = sum(char == 'G' or char == 'C' for char in seq) / len(seq)
    return gc


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('alignments', type=str)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()

