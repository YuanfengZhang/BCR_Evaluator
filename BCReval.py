# --coding: utf-8 --
"""
This version of BCReval.py is now parallelized and accelerated using Cpp module.
To use this verision of BCReval.py, please:
1. create a conda env with python>=3.13, scipy, pybind11 and zstandard (if you need it), e.g.:
a) mamba create -n pybind -c conda-forge python=3.13 scipy pybind11 zstandard
or
b) uv python install 3.13 && cd BCR_Evaluator && \
   uv python pin 3.13 && uv init && uv add pybind11

2. g++ -O3 -shared -std=c++11 -fPIC $(python3-config --includes) \
     telomere_parser.cpp \
     -o telomere_parser$(python3-config --extension-suffix)
3. run the BCREval.py. See the instructions in help message.

Modified by: Zhang Yuanfeng zhangyuanfeng1997@foxmail.com
Last update: 2025-05-26
Version: 1.4
**********************************************************************************
Original statement from https://github.com/hqyone/BCR_Evaluator:

BCR evaluator is used to evaluate bisulfite conversion ratio in WGSBS experiment
Author : Quanyuan He Ph.D
Email : hqyone@hunnu.edu.com
Insititution: School of Medicine, Hunan Normal University
Licensed under the MIT License
Last update: 2019/02/11
version: 1.3
**********************************************************************************
"""
from argparse import ArgumentParser, Namespace
from collections.abc import Iterator, Callable
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures._base import Future
from dataclasses import dataclass
from functools import lru_cache
import gzip
from itertools import islice
import logging
from pathlib import Path
from typing import Literal, TextIO
from scipy.optimize import minimize, OptimizeResult
import zstandard
from telomere_parser import telomere_strand_seq


@dataclass
class StrandConfig:
    patterns: set[str]
    n3_pattern: str
    full_methylated_unit: str
    default_seq_counts: dict[str, int]
    site1_patterns: set[str]
    site2_patterns: set[str]
    site3_patterns: set[str]


STRAND_CONFIGS: dict[str, StrandConfig] = {
    'C': StrandConfig(patterns={'TTTTAA', 'CTTTAA', 'TCTTAA', 'TTCTAA',
                                'CCTTAA', 'TCCTAA', 'CTCTAA', 'CCCTAA'},
                      n3_pattern='CCCTAA', full_methylated_unit='TTTTAA',
                      default_seq_counts={'TTTTAA': 0, 'CTTTAA': 0, 'TCTTAA': 0, 'TTCTAA': 0,
                                          'CCTTAA': 0, 'TCCTAA': 0, 'CTCTAA': 0, 'CCCTAA': 0},
                      site1_patterns={'CTTTAA', 'CCTTAA', 'CTCTAA', 'CCCTAA'},
                      site2_patterns={'TCTTAA', 'CCTTAA', 'TCCTAA', 'CCCTAA'},
                      site3_patterns={'TTCTAA', 'TCCTAA', 'CTCTAA', 'CCCTAA'}),
    'G': StrandConfig(patterns={'TTAGGG', 'TTAGGA', 'TTAGAG', 'TTAAGG',
                                'TTAGAA', 'TTAAAG', 'TTAAGA', 'TTAAAA'},
                      n3_pattern='TTAGGG', full_methylated_unit='TTAAAA',
                      default_seq_counts={'TTAGGG': 0, 'TTAGGA': 0, 'TTAGAG': 0, 'TTAAGG': 0,
                                          'TTAGAA': 0, 'TTAAAG': 0, 'TTAAGA': 0, 'TTAAAA': 0},
                      site1_patterns={'TTAAAG', 'TTAAGG', 'TTAGAG', 'TTAGGG'},
                      site2_patterns={'TTAAGA', 'TTAAGG', 'TTAGGA', 'TTAGGG'},
                      site3_patterns={'TTAGAA', 'TTAGGA', 'TTAGAG', 'TTAGGG'})}


def solve_ucr(r1: float, r2: float, r3: float) -> float:
    def objective(ucr: list[float]):
        _ucr = ucr[0]
        return (3 * _ucr * (1 - _ucr) ** 2 - r1) ** 2 +\
            (3 * _ucr ** 2 * (1 - _ucr) - r2) ** 2 +\
            (_ucr ** 3 - r3) ** 2

    if r1 + r2 != 0:
        initial_guess: list[float] = [r2 / (r1 + r2)]
    else:
        initial_guess = [r3 ** (1 / 3)]

    result: OptimizeResult = minimize(objective,
                                      initial_guess,
                                      bounds=[(0.0, 1.0)],
                                      method='L-BFGS-B')

    if not result.success:
        logging.error(f'optimization failed: {result.message}')
        return r2 / (r1 + r2) if r1 + r2 != 0 else r3 ** (1 / 3)
    return result.x[0]


class ChunkStatistics:
    def __init__(self, strand: Literal['C', 'G'] = 'C', total_reads: int = 0,
                 telomere_reads: int = 0, strand_n3: int = 0, unit_length: int = 6,
                 strand_lengths: list[int] | None = None, n3_pattern: str | None = None,
                 strand_seq_counts: dict[str, int] | None = None):

        self.strand: Literal['C'] | Literal['G'] = strand
        self.total_reads: int = total_reads
        self.telomere_reads: int = telomere_reads
        self.strand_n3: int = strand_n3
        self.unit_length: int = unit_length

        self.strand_lengths = strand_lengths if strand_lengths else [0] * 30

        if strand not in STRAND_CONFIGS:
            raise ValueError(f'Invalid strand type: {strand}. Must be "C" or "G".')

        config: StrandConfig = STRAND_CONFIGS[strand]
        self.strand_patterns: set[str] = config.patterns
        self.n3_pattern: str = n3_pattern if n3_pattern else config.n3_pattern
        self.full_methylated_unit: str = config.full_methylated_unit
        self.site1_patterns: set[str] = config.site1_patterns
        self.site2_patterns: set[str] = config.site2_patterns
        self.site3_patterns: set[str] = config.site3_patterns

        self.strand_seq_counts: dict[str, int]
        if strand_seq_counts:
            self.strand_seq_counts = strand_seq_counts.copy()
        else:
            self.strand_seq_counts = config.default_seq_counts.copy()

    def __add__(self, other: 'ChunkStatistics') -> 'ChunkStatistics':
        if self.strand == other.strand:
            return ChunkStatistics(
                strand=self.strand,
                total_reads=self.total_reads + other.total_reads,
                telomere_reads=self.telomere_reads + other.telomere_reads,
                strand_n3=self.strand_n3 + other.strand_n3,
                strand_lengths=[x + y for x, y in zip(self.strand_lengths, other.strand_lengths)],
                strand_seq_counts={
                    k: self.strand_seq_counts.get(k, 0) + other.strand_seq_counts.get(k, 0)
                    for k in self.strand_patterns})
        else:
            raise ValueError('Cannot add ChunkStatistics with different strand types.')

    def conversion_rate(self) -> dict[str, str | int | float]:
        unit_counts: int = sum(self.strand_seq_counts.values())
        completely_converted_unit_counts: int = self.strand_seq_counts[self.n3_pattern]

        r1: float
        r2: float
        r3: float
        ucr_r1r2: float
        ucr_r3: float
        optimized_ucr: float
        if unit_counts == 0:
            r1 = r2 = r3 = ucr_r1r2 = ucr_r3 = optimized_ucr = -1.0
            logging.warning(f'No valid units found for strand {self.strand}. '
                            f'Conversion rates will be set to -1.0.')
        else:
            site_counts: list[int] = [
                sum(self.strand_seq_counts[pattern] for pattern in patterns)
                for patterns in (self.site1_patterns, self.site2_patterns, self.site3_patterns)]

            (r1,
             r2,
             r3) = [_c / unit_counts if unit_counts > 0 else -1.0
                    for _c in site_counts]

            ucr_r1r2 = r2 / (r1 + r2)
            ucr_r3 = r3 ** (1 / 3)

            optimized_ucr = solve_ucr(r1=r1, r2=r2, r3=r3)

        return {
            'strand': self.strand,
            'total_reads': self.total_reads,
            'telomere_reads': self.telomere_reads,
            'completely_converted_unit_counts': completely_converted_unit_counts,
            'n3_reads': self.strand_n3,
            'strand_seq_counts': ', '.join(f'{k}:{v}'
                                           for k, v in self.strand_seq_counts.items()),
            'r1': r1, 'r2': r2, 'r3': r3,
            'optimized_conversion_rate': 1 - optimized_ucr,
            'estimated_conversion_rate': 1 - ((r1 + r2 + r3) / 3),
            'conversion_rate_r1r2': 1 - ucr_r1r2,
            'conversion_rate_r3': 1 - ucr_r3}


def read_fastq(f_path: Path, chunk_size: int) -> Iterator[list[str]]:
    """
    Reads a FASTQ file (optionally gzipped or zstandard-compressed) in chunks.

    Args:
        f_path (Path): Path to the FASTQ file.
        chunk_size (int): Number of reads per chunk.

    Yields:
        list[str]: A list of lines representing chunk_size reads (4 lines per read).
    """
    read_methods: dict[str, Callable[[Path], TextIO]] = {
        '.gz': lambda f: gzip.open(f, 'rt'),
        'zst': lambda f: zstandard.open(f, 'rt')  # type: ignore
    }

    with read_methods.get(f_path.suffix[: -3],
                          lambda f: open(f, 'r'))(f_path) as f:
        while True:
            chunk = list(islice(f, chunk_size * 4))
            if not chunk:
                break
            yield chunk


@lru_cache(maxsize=1024)
def best_match(seq: str) -> tuple[str, str]:
    return telomere_strand_seq(seq)


def process_chunk(chunk: list[str]) -> dict[str, ChunkStatistics]:
    chunk_stats: dict[str, ChunkStatistics] = {'C': ChunkStatistics(strand='C'),
                                               'G': ChunkStatistics(strand='G')}

    for i in range(0, len(chunk), 4):
        try:
            seq: str = chunk[i + 1].strip()
        except IndexError:
            break

        telomere_strand, match_seq = best_match(seq)
        units = len(match_seq) // 6

        stats: ChunkStatistics | None = chunk_stats.get(telomere_strand, None)
        if not stats:
            logging.warning(f'Unknown telomere strand for sequence: {seq}, skipping.')
            break

        stats.strand_lengths[units] += 1

        if units >= 6:
            stats.telomere_reads += 1

            if stats.n3_pattern in match_seq:
                stats.strand_n3 += 1

            for pos in range(0, len(match_seq), stats.unit_length):
                single_unit = match_seq[pos: pos + stats.unit_length]
                if single_unit in stats.strand_seq_counts:
                    stats.strand_seq_counts[single_unit] += 1

    return chunk_stats


def main():
    arg_parser: ArgumentParser = ArgumentParser(
        description=('Evaluate base conversion ratio in NGS methylation fastq file,\n'
                     'support plain text / gzip / zstd.\n'
                     'make sure the name of your gzipped fastq file ends with .gz, '
                     'and zstd compressed file ends with .zst.'),
        epilog='example: python -m BCReval.py -i A.fq.gz -o r1.bcr -t 4 -v info.')
    arg_parser.add_argument('-i', '--input', type=str, dest='input', required=True,
                            help='The path of single input fastq file')
    arg_parser.add_argument('-o', '--output', type=str, dest='output', required=True,
                            help='The path of output file, the parent dir will be created if not exists')
    arg_parser.add_argument('-t', '--threads', type=int, dest='threads', default=1,
                            help='The number of threads to use, default is 1')
    arg_parser.add_argument('-c', '--chunk-size', dest='chunk_size', type=int, default=1000000,
                            help='The number of reads to process in each chunk, default is 1000000')
    arg_parser.add_argument('-v', '--verbose', dest='log_level',
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')

    args: Namespace = arg_parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s - %(levelname)s - %(message)s')

    logging.debug('checking args...')
    input_f: Path = Path(args.input)
    if not input_f.exists():
        raise FileNotFoundError(f'Input file {input_f} does not exist.')

    output_f: Path = Path(args.output)
    if not output_f.parent.exists():
        logging.info(f'Creating parent directory for output file: {output_f.parent}')
        output_f.parent.mkdir(parents=True, exist_ok=True)

    threads: int = args.threads if args.threads > 0 else 1

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures: list[Future[dict[str,
                                  ChunkStatistics]]] = [executor.submit(process_chunk, chunk)
                                                        for chunk in read_fastq(f_path=input_f,
                                                                                chunk_size=args.chunk_size)]
        results: list[dict[str, ChunkStatistics]] = [future.result() for future in futures]

    agg_c_stats: dict[str,
                      str | int | float] = sum((result['C'] for result in results),
                                               start=ChunkStatistics(strand='C')).conversion_rate()

    agg_g_stats: dict[str,
                      str | int | float] = sum((result['G'] for result in results),
                                               start=ChunkStatistics(strand='G')).conversion_rate()

    with open(output_f, 'w+') as out_io:
        out_io.write('input_fname\tstrand\ttotal_reads\ttelomere_reads\tcompletely_converted_unit_counts\t'
                     'n3_reads\tstrand_seq_counts\tr1\tr2\tr3\toptimized_conversion_rate\t'
                     'estimated_conversion_rate\tconversion_rate_r1r2\tconversion_rate_r3\n')
        for agg_stats in (agg_c_stats, agg_g_stats):
            out_io.write('\t'.join([input_f.name] + [
                str(agg_stats[key])
                for key in ['strand', 'total_reads', 'telomere_reads',
                            'completely_converted_unit_counts', 'n3_reads',
                            'strand_seq_counts', 'r1', 'r2', 'r3',
                            'optimized_conversion_rate',
                            'estimated_conversion_rate',
                            'conversion_rate_r1r2',
                            'conversion_rate_r3']]) + '\n')


if __name__ == '__main__':
    main()
