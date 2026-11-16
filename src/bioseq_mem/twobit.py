"""Minimal reader for UCSC .2bit genomes, keeping soft-masking (repeats in lowercase).

Format: https://genome.ucsc.edu/FAQ/FAQformat.html#format7
"""

import mmap
import struct

import numpy as np

_BASES = np.frombuffer(b"TCAG", dtype=np.uint8)


class TwoBit:
    def __init__(self, path: str):
        self._f = open(path, "rb")
        self._m = mmap.mmap(self._f.fileno(), 0, access=mmap.ACCESS_READ)
        sig, version, count, _ = struct.unpack_from("<IIII", self._m, 0)
        if sig != 0x1A412743 or version != 0:
            raise ValueError(f"{path}: not a version-0 little-endian 2bit file")
        self.offsets, pos = {}, 16
        for _ in range(count):
            n = self._m[pos]
            name = self._m[pos + 1 : pos + 1 + n].decode()
            (self.offsets[name],) = struct.unpack_from("<I", self._m, pos + 1 + n)
            pos += 1 + n + 4
        self._meta = {}

    def _record(self, chrom: str):
        if chrom not in self._meta:
            pos = self.offsets[chrom]
            size, n_count = struct.unpack_from("<II", self._m, pos)
            pos += 8
            n_starts = np.frombuffer(self._m, "<u4", n_count, pos)
            n_sizes = np.frombuffer(self._m, "<u4", n_count, pos + 4 * n_count)
            pos += 8 * n_count
            (m_count,) = struct.unpack_from("<I", self._m, pos)
            pos += 4
            m_starts = np.frombuffer(self._m, "<u4", m_count, pos)
            m_sizes = np.frombuffer(self._m, "<u4", m_count, pos + 4 * m_count)
            pos += 8 * m_count + 4  # + reserved
            self._meta[chrom] = (size, n_starts, n_sizes, m_starts, m_sizes, pos)
        return self._meta[chrom]

    def size(self, chrom: str) -> int:
        return self._record(chrom)[0]

    def get(self, chrom: str, start: int, end: int) -> str:
        """Bases [start, end), 0-based; repeats lowercase, gaps N."""
        size, n_starts, n_sizes, m_starts, m_sizes, dna = self._record(chrom)
        end = min(end, size)
        b0, b1 = start // 4, (end + 3) // 4
        packed = np.frombuffer(self._m, np.uint8, b1 - b0, dna + b0)
        codes = np.stack([(packed >> s) & 3 for s in (6, 4, 2, 0)], axis=1).ravel()
        seq = _BASES[codes[start - b0 * 4 : end - b0 * 4]].copy()
        for starts, sizes, fix in ((m_starts, m_sizes, "lower"), (n_starts, n_sizes, "N")):
            i = np.searchsorted(starts, end)  # blocks starting before end
            for s, n in zip(starts[:i][-64:], sizes[:i][-64:]):  # blocks are sorted and short
                a, b = max(s, start), min(s + n, end)
                if a < b:
                    sl = slice(a - start, b - start)
                    seq[sl] = (seq[sl] + 32) if fix == "lower" else ord("N")
        return seq.tobytes().decode()
