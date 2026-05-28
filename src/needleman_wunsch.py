class NeedlemanWunschDP:
    # Implementasi global alignment DP.
    def __init__(
        self,
        seq1: str,
        seq2: str,
        match_score: int = 2,
        mismatch_score: int = -1,
        gap_penalty: int = -1,
    ) -> None:
        self.seq1: str = seq1.upper()
        self.seq2: str = seq2.upper()
        self.match_score: int = match_score
        self.mismatch_score: int = mismatch_score
        self.gap_penalty: int = gap_penalty

        self.m: int = len(seq1) + 1
        self.n: int = len(seq2) + 1

        self.dp_matrix: list[list[int]] = self._initialize_dp_matrix()
        self.traceback_matrix: list[list[str]] | None = None

    def _initialize_dp_matrix(self) -> list[list[int]]:
        # Isi baris kolom awal.
        dp_matrix: list[list[int]] = [[0 for _ in range(self.n)] for _ in range(self.m)]

        for j in range(self.n):
            dp_matrix[0][j] = j * self.gap_penalty

        for i in range(self.m):
            dp_matrix[i][0] = i * self.gap_penalty

        return dp_matrix

    def display_dp_matrix(self, max_rows: int = 12, max_cols: int = 12) -> None:
        # Tampilkan potongan matriks DP.
        print(f"\nDP Matrix ({min(max_rows, self.m)} x {min(max_cols, self.n)} dari {self.m} x {self.n}):\n")

        header = "      "
        for j in range(min(max_cols, self.n)):
            if j == 0:
                header += "  -"
            else:
                header += f"  {self.seq2[j-1]}"
        print(header)

        for i in range(min(max_rows, self.m)):
            if i == 0:
                row_label = "  -"
            else:
                row_label = f"  {self.seq1[i-1]}"

            row_str = f"{row_label:>4}"
            for j in range(min(max_cols, self.n)):
                row_str += f"{self.dp_matrix[i][j]:>4}"

            print(row_str)

        if max_rows < self.m or max_cols < self.n:
            print(f"\n  [Sisa baris/kolom tidak ditampilkan]")

    def get_matrix_stats(self) -> dict[str, int]:
        # Statistik dimensi matriks.
        flat_matrix: list[int] = [val for row in self.dp_matrix for val in row]

        stats: dict[str, int] = {
            'rows': self.m,
            'cols': self.n,
            'total_cells': self.m * self.n,
            'min_value': min(flat_matrix),
            'max_value': max(flat_matrix),
            'seq1_length': len(self.seq1),
            'seq2_length': len(self.seq2),
        }

        return stats

    def export_to_json(
        self,
        filepath: str,
        provenance: dict[str, dict[str, str]] | None = None,
    ) -> None:
        # Ekspor stage 1 JSON.
        import os
        import json
        from datetime import datetime

        os.makedirs(os.path.dirname(filepath), exist_ok=True)

        data = {
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'stage': 'STAGE_1_INITIALIZATION',
                'seq1_length': len(self.seq1),
                'seq2_length': len(self.seq2),
                'matrix_dimension': f"{self.m}x{self.n}",
                'dataset_provenance': provenance or {},
            },
            'parameters': {
                'match_score': self.match_score,
                'mismatch_score': self.mismatch_score,
                'gap_penalty': self.gap_penalty,
            },
            'sequences': {
                'seq1_preview': self.seq1[:100] + ('...' if len(self.seq1) > 100 else ''),
                'seq2_preview': self.seq2[:100] + ('...' if len(self.seq2) > 100 else ''),
            },
            'statistics': self.get_matrix_stats(),
        }

        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)

        print(f"Export: {filepath}")

    def fill_matrix(self) -> int:
        # Isi matriks DP.
        self.traceback_matrix = [['' for _ in range(self.n)] for _ in range(self.m)]

        # Arah awal baris kolom.
        for j in range(1, self.n):
            self.traceback_matrix[0][j] = 'L'
        for i in range(1, self.m):
            self.traceback_matrix[i][0] = 'U'

        for i in range(1, self.m):
            for j in range(1, self.n):
                score = self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_score

                diagonal = self.dp_matrix[i - 1][j - 1] + score
                up       = self.dp_matrix[i - 1][j]     + self.gap_penalty
                left     = self.dp_matrix[i][j - 1]     + self.gap_penalty

                best = max(diagonal, up, left)
                self.dp_matrix[i][j] = best

                # Tie-breaker prioritas diagonal.
                if best == diagonal:
                    self.traceback_matrix[i][j] = 'D'
                elif best == up:
                    self.traceback_matrix[i][j] = 'U'
                else:
                    self.traceback_matrix[i][j] = 'L'

        return self.dp_matrix[self.m - 1][self.n - 1]

    def traceback(self) -> tuple[str, str, str]:
        # Traceback dari kanan bawah.
        if self.traceback_matrix is None:
            raise ValueError("Traceback matrix kosong. Jalankan fill_matrix() terlebih dahulu.")

        aligned_seq1: list[str] = []
        aligned_seq2: list[str] = []
        alignment_bar: list[str] = []

        i, j = self.m - 1, self.n - 1

        while i > 0 or j > 0:
            direction = self.traceback_matrix[i][j]

            if direction == 'D':
                aligned_seq1.append(self.seq1[i - 1])
                aligned_seq2.append(self.seq2[j - 1])
                alignment_bar.append('|' if self.seq1[i - 1] == self.seq2[j - 1] else '*')
                i -= 1
                j -= 1
            elif direction == 'U':
                aligned_seq1.append(self.seq1[i - 1])
                aligned_seq2.append('-')
                alignment_bar.append(' ')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(self.seq2[j - 1])
                alignment_bar.append(' ')
                j -= 1

        seq1_str = ''.join(reversed(aligned_seq1))
        seq2_str = ''.join(reversed(aligned_seq2))
        bar_str = ''.join(reversed(alignment_bar))

        return seq1_str, seq2_str, bar_str

    def get_alignment_stats(
        self,
        aligned_seq1: str,
        aligned_seq2: str,
        alignment_bar: str,
    ) -> dict[str, int | float]:
        # Hitung statistik alignment.
        matches    = alignment_bar.count('|')
        mismatches = alignment_bar.count('*')
        gaps       = alignment_bar.count(' ')
        length     = len(alignment_bar)
        identity   = (matches / length * 100) if length > 0 else 0.0

        return {
            'alignment_length': length,
            'matches':    matches,
            'mismatches': mismatches,
            'gaps':       gaps,
            'identity_pct': round(identity, 4),
            'alignment_score': self.dp_matrix[self.m - 1][self.n - 1],
        }

    def display_alignment(
        self,
        aligned_seq1: str,
        aligned_seq2: str,
        alignment_bar: str,
        line_width: int = 60,
    ) -> None:
        # Cetak alignment ke terminal.
        print(f"\nHasil Alignment (total {len(alignment_bar)} karakter):\n")
        for start in range(0, len(alignment_bar), line_width):
            end = start + line_width
            print(f"  Seq1 : {aligned_seq1[start:end]}")
            print(f"         {alignment_bar[start:end]}")
            print(f"  Seq2 : {aligned_seq2[start:end]}")
            print()

    def export_alignment_to_json(
        self,
        filepath: str,
        aligned_seq1: str,
        aligned_seq2: str,
        alignment_bar: str,
    ) -> None:
        # Ekspor alignment ke JSON.
        import os
        import json
        from datetime import datetime

        os.makedirs(os.path.dirname(filepath), exist_ok=True)

        stats = self.get_alignment_stats(aligned_seq1, aligned_seq2, alignment_bar)

        data = {
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'stage': 'STAGE_2_ALIGNMENT',
                'seq1_length': len(self.seq1),
                'seq2_length': len(self.seq2),
                'matrix_dimension': f"{self.m}x{self.n}",
            },
            'parameters': {
                'match_score':    self.match_score,
                'mismatch_score': self.mismatch_score,
                'gap_penalty':    self.gap_penalty,
            },
            'alignment_statistics': stats,
            'alignment_preview': {
                'seq1': aligned_seq1[:200] + ('...' if len(aligned_seq1) > 200 else ''),
                'bar':  alignment_bar[:200] + ('...' if len(alignment_bar) > 200 else ''),
                'seq2': aligned_seq2[:200] + ('...' if len(aligned_seq2) > 200 else ''),
            },
        }

        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)

        print(f"Export: {filepath}")
