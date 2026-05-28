class FASTAReader:
    # Reader file FASTA.
    def __init__(self, filepath: str) -> None:
        self.filepath: str = filepath
        self.header: str | None = None
        self.sequence: str | None = None
        self.sequence_length: int = 0

    def read_fasta(self) -> tuple[str, str]:
        # Baca header dan sekuens.
        try:
            with open(self.filepath, 'r') as file:
                lines = file.readlines()
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"File tidak ditemukan: {self.filepath}") from exc

        if not lines:
            raise ValueError("File FASTA kosong")

        if not lines[0].startswith('>'):
            raise ValueError("File FASTA tidak valid")

        header = lines[0].strip()
        sequence = ''.join(line.strip() for line in lines[1:] if line.strip())

        self.header = header
        self.sequence = sequence
        self.sequence_length = len(sequence)

        return header, sequence

    def get_subsequence(self, start: int, end: int) -> str:
        # Ambil subsequence 1-based.
        if self.sequence is None:
            raise ValueError("Sequence belum dibaca")

        start_idx = start - 1
        end_idx = end

        if start_idx < 0 or end_idx > len(self.sequence) or start > end:
            raise ValueError(f"Posisi diluar range: {start}-{end}")

        return self.sequence[start_idx:end_idx]
