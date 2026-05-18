import os
from utils import FASTAReader
from needleman_wunsch import NeedlemanWunschDP

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_PENALTY = -1

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATASET_DIR = os.path.join(BASE_DIR, 'dataset')
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')

SEQUENCE_1973_PATH = os.path.join(DATASET_DIR, 'sequence_1973.fasta')
SEQUENCE_2020_PATH = os.path.join(DATASET_DIR, 'sequence_2020.fasta')

ENVELOPE_START = 499
ENVELOPE_END = 1977
OUTPUT_JSON = os.path.join(OUTPUT_DIR, 'stage_1_initialization.json')


def main():
    print("DETEKSI MUTASI VIRUS ENDEMIK - ORANG 1".center(80))
    print("Needleman-Wunsch Initialization".center(80))
    
    print("\n[TAHAP 1] Membaca File FASTA")
    print("-" * 80)
    
    reader_1973 = FASTAReader(SEQUENCE_1973_PATH)
    header_1973, sequence_1973 = reader_1973.read_fasta()
    print(f"Sequence 1973: {len(sequence_1973)} bp")
    
    reader_2020 = FASTAReader(SEQUENCE_2020_PATH)
    header_2020, sequence_2020 = reader_2020.read_fasta()
    print(f"Sequence 2020: {len(sequence_2020)} bp")
    
    print("\n[TAHAP 2] Slicing Gen Envelope")
    print("-" * 80)
    print(f"Source: Sequence 1973, Basis {ENVELOPE_START}-{ENVELOPE_END}")
    
    envelope_1973 = reader_1973.get_subsequence(ENVELOPE_START, ENVELOPE_END)
    envelope_2020 = sequence_2020
    
    print(f"Envelop 1973: {len(envelope_1973)} bp")
    print(f"Envelop 2020: {len(envelope_2020)} bp")
    
    print("\n[TAHAP 3] Inisialisasi Matriks DP")
    print("-" * 80)
    
    nw = NeedlemanWunschDP(
        seq1=envelope_1973,
        seq2=envelope_2020,
        match_score=MATCH_SCORE,
        mismatch_score=MISMATCH_SCORE,
        gap_penalty=GAP_PENALTY
    )
    
    print("\n[TAHAP 4] Hasil Inisialisasi")
    print("-" * 80)
    
    stats = nw.get_matrix_stats()
    print(f"\nMatrix Stats:")
    print(f"  Dimensi: {stats['rows']} x {stats['cols']}")
    print(f"  Total cells: {stats['total_cells']:,}")
    print(f"  Min/Max: {stats['min_value']} / {stats['max_value']}")
    
    nw.display_dp_matrix(max_rows=12, max_cols=12)
    
    nw.export_to_json(OUTPUT_JSON)
    
    return nw, envelope_1973, envelope_2020


if __name__ == "__main__":
    nw, seq1, seq2 = main()