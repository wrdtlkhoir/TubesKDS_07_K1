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
OUTPUT_JSON            = os.path.join(OUTPUT_DIR, 'stage_1_initialization.json')
OUTPUT_ALIGNMENT_JSON  = os.path.join(OUTPUT_DIR, 'stage_2_alignment.json')
OUTPUT_ANALYSIS_JSON   = os.path.join(OUTPUT_DIR, 'stage_3_analysis.json')
OUTPUT_MUTATIONS_CSV   = os.path.join(OUTPUT_DIR, 'stage_3_mutation_positions.csv')
OUTPUT_REPORT_TXT      = os.path.join(OUTPUT_DIR, 'stage_3_report_notes.txt')


def build_mutation_events(aligned_seq1, aligned_seq2):
    events = []
    pos1 = 0
    pos2 = 0

    for alignment_index, (base1, base2) in enumerate(zip(aligned_seq1, aligned_seq2), start=1):
        if base1 != '-':
            pos1 += 1
        if base2 != '-':
            pos2 += 1

        if base1 == base2:
            continue

        if base1 == '-':
            event_type = 'insertion_in_2020'
            ref_position = None
            query_position = pos2
        elif base2 == '-':
            event_type = 'deletion_in_2020'
            ref_position = pos1
            query_position = None
        else:
            event_type = 'substitution'
            ref_position = pos1
            query_position = pos2

        events.append({
            'alignment_index': alignment_index,
            'event_type': event_type,
            'position_1973_envelope': ref_position,
            'position_2020': query_position,
            'base_1973': base1,
            'base_2020': base2,
        })

    return events


def format_alignment_blocks(aligned_seq1, aligned_seq2, alignment_bar, line_width=70, max_blocks=4):
    blocks = []
    for start in range(0, len(alignment_bar), line_width):
        if len(blocks) >= max_blocks:
            break
        end = start + line_width
        blocks.append(
            f"Pos {start + 1}-{min(end, len(alignment_bar))}\n"
            f"1973 : {aligned_seq1[start:end]}\n"
            f"       {alignment_bar[start:end]}\n"
            f"2020 : {aligned_seq2[start:end]}"
        )
    return "\n\n".join(blocks)


def export_stage_3_analysis(stats, mutation_events, aligned_seq1, aligned_seq2, alignment_bar):
    import csv
    import json
    from datetime import datetime

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    substitutions = sum(1 for event in mutation_events if event['event_type'] == 'substitution')
    insertions = sum(1 for event in mutation_events if event['event_type'] == 'insertion_in_2020')
    deletions = sum(1 for event in mutation_events if event['event_type'] == 'deletion_in_2020')

    analysis = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'stage': 'STAGE_3_ANALYSIS',
            'comparison': 'Dengue 1973 envelope gene vs Dengue 2020 envelope protein E region',
        },
        'quantitative_summary': {
            **stats,
            'mutation_events': len(mutation_events),
            'substitutions': substitutions,
            'insertions_in_2020': insertions,
            'deletions_in_2020': deletions,
        },
        'mutation_positions': mutation_events,
        'full_alignment': {
            'seq1_1973': aligned_seq1,
            'bar': alignment_bar,
            'seq2_2020': aligned_seq2,
        },
    }

    with open(OUTPUT_ANALYSIS_JSON, 'w', encoding='utf-8') as file:
        json.dump(analysis, file, indent=2)

    with open(OUTPUT_MUTATIONS_CSV, 'w', newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=[
            'alignment_index',
            'event_type',
            'position_1973_envelope',
            'position_2020',
            'base_1973',
            'base_2020',
        ])
        writer.writeheader()
        writer.writerows(mutation_events)

    report_notes = "\n".join([
        "# Ringkasan Hasil Alignment Needleman-Wunsch",
        "",
        "# Dataset",
        f"Referensi 1973: gen Envelope Dengue hasil slicing posisi 499-1977 ({stats['alignment_length']} bp)",
        f"Varian 2020: region Envelope protein E Dengue ({stats['alignment_length']} bp)",
        "",
        "# Parameter Scoring",
        "Match = +2",
        "Mismatch = -1",
        "Gap = -1",
        "",
        "# Hasil Kuantitatif",
        f"Skor alignment akhir = {stats['alignment_score']}",
        f"Panjang alignment = {stats['alignment_length']} bp",
        f"Jumlah match = {stats['matches']}",
        f"Jumlah mismatch/substitusi = {substitutions}",
        f"Jumlah gap = {stats['gaps']} (insersi {insertions}, delesi {deletions})",
        f"Persentase identity = {stats['identity_pct']}%",
        f"Total posisi berbeda = {len(mutation_events)}",
        "",
        "# Interpretasi Singkat",
        f"Hasil alignment menunjukkan identity sebesar {stats['identity_pct']}%, sehingga sekuens Envelope Dengue tahun 1973 dan 2020 masih memiliki tingkat kemiripan yang tinggi.",
        f"Walaupun demikian, terdapat {len(mutation_events)} posisi berbeda pada region yang dibandingkan. Karena jumlah gap adalah {stats['gaps']}, perbedaan yang muncul didominasi oleh substitusi basa, bukan insersi atau delesi.",
        "",
        "# Cuplikan Alignment",
        format_alignment_blocks(aligned_seq1, aligned_seq2, alignment_bar),
        "",
    ])

    with open(OUTPUT_REPORT_TXT, 'w', encoding='utf-8') as file:
        file.write(report_notes)

    print(f"File JSON analisis: {OUTPUT_ANALYSIS_JSON}")
    print(f"File CSV posisi mutasi: {OUTPUT_MUTATIONS_CSV}")
    print(f"Ringkasan hasil untuk laporan: {OUTPUT_REPORT_TXT}")


def main():
    print("DETEKSI MUTASI DENGUE DENGAN NEEDLEMAN-WUNSCH".center(80))
    print("Analisis Sekuens DNA Varian Dengue".center(80))
    
    print("\n[1] Pembacaan File FASTA")
    print("-" * 80)
    
    reader_1973 = FASTAReader(SEQUENCE_1973_PATH)
    header_1973, sequence_1973 = reader_1973.read_fasta()
    print(f"Sekuens 1973: {len(sequence_1973)} bp")
    
    reader_2020 = FASTAReader(SEQUENCE_2020_PATH)
    header_2020, sequence_2020 = reader_2020.read_fasta()
    print(f"Sekuens 2020: {len(sequence_2020)} bp")
    
    print("\n[2] Slicing Gen Envelope")
    print("-" * 80)
    print(f"Region acuan dari sekuens 1973: posisi {ENVELOPE_START}-{ENVELOPE_END}")
    
    envelope_1973 = reader_1973.get_subsequence(ENVELOPE_START, ENVELOPE_END)
    envelope_2020 = sequence_2020
    
    print(f"Envelope 1973: {len(envelope_1973)} bp")
    print(f"Envelope 2020: {len(envelope_2020)} bp")
    
    print("\n[3] Inisialisasi Matriks DP")
    print("-" * 80)
    
    nw = NeedlemanWunschDP(
        seq1=envelope_1973,
        seq2=envelope_2020,
        match_score=MATCH_SCORE,
        mismatch_score=MISMATCH_SCORE,
        gap_penalty=GAP_PENALTY
    )
    
    print("\n[4] Hasil Inisialisasi")
    print("-" * 80)
    
    stats = nw.get_matrix_stats()
    print(f"\nStatistik matriks:")
    print(f"  Dimensi: {stats['rows']} x {stats['cols']}")
    print(f"  Total sel: {stats['total_cells']:,}")
    print(f"  Min/Max: {stats['min_value']} / {stats['max_value']}")
    
    nw.display_dp_matrix(max_rows=12, max_cols=12)
    
    nw.export_to_json(OUTPUT_JSON)
    print("\n[5] Pengisian Matriks DP")
    print("-" * 80)
    print("Mengisi matriks DP...")

    import time
    t0 = time.time()
    final_score = nw.fill_matrix()
    elapsed = time.time() - t0

    print(f"Pengisian matriks selesai dalam {elapsed:.2f} detik.")
    print(f"Skor alignment akhir: {final_score}")

    nw.display_dp_matrix(max_rows=12, max_cols=12)

    print("\n[6] Traceback Alignment")
    print("-" * 80)

    aligned_seq1, aligned_seq2, alignment_bar = nw.traceback()

    stats = nw.get_alignment_stats(aligned_seq1, aligned_seq2, alignment_bar)
    print(f"\nStatistik alignment:")
    print(f"  Panjang alignment : {stats['alignment_length']}")
    print(f"  Match             : {stats['matches']}")
    print(f"  Mismatch          : {stats['mismatches']}")
    print(f"  Gap               : {stats['gaps']}")
    print(f"  Identity          : {stats['identity_pct']}%")
    print(f"  Skor akhir        : {stats['alignment_score']}")

    nw.display_alignment(aligned_seq1, aligned_seq2, alignment_bar)

    nw.export_alignment_to_json(OUTPUT_ALIGNMENT_JSON, aligned_seq1, aligned_seq2, alignment_bar)

    print("\n[7] Analisis Hasil Alignment")
    print("-" * 80)
    mutation_events = build_mutation_events(aligned_seq1, aligned_seq2)
    export_stage_3_analysis(stats, mutation_events, aligned_seq1, aligned_seq2, alignment_bar)
    print(f"Total posisi berbeda: {len(mutation_events)}")
    print(f"Ringkasan hasil untuk laporan: {OUTPUT_REPORT_TXT}")

    return nw, envelope_1973, envelope_2020


if __name__ == "__main__":
    nw, seq1, seq2 = main()
