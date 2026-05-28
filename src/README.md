# Tubes KDS Kelompok 07 K1

Deteksi Mutasi Virus Endemik: Implementasi Algoritma Needleman-Wunsch pada Analisis Sekuens DNA Varian Dengue di Indonesia.

## Ringkasan

Program ini melakukan global alignment dua sekuens gen Envelope virus Dengue type 3 menggunakan algoritma Needleman-Wunsch murni Python. Sekuens referensi diambil dari NCBI accession L11425.1 (tahun 1973) dan dibandingkan dengan PQ671078.1 (tahun 2020). Hasil alignment dipakai untuk menghitung identity, jumlah substitusi, insersi, dan delesi.

## Struktur Folder

```
TubesKDS_07_K1/
  dataset/                file FASTA sekuens Dengue
    sequence_1973.fasta
    sequence_2020.fasta
  output/                 hasil tiga tahap
    stage_1_initialization.json
    stage_2_alignment.json
    stage_3_analysis.json
    stage_3_mutation_positions.csv
    stage_3_report_notes.txt
  src/
    main.py               orkestrator pipeline
    needleman_wunsch.py   class NeedlemanWunschDP
    utils.py              class FASTAReader
    pyproject.toml
    README.md
```

## Parameter Scoring

| Parameter | Nilai |
|-----------|-------|
| Match     | +2    |
| Mismatch  | -1    |
| Gap       | -1    |

Region acuan 1973 diambil dari posisi 499 sampai 1977 (1479 bp) agar selaras dengan panjang region Envelope protein E pada sekuens 2020.

## Cara Menjalankan

Menggunakan uv:

```bash
cd src
uv run python main.py
```

Menggunakan Python standar:

```bash
cd src
python main.py
```

Tidak ada dependensi eksternal. Hanya pustaka standar Python 3.13.

## Tahapan Pipeline

1. Pembacaan dua file FASTA
2. Slicing region Envelope 1973
3. Inisialisasi matriks DP
4. Tampilan statistik matriks dan ekspor stage 1
5. Pengisian matriks DP dan pencatatan skor akhir
6. Traceback dan ekspor stage 2
7. Analisis posisi mutasi dan ekspor stage 3

## Sumber Dataset

Sekuens diambil gratis dari NCBI GenBank pada https://www.ncbi.nlm.nih.gov/nuccore. Cari berdasarkan accession number, lalu unduh dalam format FASTA, dan letakkan di folder `dataset/`. Accession yang dipakai:

| File | Accession | Keterangan |
|------|-----------|------------|
| sequence_1973.fasta | L11425.1 | Dengue virus type 3, prM/M dan envelope glycoprotein E polyprotein mRNA, partial cds |
| sequence_2020.fasta | PQ671078.1 | Dengue virus type 3 isolate TMK-PT-001, envelope protein E region, partial cds |

Header lengkap kedua file ditembuskan ke field `metadata.dataset_provenance` pada `stage_1_initialization.json` dan `stage_3_analysis.json` agar laporan dapat mengutip sumber data tanpa parsing manual.

## Output Utama

| File | Isi |
|------|-----|
| stage_1_initialization.json | dimensi matriks, parameter scoring, preview sekuens |
| stage_2_alignment.json      | statistik alignment dan preview hasil |
| stage_3_analysis.json       | ringkasan kuantitatif dan daftar event mutasi |
| stage_3_mutation_positions.csv | tabel posisi mutasi per basa |
| stage_3_report_notes.txt    | ringkasan teks untuk laporan |
