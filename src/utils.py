import os
class FASTAReader:
    def __init__(self, filepath):
        self.filepath = filepath
        self.header = None
        self.sequence = None
        self.sequence_length = 0
    
    def read_fasta(self):
        try:
            with open(self.filepath, 'r') as file:
                lines = file.readlines()
            
            if not lines:
                raise ValueError("File FASTA kosong")
            
            if lines[0].startswith('>'):
                self.header = lines[0].strip()
            else:
                raise ValueError("File FASTA tidak valid")
            
            self.sequence = ''.join(line.strip() for line in lines[1:] if line.strip())
            self.sequence_length = len(self.sequence)
            
            return self.header, self.sequence
        
        except FileNotFoundError:
            raise FileNotFoundError(f"File tidak ditemukan: {self.filepath}")
        except Exception as e:
            raise Exception(f"Error membaca file FASTA: {str(e)}")
    
    def get_subsequence(self, start, end):
        if self.sequence is None:
            raise ValueError("Sequence belum dibaca")
        
        start_idx = start - 1
        end_idx = end
        
        if start_idx < 0 or end_idx > len(self.sequence) or start > end:
            raise ValueError(f"Posisi diluar range: {start}-{end}")
        
        return self.sequence[start_idx:end_idx]
