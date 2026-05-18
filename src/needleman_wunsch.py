class NeedlemanWunschDP:
    def __init__(self, seq1, seq2, match_score=2, mismatch_score=-1, gap_penalty=-1):
        self.seq1 = seq1.upper()
        self.seq2 = seq2.upper()
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        
        self.m = len(seq1) + 1
        self.n = len(seq2) + 1
        
        self.dp_matrix = self._initialize_dp_matrix()
        self.traceback_matrix = None
    
    def _initialize_dp_matrix(self):
        dp_matrix = [[0 for _ in range(self.n)] for _ in range(self.m)]
        
        for j in range(self.n):
            dp_matrix[0][j] = j * self.gap_penalty
        
        for i in range(self.m):
            dp_matrix[i][0] = i * self.gap_penalty
        
        return dp_matrix
    
    def display_dp_matrix(self, max_rows=12, max_cols=12):
        print(f"\nDP Matrix ({min(max_rows, self.m)} x {min(max_cols, self.n)} dari {self.m} x {self.n}):\n")
        
        header = "      "
        for j in range(min(max_cols, self.n)):
            if j == 0:
                header += "  ∅"
            else:
                header += f"  {self.seq2[j-1]}"
        print(header)
        
        for i in range(min(max_rows, self.m)):
            if i == 0:
                row_label = "  ∅"
            else:
                row_label = f"  {self.seq1[i-1]}"
            
            row_str = f"{row_label:>4}"
            for j in range(min(max_cols, self.n)):
                row_str += f"{self.dp_matrix[i][j]:>4}"
            
            print(row_str)
        
        if max_rows < self.m or max_cols < self.n:
            print(f"\n  [Sisa baris/kolom tidak ditampilkan]")
    
    def get_matrix_stats(self):
        flat_matrix = [val for row in self.dp_matrix for val in row]
        
        stats = {
            'rows': self.m,
            'cols': self.n,
            'total_cells': self.m * self.n,
            'min_value': min(flat_matrix),
            'max_value': max(flat_matrix),
            'seq1_length': len(self.seq1),
            'seq2_length': len(self.seq2),
        }
        
        return stats
    
    def export_to_json(self, filepath):
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

