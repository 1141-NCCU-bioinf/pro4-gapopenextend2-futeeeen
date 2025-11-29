import sys

def parse_fasta(file_path):
    """
    Parses a FASTA file and returns a list of (header, sequence) tuples.
    Replaces the functionality of Biostrings for Python compatibility.
    """
    sequences = []
    current_header = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line[1:] # Remove >
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences.append((current_header, "".join(current_seq)))
            
    return sequences

def read_pam_matrix(file_path):
    """
    Parses the PAM/BLOSUM score matrix file.
    Returns a dictionary {(aa1, aa2): score} and the list of amino acids.
    """
    scores = {}
    amino_acids = []
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        # Find the header line (starts with a space or first AA)
        header_idx = 0
        for i, line in enumerate(lines):
            if line.strip().startswith('#') or not line.strip():
                continue
            # Assume the first non-comment line with multiple chars is header
            parts = line.split()
            if len(parts) > 0:
                amino_acids = parts
                header_idx = i + 1
                break
                
        # Parse rows
        for i in range(header_idx, len(lines)):
            parts = lines[i].split()
            if not parts:
                continue
            row_aa = parts[0]
            row_scores = parts[1:]
            
            for j, score in enumerate(row_scores):
                if j < len(amino_acids):
                    col_aa = amino_acids[j]
                    scores[(row_aa, col_aa)] = int(score)
                    
    return scores

def alignment(input_path, score_path, output_path, aln, gap_open, gap_extend):
    """
    Performs pairwise alignment (Global or Local) with affine gap penalties.
    """
    
    # 1. Load Data
    seq_data = parse_fasta(input_path)
    if len(seq_data) < 2:
        print("Error: Input FASTA must contain at least two sequences.")
        return
        
    name1, seq1 = seq_data[0]
    name2, seq2 = seq_data[1]
    
    n = len(seq1)
    m = len(seq2)
    
    score_matrix = read_pam_matrix(score_path)
    
    # Constants for infinity
    # Using a large negative number that won't underflow easily like float('-inf') in integer math
    MIN_SCORE = -999999999 

    # 2. Initialize Matrices
    # We store tuples: (score, alignment_length) to handle the tie-breaking requirement.
    # M: Match/Mismatch state
    # X: Gap in seq1 (seq2 has character, seq1 has gap)
    # Y: Gap in seq2 (seq1 has character, seq2 has gap)
    
    M = [[(MIN_SCORE, 0) for _ in range(m + 1)] for _ in range(n + 1)]
    X = [[(MIN_SCORE, 0) for _ in range(m + 1)] for _ in range(n + 1)]
    Y = [[(MIN_SCORE, 0) for _ in range(m + 1)] for _ in range(n + 1)]
    
    # Traceback matrices: 0=Stop/End, 1=Diagonal(Match), 2=Up(GapY), 3=Left(GapX)
    # We need separate traceback pointers for affine states, but for simplicity,
    # we will deduce direction based on scores during backtracking or use a simplified state tracker.
    # To strictly follow "longest alignment" during backtracking, we rely on the (score, length) stored.
    
    is_global = (aln.lower() == 'global')
    
    # Initialization
    M[0][0] = (0, 0)
    
    if is_global:
        # Global initialization
        for i in range(1, n + 1):
            # Cost to open at start + extend
            cost = gap_open + (i - 1) * gap_extend
            M[i][0] = (cost, i)
            Y[i][0] = (cost, i) # Technically implies gap in seq2
            X[i][0] = (MIN_SCORE, 0) # Cannot have gap in seq1 at col 0 boundaries in global usually

        for j in range(1, m + 1):
            cost = gap_open + (j - 1) * gap_extend
            M[0][j] = (cost, j)
            X[0][j] = (cost, j) # Gap in seq1
            Y[0][j] = (MIN_SCORE, 0)
    else:
        # Local initialization: First row and column remain 0 (or MIN_SCORE logically, but treated as 0 in recursion)
        # Actually, in Local, M[i][0] and M[0][j] are 0.
        for i in range(n+1): M[i][0] = (0, 0)
        for j in range(m+1): M[0][j] = (0, 0)

    # Variables to track max score for Local Alignment
    max_local_score = (MIN_SCORE, 0)
    max_local_pos = (0, 0)

    # 3. DP Calculation
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s_char = seq_matrix_score(seq1[i-1], seq2[j-1], score_matrix)
            
            # --- Calculate X (Gap in Seq1 / Insertion in Seq2) ---
            # Extend existing gap in X or open new gap from M
            # Format: (score, length)
            
            # From M (Open gap): score = M_score + open, len = M_len + 1
            x_from_m = (M[i][j-1][0] + gap_open, M[i][j-1][1] + 1)
            
            # From X (Extend gap): score = X_score + extend, len = X_len + 1
            x_from_x = (X[i][j-1][0] + gap_extend, X[i][j-1][1] + 1)
            
            # In some definitions, you can move Y->X, but usually X handles horizontal, Y vertical
            X[i][j] = max(x_from_m, x_from_x)
            
            # --- Calculate Y (Gap in Seq2 / Insertion in Seq1) ---
            y_from_m = (M[i-1][j][0] + gap_open, M[i-1][j][1] + 1)
            y_from_y = (Y[i-1][j][0] + gap_extend, Y[i-1][j][1] + 1)
            
            Y[i][j] = max(y_from_m, y_from_y)
            
            # --- Calculate M (Match/Mismatch) ---
            # Match comes from diagonal (i-1, j-1)
            m_from_m = (M[i-1][j-1][0] + s_char, M[i-1][j-1][1] + 1)
            m_from_x = (X[i-1][j-1][0] + s_char, X[i-1][j-1][1] + 1)
            m_from_y = (Y[i-1][j-1][0] + s_char, Y[i-1][j-1][1] + 1)
            
            best_m = max(m_from_m, m_from_x, m_from_y)
            
            if not is_global:
                # Local alignment: Clamp negative scores to 0
                if best_m[0] < 0:
                    best_m = (0, 0)
            
            M[i][j] = best_m
            
            # Track Max for Local
            if not is_global:
                # Python tuple comparison checks score first, then length
                # This automatically satisfies: "output one with the longest alignment length"
                if M[i][j] > max_local_score:
                    max_local_score = M[i][j]
                    max_local_pos = (i, j)

    # 4. Traceback
    aln1 = []
    aln2 = []
    
    if is_global:
        i, j = n, m
        # For global, we assume we end in M state logic usually, 
        # but check which matrix has the best final score
        # (Though biologically usually ends in Match state logic)
        curr_state = 'M' 
        if X[n][m] > M[n][m] and X[n][m] > Y[n][m]: curr_state = 'X'
        elif Y[n][m] > M[n][m] and Y[n][m] > X[n][m]: curr_state = 'Y'
    else:
        i, j = max_local_pos
        curr_state = 'M' # Local starts tracking back from the max score cell in M
        print(f"Max Score: {max_local_score[0]}") # Outputting threshold/score

    while (i > 0 or j > 0):
        if not is_global and M[i][j][0] == 0 and curr_state == 'M':
            break # Stop local alignment when we hit 0
            
        if is_global and i == 0 and j == 0:
            break

        if curr_state == 'M':
            # Determine where we came from based on re-calculation
            # We look for the neighbor that could produce M[i][j]
            
            # Current score
            curr_score = M[i][j][0]
            s_char = seq_matrix_score(seq1[i-1], seq2[j-1], score_matrix)
            
            # Candidates
            # Note: We must check 'gap_open' usage carefully.
            # Usually M is formed by adding score to prev diagonal state.
            
            # Check if came from M(i-1, j-1)
            cand_m = M[i-1][j-1][0] + s_char
            # Check if came from X(i-1, j-1)
            cand_x = X[i-1][j-1][0] + s_char
            # Check if came from Y(i-1, j-1)
            cand_y = Y[i-1][j-1][0] + s_char
            
            # We assume precedence M > X > Y if scores are equal, 
            # OR we follow the length logic.
            # M[i][j] stores (score, len). Let's see which neighbor + 1 matches this len.
            target_len = M[i][j][1]
            
            # Prioritize matching both score AND length logic
            if i > 0 and j > 0 and abs(cand_m - curr_score) < 1e-9 and M[i-1][j-1][1] + 1 == target_len:
                aln1.append(seq1[i-1])
                aln2.append(seq2[j-1])
                i -= 1
                j -= 1
                curr_state = 'M'
            elif i > 0 and j > 0 and abs(cand_x - curr_score) < 1e-9 and X[i-1][j-1][1] + 1 == target_len:
                aln1.append(seq1[i-1])
                aln2.append(seq2[j-1])
                i -= 1
                j -= 1
                curr_state = 'X'
            elif i > 0 and j > 0 and abs(cand_y - curr_score) < 1e-9 and Y[i-1][j-1][1] + 1 == target_len:
                aln1.append(seq1[i-1])
                aln2.append(seq2[j-1])
                i -= 1
                j -= 1
                curr_state = 'Y'
            else:
                # Fallback for boundaries or if floating point logic gets weird (unlikely with ints)
                # Force global boundary handling
                if i > 0 and j > 0:
                    # Prefer match
                    aln1.append(seq1[i-1])
                    aln2.append(seq2[j-1])
                    i -= 1
                    j -= 1
                    curr_state = 'M' 
                elif j > 0: # Only col remains
                    curr_state = 'X'
                elif i > 0: # Only row remains
                    curr_state = 'Y'

        elif curr_state == 'X':
            # X[i][j] comes from M[i][j-1] + open OR X[i][j-1] + extend
            curr_score = X[i][j][0]
            target_len = X[i][j][1]
            
            cand_m_open = M[i][j-1][0] + gap_open
            cand_x_extend = X[i][j-1][0] + gap_extend
            
            if abs(cand_x_extend - curr_score) < 1e-9 and X[i][j-1][1] + 1 == target_len:
                aln1.append('-')
                aln2.append(seq2[j-1])
                j -= 1
                curr_state = 'X'
            else:
                aln1.append('-')
                aln2.append(seq2[j-1])
                j -= 1
                curr_state = 'M'

        elif curr_state == 'Y':
            # Y[i][j] comes from M[i-1][j] + open OR Y[i-1][j] + extend
            curr_score = Y[i][j][0]
            target_len = Y[i][j][1]
            
            cand_m_open = M[i-1][j][0] + gap_open
            cand_y_extend = Y[i-1][j][0] + gap_extend
            
            if abs(cand_y_extend - curr_score) < 1e-9 and Y[i-1][j][1] + 1 == target_len:
                aln1.append(seq1[i-1])
                aln2.append('-')
                i -= 1
                curr_state = 'Y'
            else:
                aln1.append(seq1[i-1])
                aln2.append('-')
                i -= 1
                curr_state = 'M'
                
    # 5. Output Results
    final_seq1 = "".join(reversed(aln1))
    final_seq2 = "".join(reversed(aln2))
    
    with open(output_path, 'w') as f:
        f.write(f">{name1}\n")
        f.write(f"{final_seq1}\n")
        f.write(f">{name2}\n")
        f.write(f"{final_seq2}\n")
        
    if is_global:
         print(f"Global Alignment Score: {M[n][m][0]}")

def seq_matrix_score(a, b, matrix):
    # Handle cases where key might be (A, B) or (B, A)
    if (a, b) in matrix:
        return matrix[(a, b)]
    elif (b, a) in matrix:
        return matrix[(b, a)]
    else:
        # Fallback for unknown characters (optional)
        return -1 # Penalty for unknown

if __name__ == "__main__":
    # Example usage for testing if run directly
    #alignment("examples/test.fasta", "examples/pam250.txt", "result_global.fasta", "global", -10, -2)
    pass