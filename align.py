import sys

#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.
    """
    epsilon = 10 ** (-6)
    return abs(a - b) < epsilon


#### ------- CLASSES ------- ####

class MatchMatrix(object):
    """
    Match matrix class stores the scores of matches in a data structure
    """

    def __init__(self):
        # Sequence A = rows
        # Sequence B = cols

        self.matrix = {}

    def set_score(self, a, b, score):
        """
        Updates or adds a score for a specified match

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
           score = the score to set it for
        """
        key = str(a) + str(b)
        self.matrix[key] = float(score)

    def get_score(self, a, b):
        """
        Returns the score for a particular match, where a is the
        character from sequence a and b is from sequence b.

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
        Returns:
           the score of that match
        """
        key = str(a) + str(b)
        return float(self.matrix[key])


class ScoreEntry(object):
    """
    Declare ScoreEntry object, which stores properties of individual cells in a given ScoreMatrix
    """
    def __init__(self, name, row, col, previous, score):
        """
        Initialize variables
        Args:
            name: store the name of the ScoreMatrix that the ScoreEntry is in
            row: store the row in the ScoreMatrix that the ScoreEntry is in
            col: store the col in the ScoreMatrix that the ScoreEntry is in
            previous: store the references to the previous ScoreEntries that set the value of the given ScoreEntry
            score: store the score of the given ScoreEntry
        """
        self.name = name
        self.row = row
        self.col = col
        self.prev = previous  # Tuple of format (Name, row, col)
        self.score = score


class ScoreMatrix(object):
    """
    Object to store a score matrix, which generated during the alignment process. The score matrix consists of a 2-D array of
    ScoreEntries that are updated during alignment and used to output the maximum alignment.
    """

    def __init__(self, name, nrow, ncol):
        """

        Initialize variables
        Args:
            name: identifier for the score matrix - Ix, Iy, or M
            nrow: number of rows in given ScoreMatrix
            ncol: number of cols in given ScoreMatrix

        """
        self.name = name
        self.nrow = nrow
        self.ncol = ncol
        self.score_matrix = list(list())
        for i in range(nrow):
            curr_row = []
            for j in range(ncol):
                curr_row.append(ScoreEntry(name, i, j, list(), None))
            self.score_matrix.append(curr_row)

    def get_score(self, row, col):
        """

        Args:
            row: row of ScoreEntry
            col: col of ScoreEntry

        Returns: ScoreEntry score of the row and col at

        """
        return self.score_matrix[row][col].score

    def set_score(self, row, col, score):
        """

        Args:
            row: row of ScoreEntry
            col: col of ScoreEntry
            score: score of given ScoreEntry

        Sets the score of the given ScoreEntry

        """
        self.score_matrix[row][col].score = float(score)

    def get_pointers(self, row, col):

        """
        Returns the indices of the entries that are pointed to
        """
        return [prev_tuple for prev_tuple in self.score_matrix[row][col].prev]

    def set_pointers(self, row, col, pointer_to_add: ScoreEntry):
        """

        Args:
            row: row of given ScoreEntry
            col: col of given ScoreEntry
            pointer_to_add: reference to ScoreEntry to add as a pointer

        """
        self.score_matrix[row][col].prev.append((pointer_to_add.name, pointer_to_add.row, pointer_to_add.col))

    def print_scores(self):
        """
        Returns a nicely formatted string containing the scores in the score matrix.
        """
        for i in range(len(self.score_matrix)):
            for j in range(len(self.score_matrix[0])):
                print("{:.1f}".format(self.score_matrix[i][j].score) + ",", end=" ")
            print("\n")

    def print_pointers(self):
        """
        Returns a nicely formatted string containing the pointers for each entry in the score matrix.
        """
        for row in range(len(self.score_matrix)):
            for col in range(len(self.score_matrix[0])):
                print(f"row: {row}, col: {col} \t")
                for (name, prev_row, prev_col) in self.score_matrix[row][col].prev:
                    print(f"({name}, {prev_row}, {prev_col}), ")
                print("\n")


class AlignmentParameters(object):
    """
    Object to hold a set of alignment parameters from an input file.
    """

    def __init__(self):
        # default values for variables that are filled in by reading
        # the input alignment file
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.alphabet_a = ""
        self.alphabet_b = ""
        self.len_alphabet_a = 0
        self.len_alphabet_b = 0
        self.match_matrix = MatchMatrix()

    def load_params_from_file(self, input_file):
        """
        Reads the parameters from an input file and stores in the object

        Input:
           input_file = specially formatted alignment input file
        """
        inputfile = open(input_file, 'r')
        self.seq_a = inputfile.readline().strip()
        self.seq_b = inputfile.readline().strip()
        if int(inputfile.readline().strip()) == 0:
            self.global_alignment = True
        self.dx, self.ex, self.dy, self.ey = [float(val) for val in inputfile.readline().strip().split()]
        self.len_alphabet_a = int(inputfile.readline().strip())
        self.alphabet_a = inputfile.readline().strip()
        self.len_alphabet_b = int(inputfile.readline().strip())
        self.alphabet_b = inputfile.readline().strip()

        for i in range(self.len_alphabet_a):
            for j in range(self.len_alphabet_b):
                _, _, a_letter, b_letter, match_score = inputfile.readline().strip().split()
                self.match_matrix.set_score(a_letter, b_letter, float(match_score))
        inputfile.close()


class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by using "align()"
    """

    def __init__(self, input_file, output_file):
        """
        Input:
            input_file = file with the input for running an alignment
            output_file = file to write the output alignments to
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters()

        self.m_matrix = []
        self.ix_matrix = []
        self.iy_matrix = []

    def align(self):
        """
        Main method for running alignment.
        """

        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        # populate the score matrices based on the input parameters
        self.populate_score_matrices()

        # find traceback start
        max_val, max_loc = self.find_traceback_start()

        # perform a traceback and write the output to an output file
        aligned_sequence_list = self.traceback(max_loc)
        self.write_output(aligned_sequence_list, max_val)

    def populate_score_matrices(self):
        """
        Method to populate the score matrices based on the data in align_params.
        """

        # Re-define score matrices based on params
        self.m_matrix = ScoreMatrix("M", len(self.align_params.seq_a) + 1, len(self.align_params.seq_b) + 1)
        self.ix_matrix = ScoreMatrix("Ix", len(self.align_params.seq_a) + 1, len(self.align_params.seq_b) + 1)
        self.iy_matrix = ScoreMatrix("Iy", len(self.align_params.seq_a) + 1, len(self.align_params.seq_b) + 1)

        # Local alignment
        if not self.align_params.global_alignment:
            # Check to make sure match matrix has negative scores
            all_positive = True
            for key in self.align_params.match_matrix.matrix:
                if self.align_params.match_matrix.matrix[key] < 0:
                    all_positive = False
            if all_positive:
                raise ValueError("Error: cannot perform local alignment "
                                 "with only positive values in the match matrix.")
        for i in range(len(self.m_matrix.score_matrix)):
            for j in range(len(self.m_matrix.score_matrix[0])):
                self.update(i, j)

    def update(self, row, col):
        """
        Method to update the matrices at a given row and column index.

        Input:
           row = the row index to update
           col = the column index to update
        """
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)

    def update_m(self, row, col):
        """
        Update a given ScoreEntry in M score matrix using DP
        Args:
            row: the row corresponding to a given ScoreEntry
            col: the col corresponding to a given ScoreEntry

        """
        if row == 0 or col == 0:
            self.m_matrix.score_matrix[row][col].score = 0.0
        else:
            a_char = self.align_params.seq_a[row - 1]
            b_char = self.align_params.seq_b[col - 1]
            # Calculate the score of the given match
            curr_match = self.align_params.match_matrix.get_score(a_char, b_char)
            # Calculate the maximum score from the score matrices
            match = self.m_matrix.score_matrix[row - 1][col - 1].score + curr_match
            a_gap = self.ix_matrix.score_matrix[row - 1][col - 1].score + curr_match
            b_gap = self.iy_matrix.score_matrix[row - 1][col - 1].score + curr_match
            max_val = max(match, max(a_gap, b_gap))
            # Check for local alignment
            if self.align_params.global_alignment:
                self.m_matrix.score_matrix[row][col].score = max_val
            else:
                self.m_matrix.score_matrix[row][col].score = max_val if float(max_val) > 0 else 0

            if fuzzy_equals(max_val, match):
                self.m_matrix.set_pointers(row, col, self.m_matrix.score_matrix[row - 1][col - 1])
            if fuzzy_equals(max_val, a_gap):
                self.m_matrix.set_pointers(row, col, self.ix_matrix.score_matrix[row - 1][col - 1])
            if fuzzy_equals(max_val, b_gap):
                self.m_matrix.set_pointers(row, col, self.iy_matrix.score_matrix[row - 1][col - 1])

    def update_ix(self, row, col):
        """
        Update a given ScoreEntry in the Ix score matrix using DP
        Args:
            row: the row corresponding to a given ScoreEntry
            col: the col corresponding to a given ScoreEntry

        """
        if row == 0 or col == 0:
            self.ix_matrix.score_matrix[row][col].score = 0.0
        else:
            # Calculate the maximum score from the score matrices
            open_a_gap = self.m_matrix.score_matrix[row - 1][col].score - self.align_params.dy
            extend_a_gap = self.ix_matrix.score_matrix[row - 1][col].score - self.align_params.ey
            max_val = max(open_a_gap, extend_a_gap)
            # Check for local alignment
            if self.align_params.global_alignment:
                self.ix_matrix.score_matrix[row][col].score = max_val
            else:
                self.ix_matrix.score_matrix[row][col].score = max_val if float(max_val) > 0 else 0

            if fuzzy_equals(max_val, open_a_gap):
                self.ix_matrix.set_pointers(row, col, self.m_matrix.score_matrix[row - 1][col])
            if fuzzy_equals(max_val, extend_a_gap):
                # Get tuple of (Name, row, col)
                self.ix_matrix.set_pointers(row, col, self.ix_matrix.score_matrix[row - 1][col])

    def update_iy(self, row, col):
        """
        Update a given ScoreEntry in the Iy score matrix using DP
        Args:
            row: the row corresponding to a given ScoreEntry
            col: the col corresponding to a given ScoreEntry

        """
        if row == 0 or col == 0:
            self.iy_matrix.score_matrix[row][col].score = 0.0
        else:
            # Calculate the maximum score from the score matrices
            open_b_gap = self.m_matrix.score_matrix[row][col - 1].score - self.align_params.dx
            extend_b_gap = self.iy_matrix.score_matrix[row][col - 1].score - self.align_params.ex
            max_val = max(open_b_gap, extend_b_gap)
            # Check for local alignment
            if self.align_params.global_alignment:
                self.iy_matrix.score_matrix[row][col].score = max_val
            else:
                self.iy_matrix.score_matrix[row][col].score = max_val if float(max_val) > 0 else 0
            if fuzzy_equals(max_val, open_b_gap):
                # Get tuple of (Name, row, col)
                self.iy_matrix.set_pointers(row, col, self.m_matrix.score_matrix[row][col - 1])
            if fuzzy_equals(max_val, extend_b_gap):
                self.iy_matrix.set_pointers(row, col, self.iy_matrix.score_matrix[row][col - 1])

    def find_traceback_start(self):
        """
        Finds the location to start the traceback.
        """
        max_val = float('-inf')
        max_loc = set()
        if self.align_params.global_alignment:
            # Look only in last row or column of M
            # Last row
            row = len(self.m_matrix.score_matrix) - 1
            for col in range(0, len(self.m_matrix.score_matrix[0])):
                if fuzzy_equals(self.m_matrix.score_matrix[row][col].score, max_val):
                    max_loc.add((row, col))
                elif self.m_matrix.score_matrix[row][col].score > max_val:
                    max_val = self.m_matrix.score_matrix[row][col].score
                    max_loc = set([(row, col)])
            # Last col but exclude last row/col because it has already been checked
            col = len(self.m_matrix.score_matrix[0]) - 1
            for row in range(0, len(self.m_matrix.score_matrix) - 1):
                if fuzzy_equals(self.m_matrix.score_matrix[row][col].score, max_val):
                    max_loc.add((row, col))
                elif self.m_matrix.score_matrix[row][col].score > max_val:
                    max_val = self.m_matrix.score_matrix[row][col].score
                    max_loc = set([(row, col)])
        else:
            # Find max value/locations in M matrix
            # We can't start a traceback in Ix or Iy because then the alignment would be end-gapped
            for row in range(0, len(self.m_matrix.score_matrix)):
                for col in range(0, len(self.m_matrix.score_matrix[0])):
                    max_val = max(max_val, self.m_matrix.score_matrix[row][col].score)
            for row in range(0, len(self.m_matrix.score_matrix)):
                for col in range(0, len(self.m_matrix.score_matrix[0])):
                    if fuzzy_equals(self.m_matrix.score_matrix[row][col].score, max_val):
                        max_loc.add((row, col))
        return float(max_val), max_loc

    def traceback(self, traceback_starts):
        """
        Performs a traceback for each traceback start.

        traceback_starts: Stores a list of row, col tuples that represent the ScoreEntries to start the traceback from

        Returns: A list of alignment sequences (tuples) in the form (a_alignment, b_alignment)

        """
        aligned_sequence_list = []
        if self.align_params.global_alignment:
            for traceback_start in traceback_starts:
                self.rec_traceback(self.m_matrix.score_matrix[traceback_start[0]][traceback_start[1]], "", "", "M",
                                   aligned_sequence_list)
        else:
            # Local alignment ends when you reach a 0
            # If you see a 0 in Ix/Iy AND M only explore M
            for traceback_start in traceback_starts:
                self.rec_traceback(self.m_matrix.score_matrix[traceback_start[0]][traceback_start[1]], "", "", "M",
                                   aligned_sequence_list)
        return aligned_sequence_list

    def rec_traceback(self, curr_score_entry, a_aligned_sequence, b_aligned_sequence, matrix_type,
                      aligned_sequence_list):
        """
        Recursive traceback helper
        Args:
            curr_score_entry: current ScoreEntry to traverse through
            a_aligned_sequence: current a sequence alignment being built up
            b_aligned_sequence: current b sequence alignment being built up
            matrix_type: Identity of the score matrix - M, Ix, Iy
            aligned_sequence_list: a list of all possible alignments being built up

        Returns: All possible alignments from the traceback starts

        """
        # Global
        if self.align_params.global_alignment:
            if curr_score_entry.row == 0 or curr_score_entry.col == 0:
                # We don't want to finish an alignment in Ix or Iy because then it will be start gapped
                if matrix_type == "M":
                    aligned_sequence_list.append((a_aligned_sequence[::-1], b_aligned_sequence[::-1]))
                return
        else:
            # Local alignment
            if curr_score_entry.score == 0:
                # We don't want to finish an alignment in Ix or Iy because then it will be start gapped
                if matrix_type == "M":
                    aligned_sequence_list.append((a_aligned_sequence[::-1], b_aligned_sequence[::-1]))
                return
        for pointer in curr_score_entry.prev:
            # Check what matrix curr_score_entry is in
            if curr_score_entry.name == "M":
                a_char = self.align_params.seq_a[curr_score_entry.row - 1]
                b_char = self.align_params.seq_b[curr_score_entry.col - 1]
            elif curr_score_entry.name == "Ix":
                b_char = "_"
                a_char = self.align_params.seq_a[curr_score_entry.row - 1]
            else:
                b_char = self.align_params.seq_b[curr_score_entry.col - 1]
                a_char = "_"
            if pointer[0] == "M":
                self.rec_traceback(self.m_matrix.score_matrix[pointer[1]][pointer[2]],
                                   a_aligned_sequence + a_char, b_aligned_sequence + b_char, "M",
                                   aligned_sequence_list)
            if pointer[0] == "Ix":
                self.rec_traceback(self.ix_matrix.score_matrix[pointer[1]][pointer[2]],
                                   a_aligned_sequence + a_char, b_aligned_sequence + b_char, "Ix",
                                   aligned_sequence_list)
            if pointer[0] == "Iy":
                self.rec_traceback(self.iy_matrix.score_matrix[pointer[1]][pointer[2]],
                                   a_aligned_sequence + a_char, b_aligned_sequence + b_char, "Iy",
                                   aligned_sequence_list)

    def write_output(self, aligned_sequence_list, max_val):
        """

        Args:
            aligned_sequence_list: A list of tuples of (a_alignments, b_alignments) to write to the output file
            max_val: Maximum value amongst the score matrices

        """
        output_file = open(self.output_file, 'w')
        output_file.write("{:.1f}".format(max_val) + "\n")
        for (a_alignment, b_alignment) in aligned_sequence_list:
            output_file.write("\n")
            output_file.write(a_alignment + '\n')
            output_file.write(b_alignment + '\n')
        output_file.close()


def main():
    # check that the file is being properly used
    if (len(sys.argv) != 3):
        print("Please specify an input file and an output file as args.")
        return

    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()


if __name__ == "__main__":
    main()
