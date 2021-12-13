# Positions to check for mutations from the instructions PDF
position_list = [134, 443, 769, 955, 990, 1051, 1078, 1941, 2138, 2638, 3003]


###############
# IMPORT DATA #
###############
def import_data(file_name):
    import_full = open(file_name, 'r')
    import_full = import_full.readlines()  # Read data as list of lines

    # Find lines with information about sequence length
    lines_length = [ll for ll, lines in enumerate(import_full) if lines.find("Length=") > -1]

    # Find first line with sequence information to get length of sequences in single line
    for line in range(0, len(import_full)):
        find_check = import_full[line].find("Query_")
        if find_check > -1:
            first_query_line = line + 1
            break
    length_sequence_line = len(import_full[first_query_line].split()[2])

    # Create list of amount of lines in a single dataset
    n_data_lines = []
    i = 0
    while i < len(lines_length) - 1:
        n_data_lines.append(lines_length[i + 1] - lines_length[i])
        i += 1

    # Get rid of "\n"s (not necessary but makes data look nicer in the variables)
    import_full = [import_full[i].strip() for i in range(0, len(import_full))]

    # Get sequence length data
    length_data = [import_full[i].split("=")[1] for i in lines_length]

    # Get amount of lines of sequence per patient
    n_sequence_lines = [((int(length_data[i]) // length_sequence_line) * 4) for i in range(0, len(length_data))]

    # Create lists of full sequence/hap_a/hap_b lines
    seqs_full_lists = []  # List of sequences
    hap_a_full_lists = []  # List of mutations in HAP_A
    hap_b_full_lists = []  # List of mutations in HAP_B
    amount_sets = len(length_data)  # Amount of total sequences
    x = 0
    while x < amount_sets:
        for i in range(0, len(import_full), n_data_lines[x]):
            start_seq = i + 24
            seq_single = [import_full[j] for j in range(start_seq, start_seq + n_sequence_lines[x], 4)]
            hap_a_single = [import_full[j + 1] for j in range(start_seq, start_seq + n_sequence_lines[x], 4)]
            hap_b_single = [import_full[j + 2] for j in range(start_seq, start_seq + n_sequence_lines[x], 4)]
            seqs_full_lists.append(seq_single)
            hap_a_full_lists.append(hap_a_single)
            hap_b_full_lists.append(hap_b_single)
            x += 1
    return length_data, seqs_full_lists, hap_a_full_lists, hap_b_full_lists, amount_sets


# Import .dat file and extract data
lengths, seqs, hap_a, hap_b, n_sets = import_data("stu_30.dat")


###############################################
# SPLIT & JOIN HAP_A, HAP_B and SEQUENCE DATA #
###############################################
def combine_data(hap_a_lists, hap_b_lists, sequence_lists):
    hap_a_combined = []
    hap_b_combined = []
    seqs_combined = []
    for i in range(0, len(hap_a_lists)):
        comb_a = []
        comb_b = []
        comb_seqs = []
        j = 0
        while j < len(hap_a_lists[i]):
            comb_a.append(hap_a_lists[i][j].split()[2])  # Get only mutation data
            comb_b.append(hap_b_lists[i][j].split()[2])  # Get only mutation data
            comb_seqs.append(sequence_lists[i][j].split()[2])  # Get only sequence/query data
            j += 1
        hap_a_combined.append(comb_a)
        hap_b_combined.append(comb_b)
        seqs_combined.append(comb_seqs)
    hap_a_joined_list = [None] * len(hap_a_combined)
    hap_b_joined_list = [None] * len(hap_b_combined)
    seqs_joined_list = [None] * len(seqs_combined)
    for i in range(0, len(hap_a_combined)):
        hap_a_joined_list[i] = ''.join(hap_a_combined[i])  # Combine hap_a data into single list element
        hap_b_joined_list[i] = ''.join(hap_b_combined[i])  # Combine hap_b data into single list element
        seqs_joined_list[i] = ''.join(seqs_combined[i])  # # Combine sequence data into single list element
    return hap_a_joined_list, hap_b_joined_list, seqs_joined_list


# Extract data from full lines and combine to be able to easily loop through all positions for comparison
hap_a_joined, hap_b_joined, seqs_joined = combine_data(hap_a, hap_b, seqs)


#######################
# CHECK FOR MUTATIONS #
#######################
def check_for_mutations(positions, amount_sets, hap_a_data, hap_b_data, sequence_data):
    hits_all = []
    x = 0
    while x < amount_sets:  # Counter to go through all the datasets
        c = 0
        y = 0
        hits_single = [0] * len(positions)
        start_position = int(sequence_data[x][0].split()[1])  # Make start position dependent on query start
        positions = [positions[i] - start_position + 1 for i in range(0, len(positions))]
        for i in positions:
            if hap_a_data[x][i - 1] != "." or hap_b_data[x][i - 1] != ".":
                hits_single[y] = 1
            if hap_a_data[x][i - 1] != "." and hap_b_data[x][i - 1] != ".":
                hits_single[y] = 2
            c += 1
            y += 1
        x += 1
        hits_all.append(hits_single)
    return hits_all


# Check for mutations and create nested list for final output
hits_all = check_for_mutations(position_list, n_sets, hap_a_joined, hap_b_joined, seqs)


################
# CREATE TABLE #
################
def create_pat_id(amount_patients, prefix):
    all_pat_ids = []
    for pat_nr in range(1, amount_patients + 1):
        if pat_nr < 10:
            pat_id = [prefix, str(0), str(pat_nr)]
        else:
            pat_id = [prefix, str(pat_nr)]
        pat_id = ''.join(pat_id)
        all_pat_ids.append(pat_id)
    return all_pat_ids


# Create list of patient descriptions from total amount of datasets + prefix
pat_ids = create_pat_id(n_sets, "pat")

print("\t\tT134A\tA443G\tG769C\tG955C\tA990C\tG1051A\tG1078T\tT1941A\tT2138C\tG2638T\tA3003T")
for i in range(0, n_sets):
    print(
        "{}\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}".format(pat_ids[i], hits_all[i][0],
                                                                                    hits_all[i][1], hits_all[i][2],
                                                                                    hits_all[i][3],
                                                                                    hits_all[i][4], hits_all[i][5],
                                                                                    hits_all[i][6],
                                                                                    hits_all[i][7],
                                                                                    hits_all[i][8], hits_all[i][9],
                                                                                    hits_all[i][10]))
