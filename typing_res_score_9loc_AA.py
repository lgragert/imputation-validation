#!/usr/bin/env python

import gzip
from collections import defaultdict
import sys
import hlagenie
genie = hlagenie.init("3510")

pop = sys.argv[1]  # Loop through population groups in slurm script: AFA ASN CAU HIS MLT NAM
pops = [pop]  # ['AFA', 'ASN', 'CAU', 'HIS', 'MLT', 'NAM']


loci = ["A", "B", "C", "DRB345", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]

full_start_pos = {
    "A" : 1,
    "B" : 1,
    "C" : 1,
    "DRB345" : 1,
    "DRB1" : 1,
    "DQA1" : 1,
    "DQB1" : 1,
    "DPA1" : 1,
    "DPB1" : 1,
}

full_end_pos = {
    "A": 341,
    "B": 338,
    "C": 342,
    "DRB1": 237,
    "DRB345": 237,
    "DQA1": 232,
    "DQB1": 229,
    "DPA1": 229,
    "DPB1": 229,
}


# Utility function to create dictionary
def multi_dict(K, type):
    if K == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: multi_dict(K-1, type))


happair_probs = {}  # HLA probability distribution
happair_hla = {}  # HLA haplotype pair distribution
happair_id_total = {}  # cumulative genotype frequency total
subject_ID_ethnicity_study = {}
for popn in pops:
    # load imputation output file and load probabilities for all subjects
    impute_outfilename = "./impute." + popn + ".csv.gz"
    impute_outfile = gzip.open(impute_outfilename, "rt")

    # compute cumulative genotype frequency totals per subject
    for line in impute_outfile:
        (subject_id, rank, hap1, hap2, prob) = line.split(',')
        if (subject_id == "PX_ID"):  # skip header row
            continue

        happair_freq = 0

        subject_ID_ethnicity_study[subject_id] = popn

        if subject_id not in happair_id_total:
            happair_id_total[subject_id] = 0
        happair_id_total[subject_id] += float(prob)

    impute_outfile.close()

    impute_outfile = gzip.open(impute_outfilename, "rt")

    # compute probabilties for each haplotype pair with normalization
    for line in impute_outfile:
        (subject_id, rank, hap1, hap2, prob) = line.split(',')
        if (subject_id == "PX_ID"):  # skip header row
            continue

        happair = hap1 + "+" + hap2
        happair_freq = float(prob)

        if subject_id not in happair_probs:
            happair_probs[subject_id] = list()

        # compute probability after renormalization
        happair_probs[subject_id].append(happair_freq / happair_id_total[subject_id])

        if subject_id not in happair_hla:
            happair_hla[subject_id] = list()
        happair_hla[subject_id].append(happair)
        
    impute_outfile.close()

exit

print ("Haplotype pair probabilities loaded")

# compute probabilities for each AA-level genotype
TRS_dict = multi_dict(3, float)  # TRS per subject per locus per position
subject_counter = 0
for subject_id in happair_hla:

    happair_index = 0
    AA_geno_probs = multi_dict(3, float)  # AA-level genotype list and probs per locus position
    for happair in happair_hla[subject_id]:

        (hap1, hap2) = happair.split('+')

        (a1, c1, b1, dr345_1, drb1_1, dqa1_1, dqb1_1, dpa1_1, dpb1_1) = hap1.split('~')
        (a2, c2, b2, dr345_2, drb1_2, dqa1_2, dqb1_2, dpa1_2, dpb1_2) = hap2.split('~')

        if a1 == "A*11:52":
            a1 = "A*11:52Q"
        if a2 == "A*11:52":
            a2 = "A*11:52Q"
        if a1 == "A*23:19Q":
            a1 = "A*23:19N"
        if a2 == "A*23:19Q":
            a2 = "A*23:19N"

        if b1 == "B*13:08Q":
            b1 = "B*13:08"
        if b2 == "B*13:08Q":
            b2 = "B*13:08"
        if b1 == "B*15:22":
            b1 = "B*35:43"
        if b2 == "B*15:22":
            b2 = "B*35:43"

        if c1 == "C*15:20":
            c1 = "C*15:27"
        if c2 == "C*15:20":
            c2 = "C*15:27"
        if c1 == "C*03:12":
            c1 = "C*03:19"
        if c2 == "C*03:12":
            c2 = "C*03:19"
        if c1 == "C*03:23":
            c1 = "C*03:23N"
        if c2 == "C*03:23":
            c2 = "C*03:23N"
        if c1 == "C*05:02":
            c1 = "C*05:09"
        if c2 == "C*05:02":
            c2 = "C*05:09"
        if c1 == "C*13:01":
            c1 = "C*12:02"
        if c2 == "C*13:01":
            c2 = "C*12:02"

        if drb1_1 == "DRB1*07:02":
            drb1_1 = "DRB1*07:01"
        if drb1_2 == "DRB1*07:02":
            drb1_2 = "DRB1*07:01"

        if dqa1_1 == "DQA1*01:07":
            dqa1_1 = "DQA1*01:07Q"
        if dqa1_2 == "DQA1*01:07":
            dqa1_2 = "DQA1*01:07Q"

        happair_prob_list = happair_probs[subject_id]

        happair_prob = happair_probs[subject_id][happair_index]
        happair_index += 1

        for locus in loci:

            if (locus == "A"):
                allele1 = a1
                allele2 = a2
            if (locus == "B"):
                allele1 = b1
                allele2 = b2
            if (locus == "C"):
                allele1 = c1
                allele2 = c2
            if (locus == "DRB345"):
                allele1 = dr345_1
                allele2 = dr345_2
            if (locus == "DRB1"):
                allele1 = drb1_1
                allele2 = drb1_2
            if (locus == "DQA1"):
                allele1 = dqa1_1
                allele2 = dqa1_2
            if (locus == "DQB1"):
                allele1 = dqb1_1
                allele2 = dqb1_2
            if (locus == "DPA1"):
                allele1 = dpa1_1
                allele2 = dpa1_2
            if (locus == "DPB1"):
                allele1 = dpb1_1
                allele2 = dpb1_2

            if subject_id not in happair_probs:
                happair_probs[subject_id] = list()

            # get appropriate position range per locus
            # lots of incomplete sequences in IMGT/HLA that are filled in - use only ARD positions
            for position in range(full_start_pos[locus], full_end_pos[locus]):
                if allele1 == "DRBX*NNNN":
                    AA1 = 'None'
                else:
                    AA1 = genie.getAA(allele1, position)

                if allele2 == "DRBX*NNNN":
                    AA2 = 'None'
                else:
                    AA2 = genie.getAA(allele2, position)

                (AA1, AA2) = sorted([AA1, AA2])  # alpha sort positions
                AA_geno = AA1 + "+" + AA2

                AA_geno_probs[locus][position][AA_geno] += happair_prob

    for locus in loci:

        for position in range(full_start_pos[locus], full_end_pos[locus]):
            TRS = 0

            for AA_geno in AA_geno_probs[locus][position]:
                AA_geno_prob = AA_geno_probs[locus][position][AA_geno]
                AA_geno_prob = round(AA_geno_prob, 10)

                TRS = TRS + (AA_geno_prob * AA_geno_prob)
                TRS_increment = AA_geno_prob * AA_geno_prob
            TRS_dict[subject_id][locus][position] = TRS

    subject_counter += 1
    if ((subject_counter % 10000) == 0):
        print("Number of subjects with TRS computed: " + str(subject_counter))

# print out summary of all TRS values
TRS_average = multi_dict(4, float)  # Average TRS per donor/recip per race per locus per position


nsubjects = len(subject_ID_ethnicity_study)
print ("Number of subjects in the impute." + pop + ".csv.gz file: " + str(nsubjects))

# compute the number of subjects by category to get denominators for averages
nsubject_ethnicity = defaultdict(int)
nsubject_donor_ethnicity = defaultdict(int)
nsubject_recip_ethnicity = defaultdict(int)
nsubject_donor = 0
nsubject_recip = 0

for subject_id in subject_ID_ethnicity_study:

    # ethnicity
    subject_ethnicity = subject_ID_ethnicity_study[subject_id]
    nsubject_ethnicity[subject_ethnicity] += 1

    # donor vs recip
    if subject_id.startswith("D"):
        nsubject_donor += 1
        nsubject_donor_ethnicity[subject_ethnicity] += 1
    else:
        nsubject_recip += 1
        nsubject_recip_ethnicity[subject_ethnicity] += 1

ethnicity_list = list(nsubject_ethnicity.keys())

print (ethnicity_list)
print ("Number of donors: " + str(nsubject_donor))
print ("Number of recipients: " + str(nsubject_recip))
print ("Number of " + pop + " donors: " + str(nsubject_donor_ethnicity[pop]))
print ("Number of " + pop + " recips: " + str(nsubject_recip_ethnicity[pop]))

# compute averages TRS per position for donors/recips
# compute averages TRS per position by race/ethnicity
# compute average TRS per position across the total multiethnic dataset
for subject_id in subject_ID_ethnicity_study:

    donor_or_recip = ""
    if subject_id.startswith("D"):
        donor_or_recip = "DONOR"
    else:
        donor_or_recip = "RECIP"

    ethnicity = subject_ID_ethnicity_study[subject_id]

    for locus in loci:
        # for locus in loci_selected:
        for position in range(full_start_pos[locus], full_end_pos[locus]):
            # for position in DRB1_AA_positions_selected:
            TRS = TRS_dict[subject_id][locus][position]

            TRS_average["SUBJECT"][ethnicity][locus][position] += (TRS / nsubject_ethnicity[ethnicity])
            if (donor_or_recip == "DONOR"):

                TRS_average["DONOR"][ethnicity][locus][position] += (TRS / nsubject_donor_ethnicity[ethnicity])
            else:
                TRS_average["RECIP"][ethnicity][locus][position] += (TRS / nsubject_recip_ethnicity[ethnicity])


# print summary table
TRS_average_out_filename = "HLA_AA_TRS_Average_" + pop + ".csv"
TRS_average_out_file = open(TRS_average_out_filename, "w")
TRS_average_out_file.write("Subject_Type,Ethnicity,Locus,AA_Position,TRS_Average\n")
for donor_or_recip in ["SUBJECT", "DONOR", "RECIP"]:
    for ethnicity in ethnicity_list:
        for locus in loci:
            for position in range(full_start_pos[locus], full_end_pos[locus]):
                print(donor_or_recip + " " + ethnicity + " " + locus + " " + str(position))

                # UNK recip has 0 cases, so TRS should be NA
                TRS = ""
                if (not TRS_average[donor_or_recip][ethnicity][locus][position]):
                    TRS = "NA"
                else:
                    TRS = str(round(TRS_average[donor_or_recip][ethnicity][locus][position], 8))

                print("Average TRS: " + TRS)
                TRS_out_string = ",".join([donor_or_recip, ethnicity, locus, str(position), TRS])
                TRS_average_out_file.write(TRS_out_string + "\n")

# TODO - create separate table per locus
# rows in each tables are donor/recip-ethnicity combo
# columns are positions

TRS_average_out_file.close()

exit()
