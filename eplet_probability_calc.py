
import pandas as pd
from collections import defaultdict


# Get top probability from the imputation files
def top_impute_df(top_sl_impute, geno_dict, which_impute):
    top_singleloc = pd.DataFrame()
    for id, values in geno_dict.items():
        top_genotype = max(values, key=values.get)  # only want the haplotype with the highest frequency for now
        top_freq = values[top_genotype]             # only want the highest frequency corresponding with the haplotype
        if which_impute == 'DRDQ':
            line = pd.DataFrame({'DRDQ_Pair': top_genotype, 'DRDQ_PairProb': top_freq}, index=[id])
        elif which_impute == 'DR':
            line = pd.DataFrame({'DR_Pair': top_genotype, 'DR_PairProb': top_freq}, index=[id])
        elif which_impute == 'DQ':
            line = pd.DataFrame({'DQ_Pair': top_genotype, 'DQ_PairProb': top_freq}, index=[id])

        top_singleloc = pd.concat([top_singleloc, line])

    top_sl_impute = pd.concat([top_sl_impute, top_singleloc], axis=1)

    return top_sl_impute

# Create probabilities for eplet MM counts, eplet MM risk categories (Wiebe threshold), and specific eplet MMs
impute_filename = 'lowres_DRDQ_impute.csv'    # Need the probabilites
impute_probs = pd.read_csv(impute_filename, header=0)

# Create a dictionary of all the genotype frequencies with their genotypes
GF_DRDQ = defaultdict(dict)
GF_DR = defaultdict(dict)
GF_DQ = defaultdict(dict)
for line in range(len(impute_probs)):
    subject_id = impute_probs.loc[line, 'ID']
    GENO_DRDQ = impute_probs.loc[line, 'DRDQ_GLString']
    GENO_DR = impute_probs.loc[line, 'DR_GLString']
    GENO_DQ = impute_probs.loc[line, 'DQ_GLString']
    DRDQ_freq = impute_probs.loc[line, 'DRDQ_freq']
    DR_freq = impute_probs.loc[line, 'DR_freq']
    DQ_freq = impute_probs.loc[line, 'DQ_freq']

    # Due to merge in NGS_impute_glstring.py script there are duplicates for the stand alone DR and DQ, so do not do else statement
    if GENO_DRDQ not in GF_DRDQ[subject_id]:
        GF_DRDQ[subject_id][GENO_DRDQ] = DRDQ_freq
    if GENO_DR not in GF_DR[subject_id]:
        GF_DR[subject_id][GENO_DR] = DR_freq
    if GENO_DQ not in GF_DQ[subject_id]:
        GF_DQ[subject_id][GENO_DQ] = DQ_freq

# Calculate the probability for each donor-recipient pairing
pair_computed = {}
DRDQ_donrec_prob = defaultdict(dict)
DR_donrec_prob = defaultdict(dict)
DQ_donrec_prob = defaultdict(dict)
for donor_ID in GF_DRDQ:
    for recip_ID in GF_DRDQ:
        if donor_ID == recip_ID:  # skip pairs with same ID
            continue

        # skip pairs where ID was already computed
        (donor_ID, recip_ID) = sorted([donor_ID, recip_ID])
        id_pair = "+".join([donor_ID, recip_ID])
        if id_pair in pair_computed:
            continue
        else:
            pair_computed[id_pair] = 1

        # loop through all possible genotypes for donor and recip pairing
        for donor_geno_DRDQ in GF_DRDQ[donor_ID]:
            for recip_geno_DRDQ in GF_DRDQ[recip_ID]:
                immunizer_geno_prob = GF_DRDQ[donor_ID][donor_geno_DRDQ]
                patient_geno_prob = GF_DRDQ[recip_ID][recip_geno_DRDQ]

                DR_geno_DRDQ = donor_geno_DRDQ + " " + recip_geno_DRDQ  # Make sure the recip-donor pairs are together
                # Probability for DRDQ together
                if DR_geno_DRDQ not in DRDQ_donrec_prob[id_pair]:
                    DRDQ_donrec_prob[id_pair][DR_geno_DRDQ] = float(immunizer_geno_prob) * float(patient_geno_prob)
        # Do same thing but for DR
        for don_geno_DR in GF_DR[donor_ID]:
            for rec_geno_DR in GF_DR[recip_ID]:
                immunizer_DR_prob = GF_DR[donor_ID][don_geno_DR]
                patient_DR_prob = GF_DR[donor_ID][don_geno_DR]
                DR_geno_DR = don_geno_DR + " " + rec_geno_DR
                if DR_geno_DR not in DR_donrec_prob[id_pair]:
                    DR_donrec_prob[id_pair][DR_geno_DR] = float(immunizer_DR_prob) * float(patient_DR_prob)
        # now for DQ
        for don_geno_DQ in GF_DQ[donor_ID]:
            for rec_geno_DQ in GF_DQ[recip_ID]:
                immunizer_DQ_prob = GF_DQ[donor_ID][don_geno_DQ]
                patient_DQ_prob = GF_DQ[donor_ID][don_geno_DQ]
                DR_geno_DQ = don_geno_DQ + " " + rec_geno_DQ
                if DR_geno_DQ not in DQ_donrec_prob[id_pair]:
                    DQ_donrec_prob[id_pair][DR_geno_DQ] = float(immunizer_DQ_prob) * float(patient_DQ_prob)

DRDQ_top_donrec = pd.DataFrame()
DRDQ_top_donrec = top_impute_df(DRDQ_top_donrec, DRDQ_donrec_prob, 'DRDQ')
DRDQ_top_donrec = top_impute_df(DRDQ_top_donrec, DR_donrec_prob, 'DR')
DRDQ_top_donrec = top_impute_df(DRDQ_top_donrec, DQ_donrec_prob, 'DQ')

# After calculating the probabilites, add them to the eplet MM CSVs with corresponding GLStrings and Pairings
lowres_eplet_filename = 'eplet_lowres_impute.csv'
lowres_eplet = pd.read_csv(lowres_eplet_filename, header=0)  # will concat the genotype probabilities to this DataFrame

truth_eplet_filename = 'eplet_truth_table.csv'
truth_eplet = pd.read_csv(truth_eplet_filename, header=0)
truth_eplet[['REC_ID', 'DON_ID']] = truth_eplet['ID'].str.split("+", expand=True)  # So we can sort on IDs and which goes first
truth_eplet = truth_eplet[['DON_ID', 'REC_ID', 'ALL_quantity', 'ALL_details', 'DRB_quantity', 'DRB_details', 'DQ_quantity', 'DQ_details']]

