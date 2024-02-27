
import pandas as pd
from collections import defaultdict
import random
import json
import requests
import time

# Create a Monte Carlo Sampling of the Pairs rather than going through each of them
# Select random pairings from the truth table
truth_filename = 'DRDQ_pairs_truth.csv'
truth_pairs = pd.read_csv(truth_filename, header=0)

# Start with 100 pairs as that is what the API can handle


# Run through all possibilities for those pairings in the imputation file
impute_file = 'DRDQ_pairs_imputation.csv'
impute_pairs = pd.read_csv(impute_file, header=0)

# Use weighted choice to select the genotype for that pairing

