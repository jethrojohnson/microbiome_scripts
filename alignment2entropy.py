'''
alignment2entropy.py - summarize information in a multiple sequence alignment
=============================================================================

:Author: Jethro Johnson
:Release: $Id$
:Date: |today|
:Tags: Metagenomics

Purpose
-------

When designing primers to work with commonly amplified regions in
metagenomic studies (e.g. 16S, 18S, ITS), it is of interest to know
which subregions contain information useful for distinguishing certain
taxa.

This script takes a multiple sequence alignment and outputs a table
containing the shannon entropy at each base position in the alignment. 

If -r/--variants is specified, output will report whether there is a
substitution (X) or indel (I) at each alignment position.

There is also the option to supply a tab-separated taxonomy file, in
which the first column corresponds to sequences in the alignment, and
subsequent columns correspond to taxonomic levels. If this file is 
present and a taxonomic level specified, then the per-base shannon
entropy is calculated separately for each unique taxa at that taxonomic
level.
'''

import cgatcore.experiment as E
import cgatcore.iotools as IOTools
import pandas as pd
import sys
import re
import numpy as np
from scipy.stats import entropy
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from collections import Counter

def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-a", "--alignment-file", dest="alignment_file",
                      help="an aligned fasta file")
    parser.add_option("-t", "--tax-table", dest="tax_table",
                      help="tab separated table of taxonomy (with headers)")
    parser.add_option("-b", "--baseline", dest="base_sequence",
                      help="name of the baseline sequence against which to"
                      " base comparisons")
    parser.add_option("-l", "--tax-level", dest="tax_level",
                      help="taxonomic level at which to summarize entropy, must"
                      " correspond to a header in the taxonomy file")
    parser.add_option("-m", "--min-occurence", dest="min_replicates",
                       help="minimum number of sequences that must occur"
                       " for entropy to be calculated for a particular taxon")
    parser.add_option("-s", "--sliding-window", dest="sliding_window",
                       help="option supply window for smoothing entropy")
    parser.add_option("-o", "--out-table", dest="out_table",
                       help="location of table to be written out")
    parser.add_option("-r", "--variants", dest="variants", action="store_true",
                      help="if selected, output indels & substitutions rather than entropy")
    parser.set_defaults(base_sequence=None,
                        tax_table=None,
                        tax_level="Phylum",
                        min_replicates=5,
                        sliding_window=None,
                        )

    (options, args) = E.start(parser)

    # read the infiles
    if options.tax_table:
        tax_df = pd.read_table(IOTools.open_file(options.tax_table),
                               header=0,
                               index_col=0)

    align_fa = AlignIO.read(IOTools.open_file(options.alignment_file),
                            format='fasta')

    # fetch the base sequence, if there is one
    if options.base_sequence:
        ids = [x.id for x in align_fa]
        base = align_fa[ids.index(options.base_sequence)]
        index = [i for i, b in enumerate(base.seq.__str__()) if b != '-']
    else:
        index = np.arange(align_fa.get_alignment_length())
    
    # set up dataframe to contain entropy calculations
    if options.tax_table and options.tax_level:
        columns = ['Entropy',]
        for tax in set(tax_df[options.tax_level]):
            n = tax_df[options.tax_level].tolist()
            if n.count(tax) >= int(options.min_replicates):
                columns.append(tax)
            else: 
                E.warn('Skipping taxon %s... too few sequences (%i) in the'
                       ' taxonomy file' % (tax, n.count(tax)))
    else:
        columns = ['Entropy',]

    df_entropy = pd.DataFrame(columns=columns, index=index)

    # write entropy to dataframe
    to_drop = []
    for tax in df_entropy.columns.tolist():
        if tax == 'Entropy':
            sequence_ids = [x.id for x in align_fa]
        else:
            sequence_ids = tax_df[tax_df[options.tax_level] \
                                      == tax].index.tolist()

        seq_records = []
        n = 0
        for sequence in align_fa:
            if sequence.id in sequence_ids:
                seq_records.append(sequence)
                n += 1
            elif options.base_sequence and \
                    sequence.id == options.base_sequence:
                seq_records.append(sequence)
            else:
                continue

        E.info('There are %i sequences in alignment for %s' % (n, tax))

        # drop taxa for which there are too few sequences in alignment
        # (...catches cases where sequences in tax file are absent in
        # alignment)
        if len(seq_records) < int(options.min_replicates):
            E.warn('Skipping taxon %s... too few sequences (%i) in the'
                   ' alignment file' % (tax, len(seq_records)))
            to_drop.append(tax)
            continue
        
        # Iterate across alignment and calculate either indel/subs or entropy
        sub_fa = MultipleSeqAlignment(seq_records)

        for i in df_entropy.index.tolist():
            string = sub_fa[:,i]
            assert len(set('ACTGU-') - set(string)) >= 0
            
            string = Counter(sub_fa[:,i])
            if options.variants:
                if len(string.keys()) > 1:
                    if len(string.keys()) > 2 and string.get('-'):
                        ent = "I,X"
                    elif len(string.keys()) == 2 and string.get('-'):
                        ent = "I"
                    else:
                        ent = "X"
                else:
                    ent = ""
                    
            else:
                ent = entropy([x/ float(len(string)) for x in string.values()])
            
            df_entropy.loc[i, tax] = ent

    # drop taxa & reset index to match positions in base sequence
    df_entropy.drop(to_drop, axis=1, inplace=True)
    df_entropy.reset_index(drop=True, inplace=True)

    # smooth entropy
    if options.sliding_window:
        if options.variants:
            raise IOError('Cannot select sliding window and variants')
        df_entropy = pd.rolling_mean(df_entropy,
                                     center=True,
                                     window=int(options.sliding_window))

    df_entropy.to_csv(IOTools.open_file(options.out_table, 'w'),
                      sep='\t',
                      index_label='Base_Position')

if __name__ == "__main__":
    sys.exit(main(sys.argv))
