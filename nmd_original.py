__author__ = 'rhythmicstar'

import gffutils
import pandas as pd


def possible_nmd(nmd_file):
    splicing_data = pd.read_csv(nmd_file, header=None, sep='\s+')
    index = pd.Index(splicing_data[3])
    event_ids = pd.Series(index, name='event_id')
    event_ids_df = event_ids.to_frame()
    event_ids_df['isoform1'] = splicing_data[3].map(lambda x: x.split('|')[0]
                                                    .split('junction:')
                                                    [1][:-2])
    event_ids_df['isoform2'] = splicing_data[3].map(lambda x: x.
                                                    split('isoform2=junction:')
                                                    [1].split('@')[0][:-2])
    event_ids_df['exon'] = splicing_data[3].map(lambda x: x.split('exon:')[1].
                                                split('@')[0][:-2])
    event_ids_df['junction'] = splicing_data[3].map(lambda x: x.
                                                    split('@junction:')[1].
                                                    split(' ')[0][:-2])
    event_ids_df['strand'] = splicing_data[5].map(lambda x: x)
    event_ids_df['type'] = ''
    event_ids_df['inclusion_nmd'] = ''
    event_ids_df['exclusion_nmd'] = ''
    event_ids_df = event_ids_df[['isoform1', 'isoform2', 'exon', 'junction',
                                 'strand', 'type', 'inclusion_nmd',
                                 'exclusion_nmd']]

    # NMD from exon inclusion
    # pos_str = pd.DataFrame([['', '', 'chr3:186506099-186506205', '', '+',
    #                          'inclusion', '', '']], columns=list(
    #     ['isoform1', 'isoform2', 'exon', 'junction', 'strand', 'type',
    #      'inclusion_nmd', 'exclusion_nmd']))
    # neg_str = pd.DataFrame([['', '', 'chr8:109254341-109254575', '', '-',
    #                          'inclusion', '', '']], columns=list(
    #     ['isoform1', 'isoform2', 'exon', 'junction', 'strand', 'type',
    #      'inclusion_nmd', 'exclusion_nmd']))
    # d = pd.DataFrame([['', '', 'chr2:191777919-191778090', '', '+',
    #                    'inclusion', '', '']], columns=list(
    #     ['isoform1', 'isoform2', 'exon', 'junction', 'strand', 'type',
    #      'inclusion_nmd', 'exclusion_nmd']))

    magobhb_inc = pd.DataFrame([['', '', 'chr12:10761697-10761982', '', '+',
                                 'inclusion', '', '']], columns=list(
        ['isoform1', 'isoform2', 'exon', 'junction', 'strand', 'type',
         'inclusion_nmd', 'exclusion_nmd']))

    # asout_df = pd.DataFrame([['', '', 'chr3:186502751-186502890', '', '+',
    #                           'exclusion', '', '']], columns=list(
    #     ['isoform1', 'isoform2', 'exon', 'junction', 'strand', 'type',
    #      'inclusion_nmd', 'exclusion_nmd']))

    # e = pd.DataFrame([['', '', 'chr11:107420479-107420549', '', '-',
    #                    'exclusion', '', '']], columns=list(
    #     ['isoform1', 'isoform2', 'exon', 'junction', 'strand', 'type',
    #      'inclusion_nmd', 'exclusion_nmd']))

    # event_ids_df = event_ids_df.append([pos_str, neg_str, d, new_exon,
    #                                     asout_df, e], ignore_index=True)
    eif4a2_exc = pd.DataFrame([['', '', 'chr3:186502751-186502890', '', '+',
                                'inclusion', '', '']], columns=list(
        ['isoform1', 'isoform2', 'exon', 'junction', 'strand', 'type',
         'inclusion_nmd', 'exclusion_nmd']))
    eif4a2_inc = pd.DataFrame([['', '', 'chr3:186506099-186506205', '', '+',
                                'inclusion', '', '']], columns=list(
        ['isoform1', 'isoform2', 'exon', 'junction', 'strand', 'type',
         'inclusion_nmd', 'exclusion_nmd']))
    MAGOHB_inc = pd.DataFrame([['', '', 'chr12:10761697-10761982', '', '-',
                                'inclusion', '', '']], columns=list(
        ['isoform1', 'isoform2', 'exon', 'junction', 'strand', 'type',
         'inclusion_nmd', 'exclusion_nmd']))
    neg_cont = pd.DataFrame([['', '', 'chr3:186502353-186502486', '', '+',
                              'inclusion', '', '']], columns=list(
        ['isoform1', 'isoform2', 'exon', 'junction', 'strand', 'type',
         'inclusion_nmd', 'exclusion_nmd']))
    event_ids_df = event_ids_df.append([magobhb_inc, eif4a2_exc, eif4a2_inc,
                                        MAGOHB_inc, neg_cont],
                                       ignore_index=True)
    return event_ids_df


def include_exon_nmd(event_ids_df, dbhuman):
    for i, row in event_ids_df.iterrows():
        nmd = False
        location = row['exon']
        for index, human_exon in \
                enumerate(dbhuman.features_of_type('exon', limit=location)):
            if human_exon['transcript_type'] == ["nonsense_mediated_decay"]:
                for stop_codon in dbhuman.children(human_exon['gene_id'][0],
                                                   featuretype='stop_codon'):
                    if stop_codon_nmd(human_exon, stop_codon):
                        nmd = True
        if nmd:
            event_ids_df.ix[i, 6] = 'True'
        else:
            event_ids_df.ix[i, 6] = 'False'
    return event_ids_df


# def stop_codon_nmd(human_exon, stop_codon):
#     if stop_codon.strand == '+':
#         if stop_codon.start >= human_exon.start and stop_codon.end <= \
#                 human_exon.end and stop_codon.end + 50 < human_exon.end:
#             return True
#     if stop_codon.strand == '-':
#         if stop_codon.start >= human_exon.start and stop_codon.end <= \
#                 human_exon.end and human_exon.start + 50 < stop_codon.start:
#             return True
#     return False


def main(database, nmd_file):
    dbhuman = gffutils.FeatureDB(database, keep_order=True)
    event_ids_df = possible_nmd(nmd_file)
    # if negative correlation
    event_ids_df = include_exon_nmd(event_ids_df, dbhuman)
    # if positive correlation
    # event_ids_df = exclude_exon_nmd(event_ids_df)
    return event_ids_df
