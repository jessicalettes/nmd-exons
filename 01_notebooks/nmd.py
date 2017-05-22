from itertools import product

import pandas as pd
import gffutils
from ordered_set import OrderedSet

COLON = ':'
CDS = 'CDS'
EXCLUSION = 'Exclusion causes NMD'
EXON = 'exon'
GENE_ID = 'gene_id'
INCLUSION = 'Inclusion causes NMD'
MINUS = '-'
NONSENSE_MEDIATED_DECAY = 'nonsense_mediated_decay'
NOTFOUND = 'Exon not found'
PLUS = '+'
START = 'start'
START_CODON = 'start_codon'
STOP_CODON = 'stop_codon'
TRANSCRIPT = 'transcript'
TRANSCRIPT_ID = 'transcript_id'
TRANSCRIPT_TYPE = 'transcript_type'
UNKNOWN = 'Splicing not known to cause NMD'


class NMDExons(object):
    def __init__(self, database, exon_ids):
        """Given exons and a feature database, say which cause NMD and how

        Parameters
        ----------
        database : gffutils.FeatureDB file string
            A database of gene features create from a GTF/GFF file by gffutils
        exon_ids : list-like of strings
            Exons that you want to test whether they could be causing NMD.
            Acceptable format is "exon:chrom:start-stop:strand", e.g.
            'exon:chr1:100-200:+'
        """
        self.db = gffutils.FeatureDB(database)
        self.exon_ids = exon_ids

    def find_nmd_exons(self):
        """Given list of exons, returns if exons are NMD and how

        Returns
        -------
        nmd : Series
            Series contains exon_id string and whether inclusion or exclusion
            of that exon causes NMD
        """
        nmd = pd.Series("Not yet evaluated", index=self.exon_ids)
        counter = 0
        for exon_id in self.exon_ids:
            nmd[exon_id] = self._is_this_exon_nmd(exon_id)
            if counter % 100 == 0:
                print(counter)
            counter += 1
        return nmd

    def _is_this_exon_nmd(self, exon_id):
        """Given a single exon id, returns if it causes NMD and how

        Parameters
        ----------
        exon_id : str
            Single exon ID

        Returns
        -------
        nmd : None | 'inclusion' | 'exclusion'
            If the exon does not cause NMD, return None, else return whether
            inclusion or exclusion of that exon causes NMD
        """
        try:
            exon = self.db[exon_id]
            transcripts_with_exon = self._get_transcripts_with_exon(exon)
            all_transcripts = self._get_all_transcripts_overlapping_exon(exon)
            transcripts_without_exon = all_transcripts - transcripts_with_exon
            exon_strand = exon_id.split(COLON)[3]
            transcript_exons = self._get_exons_from_transcripts(all_transcripts, exon_strand)

            iterator = product(transcripts_without_exon, transcripts_with_exon)
            for t_without, t_with in iterator:
                try:
                    t1_exons_without = transcript_exons[t_without]
                    if len(t1_exons_without) < 5:
                        continue
                except KeyError:
                    return "t_without KeyError ({t_without})".format(t_without=t_without)

                try:
                    t2_exons_with = transcript_exons[t_with]
                except KeyError:
                    return "t_with KeyError ({t_with})".format(t_with=t_with)

                try:
                    if exon_strand == PLUS:
                        first_exon_with = t2_exons_with[2]
                        last_exon_with = t2_exons_with[-1]
                        first_exon_without = t1_exons_without[2]
                        last_exon_without = t1_exons_without[-1]

                    else:
                        first_exon_with = t2_exons_with[0]
                        last_exon_with = t2_exons_with[-3]
                        first_exon_without = t1_exons_without[0]
                        last_exon_without = t1_exons_without[-3]

                    if first_exon_with != exon and last_exon_with != exon:
                        t2_exons_with = list(t2_exons_with)
                        t2_exons_with.remove(first_exon_with)
                        t2_exons_with.remove(last_exon_with)
                        t2_exons_with = tuple(t2_exons_with)
                        t1_exons_without = list(t1_exons_without)
                        t1_exons_without.remove(first_exon_without)
                        t1_exons_without.remove(last_exon_without)
                        t1_exons_without = tuple(t1_exons_without)
                        if (exon_strand == PLUS and t1_exons_without[1] ==
                            t2_exons_with[1]) or \
                                (exon_strand == MINUS and t1_exons_without[-1]
                                    == t2_exons_with[-1]):
                            if len(set(t2_exons_with).symmetric_difference(
                                    set(t1_exons_without))) == 1:
                                t1_type = self.db[t_without][TRANSCRIPT_TYPE][0]
                                t2_type = self.db[t_with][TRANSCRIPT_TYPE][0]
                                if (t1_type != NONSENSE_MEDIATED_DECAY) & \
                                        (t2_type == NONSENSE_MEDIATED_DECAY):
                                    return INCLUSION + " (annotated)"
                                if self._inclusion_nmd(t1_exons_without,
                                                       t2_exons_with, exon_strand):
                                    return INCLUSION + " (found stop codon)"
                                if (t1_type == NONSENSE_MEDIATED_DECAY) & \
                                        (t2_type != NONSENSE_MEDIATED_DECAY):
                                    return EXCLUSION + " (annotated)"
                                if self._exclusion_nmd(t1_exons_without,
                                                       t2_exons_with, exon_strand):
                                    return EXCLUSION + " (found stop codon)"
                except KeyError:
                    return "KeyError"

        except gffutils.FeatureNotFoundError:
            return NOTFOUND
        return UNKNOWN

    def _get_transcripts_with_exon(self, exon):
        """Create set of transcript ids that contain possible NMD exon

        Parameters
        ----------
        exon : gffutils feature exon
            The exon of interest that causes inclusion, exclusion or no NMD

        Returns
        -------
        transcripts_with_exon : set
            The set of transcript ids containing the exon of interest
        """
        transcripts_with_exon = OrderedSet()
        for exon_trans in self.db.parents(exon, featuretype=TRANSCRIPT):
            if self._is_valid_transcript(exon_trans):
                transcripts_with_exon.add(exon_trans[TRANSCRIPT_ID][0])
        return transcripts_with_exon  # parent_transcripts_of_exon

    def _is_valid_transcript(self, transcript):
        contains_cds = self._contains_child(transcript, "CDS")
        contains_utr = self._contains_child(transcript, "UTR")
        contains_start_codon = self._contains_child(transcript, "start_codon")
        return contains_cds and contains_utr and contains_start_codon

    def _contains_child(self, feature, featuretype):
        children = self.db.children(feature, featuretype=featuretype)
        return sum(1 for x in children) > 0

    def _get_all_transcripts_overlapping_exon(self, exon):
        """Makes set of all transcript ids in gene containing possible NMD exon

        Parameters
        ----------
        exon : gffutils feature exon
            The exon of interest that causes inclusion, exclusion or no NMD

        Returns
        -------
        all_transcripts : set
            The set of transcript ids from gene containing the exon of interest
        """
        all_transcripts = OrderedSet()
        for trans in self.db.region(region=exon, featuretype=TRANSCRIPT):
            if self._is_valid_transcript(trans):
                all_transcripts.add(trans[TRANSCRIPT_ID][0])
        return all_transcripts  # call transcripts_from_gene_containing_exon

    def _get_exons_from_transcripts(self, all_transcripts, strand):
        """Create dictionary with tuples containing exons for all transcripts

        Parameters
        ----------
        all_transcripts : set
            The set of transcript ids from gene containing the exon of interest

        strand : str
            Stand containing possible NMD causing exon

        Returns
        -------
        transcript_exons : dict
            Mapping of transcript ids to their children exons, cds and start
            codons. The values are in tuple form with transcript start
            codon, first cds and all internal exons.
        """
        transcript_exons = dict()
        for transcript_id in all_transcripts:
            this_transcript_exons = tuple(self.db.children(transcript_id,
                                                           featuretype=EXON,
                                                           order_by=START))
            transcript_cds = tuple(self.db.children(transcript_id,
                                                    featuretype=CDS,
                                                    order_by=START))
            for start_codon in self.db.children(transcript_id,
                                                featuretype=START_CODON):
                if strand == PLUS:
                    first_cds = transcript_cds[0:1]
                    real_transcript_cds = (start_codon,) + first_cds + this_transcript_exons
                else:
                    first_cds = transcript_cds[-1]
                    real_transcript_cds = this_transcript_exons + (first_cds,) + (start_codon,)
                transcript_exons[transcript_id] = real_transcript_cds
        return transcript_exons

    def _inclusion_nmd(self, trans_without_exon, trans_with_exon, strand):
        """Given transcripts differing by one exon, determine if inclusion NMD

        Parameters
        ----------
        trans_without_exon : tuple
            Identical strand to trans_with_exon but not containing exon
        trans_with_exon : tuple
            Transcript containing possible NMD causing exon
        strand : str
            Stand containing possible NMD causing exon

        Returns
        -------
        bool
            True if inclusion causes NMD, False otherwise
        """
        trans_with_exon_set = OrderedSet(trans_with_exon)
        trans_without_exon_set = OrderedSet(trans_without_exon)
        as_exon = (trans_with_exon_set - trans_without_exon_set)[0]
        index = trans_with_exon.index(as_exon)
        if strand == PLUS:
            for iterator in range(index, len(trans_with_exon)):
                for stop in self.db.children(trans_with_exon[iterator][GENE_ID]
                                             [0], featuretype=STOP_CODON):
                    if self.stop_codon_nmd(trans_with_exon[index], stop):
                        return True
        if strand == MINUS:
            for iterator in range(0, index + 1):
                for stop in self.db.children(trans_with_exon[iterator][GENE_ID]
                                             [0], featuretype=STOP_CODON):
                    if self.stop_codon_nmd(trans_with_exon[index], stop):
                        return True
        return False

    def _exclusion_nmd(self, trans_without_exon, trans_with_exon, strand):
        """Given transcripts differing by one exon, determine if exclusion NMD

        Parameters
        ----------
        trans_without_exon : tuple
            Identical strand to trans_with_exon but not containing exon
        trans_with_exon : tuple
            Transcript containing possible NMD causing exon
        strand : str
            Stand containing possible NMD causing exon

        Returns
        -------
        bool
            True if exclusion causes NMD, False otherwise
        """
        trans_with_exon_set = OrderedSet(trans_with_exon)
        trans_without_exon_set = OrderedSet(trans_without_exon)
        as_exon = (trans_with_exon_set - trans_without_exon_set)[0]
        index = trans_with_exon.index(as_exon)
        if strand == PLUS:
            for iterator in range(index, len(trans_without_exon)):
                for stop in self.db.children(trans_without_exon[iterator]
                                             [GENE_ID][0],
                                             featuretype=STOP_CODON):
                    if self.stop_codon_nmd(trans_without_exon[index], stop):
                        return True
        else:
            for iterator in range(0, index):
                for stop in self.db.children(trans_with_exon[iterator][GENE_ID]
                                             [0], featuretype=STOP_CODON):
                    if self.stop_codon_nmd(trans_with_exon[index], stop):
                        return True
        return False

    @staticmethod
    def stop_codon_nmd(exon, stop_codon):
        """Given an exon and stop_codon, check if exon is NMD

        Parameters
        ----------
        exon : exon gffutils feature object
            The exon to check NMD in
        stop_codon : stop_codon gffutils feature object
            The stop codon in the transcript containing NMD

        Returns
        -------
        bool
            True if exon is NMD and False if not
        """
        if stop_codon.strand == PLUS:
            if stop_codon.start >= exon.start and stop_codon.end <= exon.end \
                    and stop_codon.end + 50 < exon.end:
                return True
        else:
            if stop_codon.start >= exon.start and stop_codon.end <= exon.end \
                    and exon.start + 50 < stop_codon.start:
                return True
        return False
