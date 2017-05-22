import gffutils
import pandas.util.testing as pdt
import pandas as pd
import pytest


@pytest.fixture
def database():
    return '/Users/rhythmicstar/projects/exon_evolution//gencode.v19.' \
           'annotation.outrigger.nmdtest.gtf.db'


@pytest.fixture
def exon_ids():
    return ('exon:chr10:101510126-101510153:+',
            'exon:chr7:33075546-33075600:-',
            'exon:chr12:42778742-42778798:+',
            'exon:chr8:29931393-29931571:-',
            'exon:chr3:186502353-186502486:+',
            'exon:chr3:42661508-42661535:+',
            'exon:chr14:20785953-20786133:-')


@pytest.fixture
def single_exon_id():
    return 'exon:chr3:186502751-186502890:+'


@pytest.fixture
def is_exon_nmd():
    return 'Exclusion causes NMD (annotated)', \
           'Inclusion causes NMD (annotated)', \
           'Exclusion causes NMD (found stop codon)', \
           'Inclusion causes NMD (found stop codon)', \
           'Splicing not known to cause NMD', \
           'Splicing not known to cause NMD', \
           'Splicing not known to cause NMD'


@pytest.fixture()
def stop_codon_exon_ids():
    return [('stop_codon:chr3:186506106-186506108:+',
             'exon:chr3:186506099-186506205:+', True),
            ('stop_codon:chr17:44101535-44101537:+',
            'exon:chr17:44101322-44101549:+', False)]


@pytest.fixture()
def parent_transcripts_of_exon():
    return {'ENST00000323963.5', 'ENST00000498746.1', 'ENST00000440191.2',
            'ENST00000425053.1'}


@pytest.fixture()
def all_transcripts_of_exon():
    return {'ENST00000497177.1', 'ENST00000429589.1', 'ENST00000486805.1',
            'ENST00000466362.1', 'ENST00000441007.1', 'ENST00000475653.1',
            'ENST00000465792.1', 'ENST00000440191.2', 'ENST00000445596.1',
            'ENST00000465222.1', 'ENST00000425053.1', 'ENST00000492144.1',
            'ENST00000494445.1', 'ENST00000467585.1', 'ENST00000498746.1',
            'ENST00000461021.1', 'ENST00000465032.1', 'ENST00000443963.1',
            'ENST00000426808.1', 'ENST00000495049.1', 'ENST00000491473.1',
            'ENST00000323963.5', 'ENST00000496382.1', 'ENST00000468362.1',
            'ENST00000475409.1', 'ENST00000465267.1', 'ENST00000485101.1',
            'ENST00000356531.5'}


@pytest.fixture()
def inc_trans_without_exon_ids():
    return ['start_codon:chr11:57480091-57480093:+',
            'CDS:chr11:57480091-57480279:+:0',
            'exon:chr11:57505080-57505140:+',
            'exon:chr11:57505826-57505902:+',
            'exon:chr11:57506136-57506242:+',
            'exon:chr11:57506446-57506511:+',
            'exon:chr11:57506603-57506732:+']



@pytest.fixture()
def inc_trans_with_exon_ids():
    return ['start_codon:chr11:57480091-57480093:+',
            'CDS:chr11:57480091-57480279:+:0',
            'exon:chr11:57505080-57505140:+',
            'exon:chr11:57505385-57505498:+',
            'exon:chr11:57505826-57505902:+',
            'exon:chr11:57506136-57506242:+',
            'exon:chr11:57506446-57506511:+',
            'exon:chr11:57506603-57506732:+']


@pytest.fixture()
def exc_trans_without_exon_ids():
    return ['start_codon:chr12:42729705-42729707:+',
            'CDS:chr12:42729705-42729776:+:0',
            'exon:chr12:42729685-42729776:+',
            'exon:chr12:42745687-42745851:+',
            'exon:chr12:42748963-42749024:+',
            'exon:chr12:42768665-42768876:+',
            'exon:chr12:42781258-42781337:+',
            'exon:chr12:42787372-42787491:+',
            'exon:chr12:42792656-42792796:+']

@pytest.fixture()
def exc_trans_with_exon_ids():
    return ['start_codon:chr12:42729705-42729707:+',
            'CDS:chr12:42729705-42729776:+:0',
            'exon:chr12:42729685-42729776:+',
            'exon:chr12:42745687-42745851:+',
            'exon:chr12:42748963-42749024:+',
            'exon:chr12:42768665-42768876:+',
            'exon:chr12:42778742-42778798:+',
            'exon:chr12:42781258-42781337:+',
            'exon:chr12:42787372-42787491:+',
            'exon:chr12:42792656-42792796:+']


@pytest.fixture()
def true_exc_nmd():
    return True


@pytest.fixture()
def true_inc_nmd():
    return True


@pytest.fixture()
def strand_true_exc_nmd():
    return '+'


@pytest.fixture()
def strand_true_inc_nmd():
    return '+'


@pytest.fixture()
def true_dict():
    return {'ENST00000443963.1': ['start_codon:chr3:186501400-186501402:+',
                                  'CDS:chr3:186501400-186501428:+:0',
                                  'exon:chr3:186501386-186501428:+',
                                  'exon:chr3:186502218-186502266:+',
                                  'exon:chr3:186502353-186502485:+',
                                  'exon:chr3:186503672-186503840:+',
                                  'exon:chr3:186503953-186504062:+',
                                  'exon:chr3:186504291-186504434:+',
                                  'exon:chr3:186504916-186505053:+',
                                  'exon:chr3:186505284-186505373:+',
                                  'exon:chr3:186505592-186505671:+',
                                  'exon:chr3:186506914-186507686:+'],
            'ENST00000441007.1': ['start_codon:chr3:186501400-186501402:+',
                                  'CDS:chr3:186501400-186501428:+:0',
                                  'exon:chr3:186500994-186501139:+',
                                  'exon:chr3:186501237-186501428:+',
                                  'exon:chr3:186502221-186502266:+',
                                  'exon:chr3:186502353-186502448:+'],
            'ENST00000498746.1': ['start_codon:chr3:186502243-186502245:+',
                                  'CDS:chr3:186502243-186502266:+:0',
                                  'exon:chr3:186501992-186502266:+',
                                  'exon:chr3:186502353-186502485:+',
                                  'exon:chr3:186502751-186502890:+',
                                  'exon:chr3:186503672-186503702:+'],
            'ENST00000429589.1': ['start_codon:chr3:186501400-186501402:+',
                                  'CDS:chr3:186501400-186501428:+:0',
                                  'exon:chr3:186501366-186501428:+',
                                  'exon:chr3:186502221-186502266:+',
                                  'exon:chr3:186502353-186502485:+',
                                  'exon:chr3:186503672-186503840:+',
                                  'exon:chr3:186503953-186504062:+',
                                  'exon:chr3:186504291-186504434:+',
                                  'exon:chr3:186504916-186505053:+',
                                  'exon:chr3:186505268-186505373:+',
                                  'exon:chr3:186505592-186505671:+',
                                  'exon:chr3:186506914-186506929:+'],
            'ENST00000323963.5': ['start_codon:chr3:186501400-186501402:+',
                                  'CDS:chr3:186501400-186501428:+:0',
                                  'exon:chr3:186501336-186501428:+',
                                  'exon:chr3:186502221-186502266:+',
                                  'exon:chr3:186502353-186502485:+',
                                  'exon:chr3:186502751-186502890:+',
                                  'exon:chr3:186503672-186503840:+',
                                  'exon:chr3:186503953-186504062:+',
                                  'exon:chr3:186504291-186504434:+',
                                  'exon:chr3:186504916-186505053:+',
                                  'exon:chr3:186505284-186505373:+',
                                  'exon:chr3:186505592-186505671:+',
                                  'exon:chr3:186506914-186507689:+'],
            'ENST00000445596.1': ['start_codon:chr3:186501400-186501402:+',
                                  'CDS:chr3:186501400-186501428:+:0',
                                  'exon:chr3:186501094-186501428:+',
                                  'exon:chr3:186502221-186502266:+',
                                  'exon:chr3:186502353-186502485:+',
                                  'exon:chr3:186502751-186502836:+'],
            'ENST00000425053.1': ['start_codon:chr3:186501400-186501402:+',
                                  'CDS:chr3:186501400-186501428:+:0',
                                  'exon:chr3:186501366-186501428:+',
                                  'exon:chr3:186502221-186502266:+',
                                  'exon:chr3:186502353-186502485:+',
                                  'exon:chr3:186502751-186502890:+',
                                  'exon:chr3:186503672-186503840:+',
                                  'exon:chr3:186503953-186504062:+',
                                  'exon:chr3:186504291-186504434:+',
                                  'exon:chr3:186504916-186505053:+',
                                  'exon:chr3:186505284-186505373:+',
                                  'exon:chr3:186505592-186505671:+',
                                  'exon:chr3:186506099-186506205:+',
                                  'exon:chr3:186506914-186507670:+'],
            'ENST00000440191.2': ['start_codon:chr3:186501400-186501402:+',
                                  'CDS:chr3:186501400-186501428:+:0',
                                  'exon:chr3:186501366-186501428:+',
                                  'exon:chr3:186502218-186502266:+',
                                  'exon:chr3:186502353-186502485:+',
                                  'exon:chr3:186502751-186502890:+',
                                  'exon:chr3:186503672-186503840:+',
                                  'exon:chr3:186503953-186504062:+',
                                  'exon:chr3:186504291-186504434:+',
                                  'exon:chr3:186504916-186505053:+',
                                  'exon:chr3:186505284-186505373:+',
                                  'exon:chr3:186505592-186505671:+',
                                  'exon:chr3:186506914-186507686:+'],
            'ENST00000426808.1': ['start_codon:chr3:186501400-186501402:+',
                                  'CDS:chr3:186501400-186501428:+:0',
                                  'exon:chr3:186501361-186501428:+',
                                  'exon:chr3:186502221-186502266:+',
                                  'exon:chr3:186502353-186502485:+',
                                  'exon:chr3:186503672-186503840:+',
                                  'exon:chr3:186503953-186504062:+',
                                  'exon:chr3:186504291-186504434:+',
                                  'exon:chr3:186504916-186505053:+',
                                  'exon:chr3:186505284-186505373:+',
                                  'exon:chr3:186505592-186505671:+',
                                  'exon:chr3:186506914-186507686:+'],
            'ENST00000356531.5': ['start_codon:chr3:186502423-186502425:+',
                                  'CDS:chr3:186502423-186502485:+:0',
                                  'exon:chr3:186501386-186501428:+',
                                  'exon:chr3:186502218-186502266:+',
                                  'exon:chr3:186502353-186502485:+',
                                  'exon:chr3:186503672-186503840:+',
                                  'exon:chr3:186503953-186504062:+',
                                  'exon:chr3:186504291-186504434:+',
                                  'exon:chr3:186504916-186505053:+',
                                  'exon:chr3:186505284-186505373:+',
                                  'exon:chr3:186505592-186505671:+',
                                  'exon:chr3:186506914-186507683:+']}


class TestNMDExons(object):

    @pytest.fixture
    def nmd_exons(self, database, exon_ids):
        from nmd import NMDExons
        return NMDExons(database, exon_ids)

    @pytest.fixture
    def inc_trans_without_exon(self, nmd_exons, inc_trans_without_exon_ids):
        return [nmd_exons.db[x] for x in inc_trans_without_exon_ids]

    @pytest.fixture
    def inc_trans_with_exon(self, nmd_exons, inc_trans_with_exon_ids):
        return [nmd_exons.db[x] for x in inc_trans_with_exon_ids]

    @pytest.fixture
    def exc_trans_without_exon(self, nmd_exons, exc_trans_without_exon_ids):
        return [nmd_exons.db[x] for x in exc_trans_without_exon_ids]

    @pytest.fixture
    def exc_trans_with_exon(self, nmd_exons, exc_trans_with_exon_ids):
        return [nmd_exons.db[x] for x in exc_trans_with_exon_ids]

    def test___init__(self, nmd_exons, exon_ids):
        assert isinstance(nmd_exons.db, gffutils.FeatureDB)
        pdt.assert_equal(nmd_exons.exon_ids, exon_ids)

    def test_find_nmd_exons(self, nmd_exons, is_exon_nmd, exon_ids):
        test = nmd_exons.find_nmd_exons()
        true = pd.Series(is_exon_nmd, index=exon_ids)
        pdt.assert_series_equal(test, true)

    def test__is_this_exon_nmd(self, is_exon_nmd, exon_ids, nmd_exons):
        for exon_id, true in zip(exon_ids, is_exon_nmd):
            test = nmd_exons._is_this_exon_nmd(exon_id)
            assert test == true

    def test__get_transcripts_with_exon(self, database, single_exon_id,
                                        nmd_exons, parent_transcripts_of_exon):
        db = gffutils.FeatureDB(database)
        exon = db[single_exon_id]
        test = nmd_exons._get_transcripts_with_exon(exon)
        true = parent_transcripts_of_exon
        assert test == true

    def test__get_all_transcripts_overlapping_exon(self, database,
                                                   single_exon_id, nmd_exons,
                                                   all_transcripts_of_exon):
        db = gffutils.FeatureDB(database)
        exon = db[single_exon_id]
        test = nmd_exons._get_all_transcripts_overlapping_exon(exon)
        true = all_transcripts_of_exon
        assert test == true

    def test__create_dict(self, all_transcripts_of_exon, strand_true_exc_nmd,
                          nmd_exons, true_dict):
        test = nmd_exons._get_exons_from_transcripts(all_transcripts_of_exon,
                                                     strand_true_exc_nmd)
        test = dict((key, [v.id for v in values]) for key, values in
                    test.items())
        true = true_dict
        pdt.assert_dict_equal(test, true)

    def test__inclusion_nmd(self, inc_trans_without_exon, inc_trans_with_exon,
                            strand_true_inc_nmd, nmd_exons, true_inc_nmd):
        test = nmd_exons._inclusion_nmd(inc_trans_without_exon,
                                        inc_trans_with_exon,
                                        strand_true_inc_nmd)
        true = true_inc_nmd
        assert test == true

    def test__exclusion_nmd(self, exc_trans_without_exon, exc_trans_with_exon,
                            strand_true_exc_nmd, nmd_exons, true_exc_nmd):
        test = nmd_exons._exclusion_nmd(exc_trans_without_exon,
                                        exc_trans_with_exon,
                                        strand_true_exc_nmd)
        true = true_exc_nmd
        assert test == true

    def test_nmd_stop_codon_nmd(self, database, stop_codon_exon_ids,
                                nmd_exons):
        for stop_codon_exon_pair in stop_codon_exon_ids:
            db = gffutils.FeatureDB(database)
            stop = stop_codon_exon_pair[0]
            stop_codon = db[stop]
            exon = db[stop_codon_exon_pair[1]]
            true = stop_codon_exon_pair[2]
            test = nmd_exons.stop_codon_nmd(exon, stop_codon)
            assert test == true
