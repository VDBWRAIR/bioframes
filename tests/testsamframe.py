from bioframes import samframe
import mock
import pandas as pd
from pandas.util.testing import assert_series_equal, assert_frame_equal, assert_index_equal
from numpy.testing import assert_array_equal, assert_array_almost_equal
from operator import itemgetter, attrgetter as attr
from testbiopandas import mock_file
import unittest

#class TestClassicSam(unittest.TestCase):
#    def setUp(self):
#        self.samtext	=	'\n'.join(['read1\t1\tchr1\t1\t60	10M	=	1	1	TTTCGAATC	FFFFFFFFF	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq',
#                                  'read2	1	chr2	1	60	10M	=	1	1	CTTCGATC	AFFDDDDD	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq',
#                                  'read3	1	chr1	1	60	10M	=	1	1	CCGATCAA	FF@@@F@F	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq'])
#    def test_sam_to_df(self):
#        result = pc.samview_to_df(self.samtext)
#        self.assertEquals(result.columns.tolist(), self.columns)
#        self.assertEquals(result.ix[2]['QNAME'], 'read3')
#
#    def test_df_from_collection_attributes(self):
#        mocks = [mock.Mock() for i in range(5)]
#        [[mock_do() for i in range(index)] for index, mock_do in enumerate(mocks)]
#        columns = ['call_count', 'called']
#        expected = pd.DataFrame( [(i, bool(i)) for i in range(5)])
#        expected.columns = columns
#        result = bf.df_from_collection_attributes(columns, mocks)
#        assert_frame_equal(expected, result)
#

class TestSamframe(unittest.TestCase):
    def setUp(self):
        #TODO: fix this so data is accurate;  i.e.:
   # ''' assert sum(itemgetter('M', 'I', 'S', '=', 'X')) == len(seq) == len(quality), \ "cigar string M/I/S/=/X should sum to the length of the query sequence." '''


        self.samtext='\n'.join([
                         'read1	60	chr1	1	1	10M	=	1	1	TTTCGAATC	FFFFFFFFF	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq',
                         'read2	8 	chr2	1	1	3I3M4D	=	1	1	CTTCGATC	AFFDDDDD	NM:i:3	AS:i:2	XS:i:0	RG:Z:Sanger',
                         'read3	60	chr1	1	1	2D2M2I2M2=	=	1	1	CCGATCAA	FF@@@F@F	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq'])

        self.result = mock_file(samframe.load_sam, self.samtext)

    def test_cigar_scores(self):
        #TODO: with mock_open
        e_strings = pd.Series(['10M', '3I3M4D', '2D2M2I2M2='])
        e_m, e_i, e_d, e_total = pd.Series([10, 3, 4]), pd.Series([float('nan'), 3, 2]), pd.Series([float('nan'), 4, 2]), pd.Series([0, 7, 4])
        assert_series_equal(self.result.cigar_I, e_i)
        assert_series_equal(self.result.cigar_M, e_m)
        assert_series_equal(self.result.cigar_D, e_d)
        assert_series_equal(self.result.cigar_score, e_total)
        assert_series_equal(self.result.CIGAR, e_strings)

    def test_options(self):
       e_nm, e_as, e_xs, e_rg = map(pd.Series, [[3, 3, 3], [231, 2, 231], [0, 0, 0], ['MiSeq', 'Sanger', 'MiSeq']])
       df = self.result
       results = [df.NM, df.AS, df.XS, df.RG]
       map(assert_series_equal, [e_nm, e_as, e_xs, e_rg], results)

    def test_df_with_options(self):
        #df = self.result.set_index( ['QNAME', 'POS', 'RNAME'])
        sam_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"] + ['NM', 'AS', 'XS', 'RG']
        expected = pd.Series(['read3',	60,	'chr1',	1,1,	'2D2M2I2M2=',	'=',	1,	1,	'CCGATCAA',	'FF@@@F@F',	3, 231, 0, 'MiSeq']).values#, 4, 2, 2, 4)
        res = self.result.ix[2][sam_columns].values
        assert_array_equal(expected, res)
        #assert_series_equal(df.ix['read3', 1, 'chr1'][:len(expected)], expected)

    def test_flags(self):
        df = self.result
        middle_e = pd.Series([False, False, False, True], dtype=object).values
        outer_e = pd.Series([False, False, True, True], dtype=object).values
        names = ["template having multiple segments in sequencing",
        "each segment properly aligned according to the aligner",
        "segment unmapped",
        "next segment in the template unmapped"]
#        inner_actual = df[df['QNAME'] == 'read2'][names]
#        outer_actual1= df[df['QNAME'] == 'read1'][names]
#        outer_actual3= df[df['QNAME'] == 'read3'][names]
        flag8, flag60_1, flag60_2  = map(attr('values'),  map(itemgetter(names), itemgetter(1, 0, 2)(df.ix)))
        assert_array_equal(middle_e,flag8)
        assert_array_equal(outer_e, flag60_1)
        assert_array_equal(outer_e, flag60_2)
