
import nmdc
import json
from requests import Response
from unittest.mock import MagicMock, patch,  create_autospec


class MockWS():
    def get_workspace_info(self, params):
        return [100100,
                'nmdc:gold_Gs0114663',
                'nmdc', '2021-10-02T16:54:24+0000',
                6, 'a', 'n', 'unlocked',
                    {'narrative_nice_name': 'a name',
                    'cell_count': '1',
                    'searchtags': 'narrative',
                    'is_temporary': 'false',
                    'narrative': '1'
                    }
                 ]

    def list_objects(self, params):
        obj = [5,
               '1781_100329_metagenome',
               'KBaseMetagenomes.AnnotatedMetagenomeAssembly-1.0',
               '2021-10-02T16:52:14+0000',
               1,
               'nmdc',
               99947,
               'nmdc:gold_Gs0114663',
               '19f6a28d8de453ff405d6e78d9e8271c',
               1942025,
               None]
        return [obj]


class TestNMDC:
    _wsinfo = [
        99947,
        'nmdc:gold_Gs0114663',
        'nmdc', '2021-10-02T16:54:24+0000',
        6, 'a', 'n', 'unlocked',
            {'narrative_nice_name': 'a name',
            'cell_count': '1',
            'searchtags': 'narrative',
            'is_temporary': 'false',
            'narrative': '1'
            }
        ]

    samples = json.load(open('test/samples.json'))
    study = json.load(open('test/study.json'))

    def _mock_get_staged(that):
        return {}

    @patch('nmdc.Workspace')
    @patch('nmdc.requests')
    def test_find_new_data(self, mrq, mockws):
        n = nmdc.NMDC("gold:Gs0114663")
        n.ws = mockws
        n.get_study_info = MagicMock(return_value=self.study)
        mockws.get_workspace_info.return_value = self._wsinfo
        n.initialize_narrative()
        n._get_samples = MagicMock(return_value = self.samples[0:1])
        mock = [{'name': '1781_100342_assembly_contigs.fna'},
                {'name': '1781_100342_functional_annotation.gff'}]

        stage_resp = create_autospec(Response)
        stage_resp.status_code = 200
        stage_resp.json.return_value = mock
        mrq.get.return_value = stage_resp
        n.find_new_data(dryrun=True)
        # n.link_objects(dryrun=True)

    def test_make_table(self):
        n = nmdc.NMDC("gold:Gs0114663")
        n._get_samples = MagicMock(return_value = self.samples[0:1])
        n.build_sample_map()
        n.make_table()

    @patch('nmdc.requests')
    def test_make_samples(self, mrq):
        n = nmdc.NMDC("gold:Gs0114663")
        n._get_samples = MagicMock(return_value = self.samples[0:2])
        n.build_sample_map()
        n.make_samples()
        mrq.post.assert_called_once()
        with open('/tmp/gold:Gs0114663.tsv') as f:
            data = f.read()
            print(data)
