import requests
import json
import pandas as pd
import os
import sys
from biokbase.workspace.client import Workspace
from yaml import load, Loader
try:
    from biokbase.narrative.jobs.appmanager import AppManager
    from IPython.display import Javascript
except:
    pass

def _install_client():
    _base = "https://raw.githubusercontent.com/kbaseapps/GenericsAPI/master/lib/installed_clients"
    _cdir = "/tmp/clients"
    os.makedirs(_cdir, exist_ok=True)
    for i in ["SampleServiceClient.py", "baseclient.py", "authclient.py", "__init__.py"]:
        fn = os.path.join(_cdir, i)
        if os.path.exists(fn):
            continue
        url = '%s/%s' % (_base, i)
        resp = requests.get(url)
        open(fn, "w").write(resp.text)
    sys.path.append(_cdir)

try:
    from SampleServiceClient import SampleService
except:
    # Install Sample Service Client
    _install_client()
    from SampleServiceClient import SampleService


class NMDC:
    """
    NMDC Helper Class
    """
    _allowed_file_types = ["Assembly Contigs", "Functional Annotation GFF"]
    _nmdc_types = ["nmdc:MetagenomeAssembly", "nmdc:MetagenomeAnnotation"]
    _nmdc_url = "https://data.microbiomedata.org/"
    _service_url = "https://kbase.us/services"
    _staging = "%s/staging_service" % (_service_url)

    _file_types = {
                   "Assembly Contigs": {
                        "dirname": "assembly",
                        "ext": "assembly_contigs.fna"
                      },
                   "Functional Annotation GFF": {
                        "dirname": "annotation",
                        "ext": "functional_annotation.gff"
                      }}

    def __init__(self, study_id):
        self.study_id = study_id
        self.wsid = int(os.environ["KB_NARRATIVE"].split('.')[1])
        self.ws = Workspace("%s/ws" % (self._service_url))
        self.ss = SampleService("%s/sampleservice" % (self._service_url))
        self.samples = self._find_samples()
        self.oid2sid = self._build_sample_map()
        self.headers = {"Authorization": os.environ["KB_AUTH_TOKEN"]}
        m = "https://raw.githubusercontent.com/microbiomedata/nmdc_kbase/main/mapping.yaml"
        resp = requests.get(m)
        self.mapping = load(resp.text, Loader=Loader)

    # Query NMDC for samples related to this study
    def _find_samples(self):
        # TODO: Deal with more than 100 hits
        url = "%sapi/biosample/search?offset=0&limit=100" % (self._nmdc_url)
        q = {"conditions": [
                {"value": self.study_id,
                 "table": "study",
                 "op": "==",
                 "field": "study_id"}
                ],
             "data_object_filter": []
             }
        h = {
             "content-type": "application/json",
             "accept": "application/json"
             }
        resp = requests.post(url, headers=h, data=json.dumps(q))

        if resp.status_code == 200:
            return resp.json()['results']
        return None

    def _build_sample_map(self):
        # Build a map from omic id to NMDC Sample ID
        oid2sid = {}
        for sample in self.samples:
            sid = sample['id']
            for op in sample['omics_processing']:
                if op['annotations']['omics_type'] in ["Metagenome"]:
                    for act in op['omics_data']:
                        if act['type'] != "nmdc:MetagenomeAssembly":
                            continue
                        oid = act['name'].split(' ')[-1]
                        oid2sid[oid] = sid
        return oid2sid

    # Generate Sample TSV and store in Staging
    def make_samples(self):
        mapping = self.mapping["Mapping"]
        headings = mapping.keys()
        data = []
        data.append("\t".join(headings))
        for sample in self.samples:
            has_meta = False
            for op in sample['omics_processing']:
                if op['annotations']['omics_type'] == "Metagenome":
                    has_meta = True
            if not has_meta:
                continue
            row = []
            # Workaround
            if sample['env_local_scale_id'] in ['ENVO:01000621']:
                sample['env_local_scale_id'] = ""
            for c in headings:
                row.append(str(sample[mapping[c]]))
            data.append('\t'.join(row))
        url = "%s/upload" % (self._staging)
        fn = os.path.join("/tmp", "%s.tsv" % (self.study_id))
        with open(fn, "w") as f:
            f.write('\n'.join(data))
        params = {'destPath': '/'}
        files = {"uploads": open(fn, "rb")}
        resp = requests.post(url, headers=self.headers, data=params, files=files)
        return resp.json()

    # Make a table of the samples
    # Just focus on the ones with Metagenomes for now.
    def make_table(self):
        data = []
        cols = ["id", "name", "latitude", "longitude", "env_broad_scale_id",
                "env_local_scale_id", "env_medium_id", "description",
                "study_id", "depth", "collection_date", "ecosystem",
                "ecosystem_category", "ecosystem_type", "ecosystem_subtype",
                "specific_ecosystem"]

        for sample in self.samples:
            row = []
            has_meta = False
            for op in sample['omics_processing']:
                if op['annotations']['omics_type'] == "Metagenome":
                    has_meta = True
            if not has_meta:
                continue
            for col in cols:
                row.append(sample[col])
            data.append(row)
        return pd.DataFrame(data, columns=cols)

    def _get_oid(self, op):
        for act in op['omics_data']:
            if act['type'] == "nmdc:MetagenomeAssembly":
                return act['name'].split(' ')[-1]


    # This attempts to construct URLs to the data that we would like
    def _proc_metag(self, op, check=False):
        return results

    # Try to build URLs
    def get_urls(self, check=False):
        """
        Get the list of URLs for valid types.
        """
        headers = {"Range": "bytes=0-100"}  # first 100 bytes
        results = []
        for sample in self.samples:
            for op in sample['omics_processing']:
                # We just handle metagenome processing for now
                if op['annotations']['omics_type'] != "Metagenome":
                    continue
            oid = self._get_oid(op)
            urlbase = "%sdata/%s" % (self._nmdc_url, oid)
            for act in op['omics_data']:
                if act['type'] not in self._nmdc_types:
                    continue
                for do in act['outputs']:
                    if do['file_type'] not in self._file_types:
                        continue
                    ft = self._file_types[do['file_type']]
                    url = "%s/%s/%s_%s" % (urlbase, ft['dirname'],
                                           oid, ft['ext'])
                    if check:
                        resp = requests.get(url, headers=headers)
                        if resp.status_code != 206:
                            sys.stderr.write("Bad link %s\n" % (url))
                            continue
                    results.append(url)
        return results

    def build_staging_url(self):
        """
        return back list of urls for the staging app
        """
        url = "%s/list" % (self._staging)
        resp = requests.get(url, headers=self.headers)
        got_it = {}
        for fo in resp.json():
            got_it[fo["name"]] = fo

        results = []
        for url in self.get_urls():
            fn = url.split('/')[-1]
            if fn in got_it:
                continue
            results.append({"file_url": url})
        return results

    def add_staging(self):
        # Generate Staging App
        import json
        li = self.build_staging_url()
        code = """
                    Jupyter.narrative.addAndPopulateApp(
                                    "kb_uploadmethods/upload_web_file",
                                    'beta',
                                    {
                                    "download_type": "Direct Download",
                                    "urls_to_add_web_unpack": %s
                                    }
                                );
        """ % (json.dumps(li))
        #print(li)
        return Javascript(data=code, lib=None, css=None)

    def build_batch_import(self, max=None):
        # Build Batch Import
        # TODO: Make it skip existing items
        amat = 'KBaseMetagenomes.AnnotatedMetagenomeAssembly'
        done = {}
        for obj in self.ws.list_objects({'ids': [self.wsid], "type": amat}):
            done[obj[1]] = 1
        ct = 0
        batch_params = []
        ext = "_functional_annotation.gff"
        for link in self.get_urls():
            if link.endswith(ext):
                fn = link.split('/')[-1]
                ct += 1
                id = fn.replace(ext, '')
                oname = "%s_metagenome" % (id)
                if oname in done:
                    continue
                p = {
                    "fasta_file": "%s_assembly_contigs.fna" % (id),
                    "gff_file": "%s_functional_annotation.gff" % (id),
                    "genome_name": oname,
                    "source": "Other",
                    "release": "",
                    "genetic_code": None,
                    "generate_missing_genes": 1
                }
                batch_params.append(p)
                if max and ct >= max:
                    break
        return batch_params

    def make_bulk_import(self):
        return AppManager().run_app_bulk(
            [{
                "app_id": "kb_uploadmethods/import_gff_fasta_as_metagenome_from_staging",
                "tag": "release",
                "params": self.build_batch_import()
            }],
            cell_id="e3d8288a-9191-48c3-b39d-271306b0eda1",
            run_id="c4198db6-1566-4aed-83e7-92b61069affb"
        )

    # Link Objects
    def link_objects(self, dryrun=False):
        # Find the sample set
        # Note: This assume only one sample set in the workspace
        sset = self.ws.list_objects({'ids': [self.wsid],
                                     'type': 'KBaseSets.SampleSet'})[0]
        ssid = "%d/%d/%d" % (sset[6], sset[0], sset[4])

        # Build a map of sample names to sample IDs in KBase
        sample2obj = {}
        res = self.ws.get_objects2({'objects': [{'ref': ssid}]})
        for s in res['data'][0]['data']['samples']:
            sample2obj[s['name']] = s
        types = ['KBaseMetagenomes.AnnotatedMetagenomeAssembly',
                 'KBaseGenomeAnnotations.Assembly']
        for obj in self.ws.list_objects({'ids': [self.wsid]}):
            otyp = obj[2].split('-')[0]
            if otyp not in types:
                continue
            name = obj[1]
            upa = "%d/%d/%d" % (obj[6], obj[0], obj[4])
            oid = name.rstrip("_metagenome").rstrip("_metagenome.assembly")
            sample_name = self.oid2sid[oid]
            sampleid = sample2obj[sample_name]['id']
            print("Linking %s to %s" % (upa, sampleid))
            sample = self.ss.get_sample({'id': sampleid})
            if dryrun:
                continue
            self.ss.create_data_link({
                'upa': upa,
                'id': sampleid,
                'version': sample2obj[sample_name]['version'],
                'node': sample['node_tree'][0]['id'],
                'update': 1
            })


if __name__ == "__main__":
    n = NMDC("gold:Gs0114663")
    #print(n.make_samples())
    for l in n.get_urls():
        print(l)
    print(n.build_staging_url())
    sys.exit()
    print(json.dumps(n.build_batch_import(max=5), indent=2))
    print(n.make_table())
    n.link_objects(dryrun=True)
