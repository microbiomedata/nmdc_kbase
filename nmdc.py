import requests
import json
import pandas as pd
import os
import sys
from installed_clients.workspaceClient import Workspace
from installed_clients.execution_engine2Client import execution_engine2
from installed_clients.SampleServiceClient import SampleService
from transfer import transfer
from yaml import load, Loader
from uuid import uuid4
# try:
#     from biokbase.narrative.jobs.appmanager import AppManager
#     from IPython.display import Javascript
# except:
#     pass

_DEBUG = os.environ.get("DEBUG")


def _debug(text):
    if _DEBUG:
        print(text)


class File:
    """
    File Object: This captures important information about an
    NMDC file.

    Variables:
    fn: File name
    src: Source path (relative to base)
    url: data url
    """
    _nmdc_url = "https://data.microbiomedata.org/"

    def __init__(self, oid, ft):
        ext = ft['ext']
        dirname = ft['dirname']
        fn = "%s_%s" % (oid, ext)
        self.fn = fn
        self.src = '%s/%s/%s' % (oid, dirname, fn)
        self.url = "%sdata/%s/%s/%s_%s" % (self._nmdc_url,
                                           oid,
                                           dirname,
                                           oid,
                                           ext)


class Sample:
    """
    Sample Object: This holds information about the sample.

    Variables
    oid: the omic activity id
    sid: The same name
    name: The object name in KBase
    staged: T/F if all data is staged
    ftypes: The File Types for the data type
    params: The params block for the import
    file_list: list of files to stage (strings)
    files: List of File objects to stage.
    """

    _data = open("mapping.yaml").read()
    _mapping = load(_data, Loader=Loader)

    def __init__(self, oid, typ, sid=None, sinfo=None):
        self.oid = oid
        self.type = typ
        self.name = "%s_%s" % (oid, typ)
        self.sid = sid
        self.staged = False
        ftypes = self._mapping["Files"][typ]
        imprt = self._mapping["Imports"][typ]
        self.params = {
            imprt['name']: self.name
        }
        self.file_list = []
        self.files = []
        self.sample_info = {}
        for ft in ftypes:
            fo = File(oid, ftypes[ft])
            self.files.append(fo)
            self.file_list.append(fo.fn)
            self.params[ftypes[ft]['param']] = fo.fn

        if sinfo:
            for k in self._mapping['Mapping']:
                col = self._mapping['Mapping'][k]
                self.sample_info[col] = sinfo[col]


class NMDC:
    """
    NMDC Helper Class
    """
    _nmdc_url = "https://data.microbiomedata.org/"
    _service_url = "https://kbase.us/services"
    _staging = "%s/staging_service" % (_service_url)

    def __init__(self, study_id):
        """
        study_id: The NMDC study_id (e.g. gold:Gs123456)
        """
        self.study_id = study_id
        self.wsid = None
        self.ws = Workspace("%s/ws" % (self._service_url))
        self.ss = SampleService("%s/sampleservice" % (self._service_url))
        self.ee = execution_engine2("%s/ee2-nr" % (self._service_url))
        self.study = None
        # self.get_study_info()
        self.o_samples = None
        self.headers = {"Authorization": os.environ["KB_AUTH_TOKEN"]}
        if os.path.exists("mapping.yaml"):
            data = open("mapping.yaml").read()
        else:
            m = "https://raw.githubusercontent.com/" + \
                "microbiomedata/nmdc_kbase/main/mapping.yaml"
            data = requests.get(m).text
        self.mapping = load(data, Loader=Loader)

    # Query NMDC for samples related to this study
    def _get_samples(self):
        """
        Queries NMDC API to get the list of samples for the study.
        """

        # TODO: Deal with more than 100 hits
        url = "%sapi/biosample/search?offset=0&limit=100" % (self._nmdc_url)
        q = {"conditions": [
                {"value": self.study_id,
                 "table": "study",
                 "op": "==",
                 "field": "study_id"
                 }],
             "data_object_filter": []
             }
        h = {
             "content-type": "application/json",
             "accept": "application/json"
             }
        resp = requests.post(url, headers=h, data=json.dumps(q))
        if resp.status_code == 200:
            samples = resp.json()['results']
            return samples
        return None

    def build_sample_map(self):
        """
        Finds the relevant samples and creates Sample objects.
        Creates a map from the omic activity ID to the sample object.
        """

        # Right now we extract the oid from the name
        # of the metagnome assembly.
        self.oid2sample = {}
        samples = []
        for sample in self._get_samples():
            sid = sample['id']
            oid = self._get_oid(sample)
            # oid = None
            # for op in sample['omics_processing']:
            #     # oid = self._get_oid(op)
            #     if op['annotations']['omics_type'] in ["Metagenome"]:
            #         for act in op['omics_data']:
            #             if act['type'] != "nmdc:MetagenomeAssembly":
            #                 continue
            #             oid = act['name'].split(' ')[-1]
            if oid:
                sample = Sample(oid, 'metagenome', sid, sample)
                samples.append(sample)
                self.oid2sample[oid] = sample
        self.o_samples = samples

    def make_samples(self):
        """
        Create the Sample TSV file for uploading to KBase
        """

        mapping = self.mapping["Mapping"]
        headings = mapping.keys()
        data = []
        data.append("\t".join(headings))
        for sample in self.o_samples:
            row = []
            # Workaround
            sinfo = sample.sample_info
            if sinfo['env_local_scale_id'] in ['ENVO:01000621']:
                sinfo['env_local_scale_id'] = ""
            for c in headings:
                row.append(str(sinfo[mapping[c]]))
            data.append('\t'.join(row))
        url = "%s/upload" % (self._staging)
        fn = os.path.join("/tmp", "%s.tsv" % (self.study_id))
        with open(fn, "w") as f:
            f.write('\n'.join(data))
        params = {'destPath': '/'}
        files = {"uploads": open(fn, "rb")}
        resp = requests.post(url, headers=self.headers,
                             data=params, files=files)
        return resp.json()

    def submit_sample_import(self):
        """
        Incomplete: Submit the sample import for execution
        """

        fn = "%s.tsv" % (self.study_id)
        params = {
            "sample_file": fn,
            "file_format": "kbase",
            "set_name": "%s" % (self.study_id.replace(":", "_")),
            "header_row_index": None,
            "name_field": "",
            "description": self.study["name"],
            "sample_set_ref": None,
            "output_format": "",
            "taxonomy_source": "n/a",
            "num_otus": 20,
            "incl_seq": 0,
            "otu_prefix": "OTU",
            "incl_input_in_output": 1,
            "share_within_workspace": 1,
            "prevalidate": 1,
            "keep_existing_samples": 1,
            "ignore_warnings": 0
        }

        rjs = {
            "method": "sample_uploader/import_samples",
            "params": params,
            "wsid": self.wsid
        }
        return rjs

    def make_table(self):
        """
        Create a pandas data table of the samples
        """

        data = []
        cols = self.mapping['Table']['Columns']

        for sample in self.o_samples:
            row = []
            for col in cols:
                row.append(sample.sample_info[col])
            data.append(row)
        return pd.DataFrame(data, columns=cols)

    def _get_oid(self, sample):
        """
        Find the Omic activity for the sample.
        This is a bit of a hack at the moment.
        """

        oid = None
        for op in sample['omics_processing']:
            # oid = self._get_oid(op)
            if op['annotations']['omics_type'] in ["Metagenome"]:
                for act in op['omics_data']:
                    if act['type'] != "nmdc:MetagenomeAssembly":
                        continue
                    oid = act['name'].split(' ')[-1]
        return oid

    # Try to build URLs
    def get_urls(self):
        """
        Get the list of URLs for valid types.
        Deprecated
        """

        results = []
        for sample in self.o_samples:
            for f in sample.files:
                results.append(f.url)
        return results

    def get_staged_files(self):
        """
        Returns dictionary of staged files
        """

        url = "%s/list" % (self._staging)
        resp = requests.get(url, headers=self.headers)
        if resp.status_code != 200:
            sys.write.stderr(resp.text + '\n')
        staged = {}
        for fo in resp.json():
            staged[fo["name"]] = fo
        return staged

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
        for sample in self.samples:
            for f in sample.files:
                if f.fn in got_it:
                    continue
                results.append({"file_url": f.url})
        # for url in self.get_urls():
        #     fn = url.split('/')[-1]
        #     if fn in got_it:
        #         continue
        #     results.append({"file_url": url})
        return results

    # def add_staging(self):
    #     """
    #     Depcreated
    #     """

    #     # Generate Staging App
    #     import json
    #     li = self.build_staging_url()
    #     code = """
    #             Jupyter.narrative.addAndPopulateApp(
    #                 "kb_uploadmethods/upload_web_file",
    #                 'beta',
    #                 {
    #                 "download_type": "Direct Download",
    #                 "urls_to_add_web_unpack": %s
    #                 }
    #             );
    #     """ % (json.dumps(li))
    #     return Javascript(data=code, lib=None, css=None)

    def find_new_data(self, dryrun=False):
        """
        Query what is in the workspace and compare
        that against the samples assoicated with the
        study.  Come up with what needs to be staged
        and imported
        """
        if not self.o_samples:
            self.build_sample_map()

        # Let's start with what do have
        amat = 'KBaseMetagenomes.AnnotatedMetagenomeAssembly'
        done = {}
        for obj in self.ws.list_objects({'ids': [self.wsid], "type": amat}):
            done[obj[1]] = 1

        # Now let's go through the samples and see what's missing
        missing = []
        for sample in self.o_samples:
            if sample.name not in done:
                missing.append(sample)

        staged = self.get_staged_files()
        to_stage = []
        to_import = []
        for item in missing:
            ready = True
            for fo in item.files:
                if fo.fn not in staged:
                    ready = False
                    print("Staging %s" % (fo.src))
                    to_stage.append(fo.src)

            if ready:
                item.staged = True
                print("Import %s" % (item.oid))
                if not dryrun:
                    to_import.append(item)
        # TODO: Some way to see what is in flight
        if len(to_stage) > 0:
            transfer(to_stage)
        if not dryrun:
            self.submit_import(to_import)

    def submit_import(self, to_import, dryrun=False):
        """
        Create a bulk import job and append an app to the narrative.
        """
        cell_id = str(uuid4())
        run_id = str(uuid4())

        function_name = "kb_uploadmethods.upload_metagenome_fasta_gff_file"
        app_id = function_name.replace(".", "/")
        job_meta = {"cell_id": cell_id, "run_id": run_id}

        # plist is the parameter list for the batch import
        plist = []
        for sam in to_import:
            print("Import: %s" % (sam.oid))
            imprt = self.mapping['Imports'][sam.type]
            iv = {"workspace_name": self.workspace_name}
            for k, v in imprt['default_params'].items():
                iv[k] = v

            for k, v in sam.params.items():
                iv[k] = v
            param = {
                "method": function_name,
                "params": [iv],
                "app_id": app_id,
                "meta": job_meta,
            }
            plist.append(param)

        # Batch Param.. Just the workspace id
        bp = {'wsid': self.wsid}
        if dryrun:
            return
        resp = self.ee.run_job_batch(plist, bp)
        job_id = resp['batch_id']
        _debug(cell_id)
        _debug(job_id)
        if not dryrun:
            self.add_batch_cell(to_import, cell_id, run_id, job_id)

    def add_batch_cell(self, to_import, cell_id, run_id, job_id):
        """
        Append a batch input cell to the narrative.
        This requires:
        - the list of samples to import
        - The cell ID
        - The run ID (not used yet)
        - The job ID for the parent job.
        """
        cell = json.load(open('bulk_import.json'))
        # TODO run_id
        cell['metadata']['kbase']['attributes']['id'] = cell_id

        flist = []
        fpaths = []
        for item in to_import:
            # build params
            flist.extend(item.file_list)
            fpaths.append(item.params)
        blk = cell['metadata']['kbase']['bulkImportCell']
        blk['inputs']['gff_metagenome']['files'] = flist
        blk['params']['gff_metagenome']['filePaths'] = fpaths
        blk['exec']['jobs'] = {'byId': {}}
        js = self.ee.check_job({'job_id': job_id})
        job_list = js['child_jobs']
        job_list.append(job_id)
        jobs = self.ee.check_jobs({'job_ids': job_list})["job_states"]
        job_by_id = dict()
        # TODO: Can we just use the output from check_jobs?
        for job in jobs:
            job_by_id[job['job_id']] = job
        for job in job_list:
            jr = job_by_id[job]
            rec = {}
            rec['cell_id'] = cell_id
            for f in ['batch_id', 'batch_job', 'created',
                      'child_jobs', 'job_id',
                      'queued', 'retry_count', 'retry_ids',
                      'retry_saved_toggle',
                      'status', 'updated', 'user', 'wsid']:
                if f in jr:
                    rec[f] = jr[f]
            rec['job_output'] = jr.get('job_output', {})
            blk['exec']['jobs']['byId'][job] = rec
            if job == job_id:
                jstate = rec
        blk['exec']['jobState'] = jstate
        ref = self.narrative_ref
        resp = self.ws.get_objects2({"objects": [{"ref": ref}]})
        narr = resp['data'][0]
        usermeta = narr['info'][10]
        narrdata = narr['data']
        narrdata['cells'].append(cell)
        obj = {
            "objid": self.narrative_id,
            "type": "KBaseNarrative.Narrative-4.0",
            "data": narrdata,
            "meta": usermeta
            }
        if _DEBUG:
            json.dump(cell, open('debug.json', 'w'), indent=2)
        resp = self.ws.save_objects({"id": self.wsid, "objects": [obj]})

    def append_to_log(self, text, narrdata=None, save=False):
        """
        Append to text the Log cell of the narrative.
        WIP.
        """

        if not narrdata:
            # Read the narrative
            ref = self.narrative_ref
            resp = self.ws.get_objects2({"objects": [{"ref": ref}]})
            narr = resp['data'][0]
            usermeta = narr['info'][10]
            narrdata = narr['data']
        modified = False
        for cell in narrdata['cells']:
            if cell['source'].startswith("# Log"):
                cell['source'] += '%s\n' % (text)
                modified = True

        if modified and save:
            obj = {
                "objid": self.narrative_id,
                "type": "KBaseNarrative.Narrative-4.0",
                "data": narrdata,
                "meta": usermeta
                }
            resp = self.ws.save_objects({"id": self.wsid, "objects": [obj]})
        return narrdata

    # def build_batch_import(self, max=None):
    #     # Build Batch Import
    #     # May be deprecated soon
    #     amat = 'KBaseMetagenomes.AnnotatedMetagenomeAssembly'
    #     done = {}
    #     for obj in self.ws.list_objects({'ids': [self.wsid], "type": amat}):
    #         done[obj[1]] = 1
    #     ct = 0
    #     batch_params = []
    #     ext = "_functional_annotation.gff"
    #     for link in self.get_urls():
    #         if link.endswith(ext):
    #             fn = link.split('/')[-1]
    #             ct += 1
    #             id = fn.replace(ext, '')
    #             oname = "%s_metagenome" % (id)
    #             if oname in done:
    #                 continue
    #             p = {
    #                 "fasta_file": "%s_assembly_contigs.fna" % (id),
    #                 "gff_file": "%s_functional_annotation.gff" % (id),
    #                 "genome_name": oname,
    #                 "source": "Other",
    #                 "release": "",
    #                 "genetic_code": None,
    #                 "generate_missing_genes": 1
    #             }
    #             batch_params.append(p)
    #             if max and ct >= max:
    #                 break
    #     return batch_params

    # def make_bulk_import(self):
    #     """
    #     Create a bulk import cell (deprecated)
    #     """
    #     apid = "kb_uploadmethods/import_gff_fasta_as_metagenome_from_staging"
    #     return AppManager().run_app_bulk(
    #         [{
    #             "app_id": apid,
    #             "tag": "release",
    #             "params": self.build_batch_import()
    #         }],
    #         cell_id="e3d8288a-9191-48c3-b39d-271306b0eda1",
    #         run_id="c4198db6-1566-4aed-83e7-92b61069affb"
    #     )

    def _get_sample_set_map(self):
        """
        Returns a map from the sample name to the KBase
        sample ID
        """
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
        return sample2obj

    # Link Objects
    def link_objects(self, dryrun=False):
        """
        Link the imported objects to their samples.
        """
        sample2obj = self._get_sample_set_map()
        linked_types = self.mapping['Linking']
        # types = ['KBaseMetagenomes.AnnotatedMetagenomeAssembly',
        #          'KBaseGenomeAnnotations.Assembly']
        for obj in self.ws.list_objects({'ids': [self.wsid]}):
            otyp = obj[2].split('-')[0]
            if otyp not in linked_types:
                continue
            upa = "%d/%d/%d" % (obj[6], obj[0], obj[4])
            # We need to extract the omic ID from the name
            ext = linked_types[otyp]['ext']
            name = obj[1]
            # oid = name.rstrip("_metagenome").rstrip("_metagenome.assembly")
            oid = name.rstrip(ext)
            # Now we need to get the KBase Sample ID
            sample_name = self.oid2sample[oid].sid
            kbase_sid = sample2obj[sample_name]['id']
            print("Linking %s to %s" % (upa, kbase_sid))
            # Get info about the sample from KBase
            kbase_sample = self.ss.get_sample({'id': kbase_sid})
            if dryrun:
                continue
            self.ss.create_data_link({
                'upa': upa,
                'id': kbase_sid,
                'version': sample2obj[sample_name]['version'],
                'node': kbase_sample['node_tree'][0]['id'],
                'update': 1
            })

    def get_study_info(self):
        """
        Fetch the study information from NMDC
        """
        url = "%sapi/study/%s" % (self._nmdc_url, self.study_id)
        resp = requests.get(url)
        self.study = resp.json()

    def generate_markdown_header(self):
        """
        Generate the header markdown block for the narrative.
        """
        data = self.study
        publ = data['publication_doi_info']
        pubs = ""
        dataset = ""
        for p in publ:
            if 'link' in publ[p]:
                pubs += "[{title}]({link[0][URL]})\n".format(**publ[p])
            elif publ[p].get("type") == "dataset":
                dataset += "[{title}]({URL})\n".format(**publ[p])

        data['pubs'] = pubs
        data['dataset'] = dataset
        _TEMPLATE = open("template.md").read()

        text = _TEMPLATE.format(**data)

        return text

    def initialize_narrative(self):
        """
        Initialize the narrative.  If the study hasn't yet been
        initialized this will create the workspace and initalize
        the narrative with a markdown header.

        If the workspace already exists, this will populate the
        workspace name, workspace ID, Narrative ID, and narrative
        ref.
        """

        if not self.study:
            self.get_study_info()
        ws_name = "nmdc:%s" % (self.study_id.replace(":", "_"))
        try:
            resp = self.ws.get_workspace_info({"workspace": ws_name})
            self.wsid = resp[0]
            self.workspace_name = resp[1]
            self.narrative_id = resp[8]['narrative']
            self.narrative_ref = "%s/%s" % (self.wsid, self.narrative_id)
            print("Previously Initialized %d" % (self.wsid))
            return
        except:
            print("Add")
        markdown = self.generate_markdown_header()
        resp = self.ws.create_workspace({"workspace": ws_name})
        self.wsid = resp[0]
        self.workspace_name = resp[1]
        meta = {
            "narrative": "1",
            "narrative_nice_name": "%s" % (self.study['name']),
            "cell_count": 1,
            "searchtags": "narrative",
            "is_temporary": "false"
        }
        req = {"wsi": {"id": self.wsid}, "new": meta}
        self.ws.alter_workspace_metadata(req)
        narrative = json.load(open('narrative.json'))
        c0 = narrative['cells'][0]
        c0["source"] = markdown
        c0["metadata"]["kbase"]["attributes"]["title"] = self.study['name']
        narrative["metadata"]["ws_name"] = ws_name
        narrative["metadata"]["name"] = self.study['name']
        narrative["metadata"]["kbase"]["ws_name"] = ws_name
        usermeta = {
            "creator": "nmdc",
            "data_dependencies": "[]",
            "jupyter.markdown": "1",
            "is_temporary": "false",
            "job_info": "{\"queue_time\": 0, " +
                        "\"run_time\": 0, \"running\": 0, " +
                        "\"completed\": 0, \"error\": 0}",
            "format": "ipynb",
            "name": "%s" % (self.study['name']),
            "description": "",
            "type": "KBaseNarrative.Narrative",
            "ws_name": ws_name
        }
        obj = {
            "name": "narrative",
            "type": "KBaseNarrative.Narrative-4.0",
            "data": narrative,
            "meta": usermeta
            }
        resp = self.ws.save_objects({"id": self.wsid, "objects": [obj]})
        return resp


if __name__ == "__main__":
    n = NMDC("gold:Gs0114663")
    n.initialize_narrative()
    # n.append_to_log("This is a test")
    # n.find_new_data()
    n.build_sample_map()
    # n.link_objects()
    # n.get_urls()
    # print(n.make_table())
    sams = []
    for i in ["1781_100347", "1781_86097",
              "1781_86098", "1781_100341",
              "1781_100327"]:
        sams.append(n.oid2sample[i])
    n.submit_import(sams, dryrun=True)
    # cell_id = "15eaa207-c957-43a9-aa25-5ac8a6a0950f"
    # job_id = "615ba2f124131f13b740b76e"
    # n.add_batch_cell(sams, cell_id, None, job_id)

    # json.dump(data, open('samples.json', 'w'), indent=2)
    # n.initialize_narrative()
    # n.get_samples()
    # print(len(n.samples))
    # print(len(n.get_urls()))
    # print(n.submit_sample_import())
    # n.find_new_data()
