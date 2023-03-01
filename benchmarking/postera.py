import requests, os, json, gzip
import plotly.express as px
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List, Dict, Optional, Any, Union

class SinglePostera:
    """
    The Posera SDK by Rubén does a lot more,
    this is a simple one aimed at pandas.
    It does log the queries...
    
    """
    # https://api.postera.ai/api/v1/docs/#operation/api_v1_post
    def __init__(self, **options):
        self.options = options
        self.responses: List[Tuple[str, Dict[str, Union[str, float, int]], Any]] = []
        
    def dejavu(self, url: str, data: dict) -> Union[dict, None]:
        """
        """
        h = hash(tuple(data.items()))
        for u, d, r in self.responses:
            if u == url and h == hash(tuple(d.items())):
                return r
        else:
            None
            
    def load_cache(self, fn: str='responses.json'):
        extension = os.path.splitext(fn)[1]
        if extension == '.json':
            with open(fn) as fh:
                self.responses += json.load(fh)
        elif extension in ('.p', '.pkl'):
            with open(fn, 'rb') as fh:
                self.responses += pickle.load(fh)
        elif os.path.splitext(fn)[1] == '.gz':
            with gzip.open(fn, 'rb') as fh:
                self.responses += pickle.load(fh)
            
                
    def dump_cache(self, fn: str='responses.json'):
        extension = os.path.splitext(fn)[1]
        if extension == '.json':
            with open(fn, 'w') as fh:
                json.dump(self.responses, fh)
        elif extension in ('.p', '.pkl'):
            with open(fn, 'wb') as fh:
                pickle.dump(self.responses, fh)
        elif os.path.splitext(fn)[1] == '.gz':
            with gzip.open(fn, 'wb') as fh:
                pickle.dump(self.responses, fh)
            
            

    def _post(self, url: str, data: Dict[str, Union[str, float, int]]) -> dict:
        # check if done before
        past_response = self.dejavu(url, data)
        if past_response is not None:
            return past_response
        if 'smiles' not in data:
            raise ValueError('no smiles?')
        if not data['smiles'] or not isinstance(data['smiles'], str):
            return {}
        headers = {'X-API-KEY': os.environ["MANIFOLD_API_KEY"]}
        response: requests.Response = requests.post(url, headers=headers, json=data)
        if response.status_code == 429:
            raise requests.ConnectionError(response.text)
        elif not response.ok:  # error tolerant
            print(data, response.status_code, response.text)
            return {}
        response_json = response.json()
        results: dict = response_json['results'] if 'results' in response_json else response_json
        self.responses.append((url, data, results))
        return results

    def retrosynthesis(self, smiles) -> dict:
        return self._post('https://api.postera.ai/api/v1/retrosynthesis/',
                          {'smiles': smiles, **self.options})
    
    def sa_retro(self, smiles) -> dict:
        return self._post('https://api.postera.ai/api/v1/synthetic-accessibility/retrosynthesis/',
                          {'smiles': smiles, **self.options})

    def similarity(self, smiles) -> dict:
        return self._post('https://api.postera.ai/api/v1/similarity/', 
                          {'smiles': smiles, **self.options})

    def fast(self, smiles) -> dict:
        return self._post('https://api.postera.ai/api/v1/synthetic-accessibility/fast-score/', 
                          {'smiles': smiles, **self.options})
    
def extract_vendors(data: List[Dict[str, Any]], flat=False) \
                        -> Union[Dict[str, Dict[str, Union[str, float]]], Dict[str, Any]]:
    # per compound
    vendored = {}

    for mol_info in data:
        for catalogue in mol_info['catalogEntries']:
            if catalogue['catalogName'] in vendored and \
                mol_info['similarity'] < vendored[catalogue['catalogName']]['similarity']:
                continue
            vendored[catalogue['catalogName']] = dict(similarity=mol_info['similarity'], 
                                                      smiles=mol_info['smiles'],
                                                      name=catalogue['catalogId'])
    if not flat:
        return vendored
    return flatten_vendors(vendored)


def flatten_vendors(data: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    return {f'{vendor}::{k}': v for vendor in data for k, v in data[vendor].items()}