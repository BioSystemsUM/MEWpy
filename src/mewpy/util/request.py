# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
##############################################################################
Authors: Vitor Pereira
##############################################################################
"""
import json
import urllib.request
from hashlib import sha256


def process_entry(entry):

    protein = entry['primaryAccession']
    org = entry['organism']['scientificName']
    try:
        name = entry['genes'][0]['geneName']['value']
    except:
        tokens = entry['uniProtkbId'].split('_')
        name = tokens[0]
        print(f'No gene name for {protein} using uniProtkbId')
    props = {}
    props['Catalytic Activity'] = []
    # synonyms
    if 'synonyms' in entry['genes'][0].keys():
        l = entry['genes'][0]['synonyms']
        sysnames = ' '.join([a['value'] for a in l])
    else:
        sysnames = ''
    # ordered locus
    if 'orderedLocusNames' in entry['genes'][0].keys():
        l = entry['genes'][0]['orderedLocusNames']
        locus = [a['value'] for a in l]
    else:
        locus = []

    # comments
    try:
        for comment in entry['comments']:

            if comment['commentType'] == "CATALYTIC ACTIVITY":
                activity = comment['reaction']['name']
                ecnumber = ''
                try:
                    ecnumber = comment['reaction']['ecNumber']
                except Exception:
                    pass
                props['Catalytic Activity'].append((activity, ecnumber))

            elif comment['commentType'] == "ACTIVITY REGULATION":
                activity = comment['texts'][0]['value']
                props['Regulation Activity'] = activity

            elif comment['commentType'] == "PATHWAY":
                pathway = comment['texts'][0]['value']
                props['Pathway'] = pathway

            elif comment['commentType'] == "FUNCTION":
                function = comment['texts'][0]['value']
                props['Function'] = function
    except:
        print('No comments')

    # sequence
    seq = None
    mw = None
    try:
        seq = entry["sequence"]["value"]
    except:
        pass
    try:
        mw = float(entry["sequence"]["molWeight"])
    except:
        pass
    return {"organism": org,
            "protein": protein,
            "gene": name,
            'locus': locus,
            'synonyms': sysnames,
            'properties': props,
            'seq': seq,
            'mw': mw,
            }


def retreive(data, organism=None):
    """Retreives a protein function, pathways, ...

    Args:
        data ([type]): The protein data json

    Returns:
        A dictionary or a json string
    """
    d = json.loads(data)
    results = d['results']
    rid = 0
    if organism:
        for idx, entry in enumerate(results):
            n = entry['organism']['scientificName']
            if organism in n:
                rid = idx
                break
    entry = results[rid]
    return process_entry(entry)


def retreive_gene(gene, organism=None):
    """retrieves uniprot data using a gene name

    :param gene: gene name
    :type proteinid: str
    :param organism: the organism name
    :type organism: str, optional
    :return: the organism name, the protein assertion, the gene name, genes synonims, properties
    :rtype: tuple
    """
    gene = gene.strip()
    gn = gene.replace(' ', '+')
    url = f'https://rest.uniprot.org/uniprotkb/search?query=(gene:{gn})%20AND%20(reviewed:true)&format=json'
    with urllib.request.urlopen(url) as response:
        data = response.read().decode('ascii')
    data = retreive(data, organism)
    data['gene'] = gene
    return data


def retreive_protein(proteinid):
    """retrieves uniprot data using a protein assertion

    :param proteinid: protein assertion
    :type proteinid: str
    :return: the organism name, the protein assertion, the gene name, genes synonims, properties
    :rtype: tuple
    """
    url = f"https://rest.uniprot.org/uniprotkb/{proteinid}?format=json"
    with urllib.request.urlopen(url) as response:
        data = response.read().decode('ascii')
    entry = json.loads(data)
    return process_entry(entry)


def brenda_query(user, password, ecNumber, organism=None, field='KCAT'):
    org= "" if organism is None else organism
    wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    try:
        from zeep import Client
        client = Client(wsdl)
    except ImportError:
        raise Exception("zeep library is required.")
        
    passwd = sha256(password.encode("utf-8")).hexdigest()
    parameters = (user, passwd,
                  "ecNumber*"+ecNumber,
                  "organism*"+org
                  )

    if field == 'KCAT':
        resultString = client.service.getTurnoverNumber(*parameters)
    elif field == 'SEQ':
        resultString = client.service.getSequence(*parameters)
    elif field == 'MW':
        resultString = client.service.getMolecularWeight(*parameters)
    else:
        raise ValueError(f"{field} unknown field.")

    return resultString


def get_smiles(name):
    try:
        import requests
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/CanonicalSMILES/TXT" % name
        req = requests.get(url)
        if req.status_code != 200:
            smiles = None
        else:
            smiles = req.content.splitlines()[0].decode()
    except:
        smiles = None

    return smiles
