# python3
# ################################################################## #
# EM_biopython_extras.py (C) June-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project: geneLength
#
# Purpose: module to deal with some extra functions for biopython
#
# ################################################################## #
import sys
#from Bio import SeqIO
import re

def parse_uniprot_header(record):
    """
        From a record (uniprot reference proteoma fasta) from seqIO.parse => parse the header

        Parameters
        ----------
        record: object obtained from seqIO.parse (biopython)
            It is an entry from a multifasta

        Output:
        ------
        db: str
            database:sp or tr
        UniprotID: str
        EntryName: str
        ProteinName: str
        OS: str
            OrganiSm name
        OX: str
            Organism Identifier
        GN: str
            Gene Name
        PE: str
            Protein Evidence
        SV: str
            Sequence Variant
    """
    id=record.id
    description=record.description
    if 0:
        print("id:(",id,")")
        print("description:(",description,")")

    db, UniprotID, EntryName, ProteinName = "None", "None", "None", "None"
    OS,OX,GN,PE,SV = "None", "None", "None", "None", "None"
    db=str(id.split('|')[0])
    UniprotID=str(id.split('|')[1])
    EntryName=str(id.split('|')[2])

    m=re.search(r'\s(.*)\sOS=',description)
    ProteinName=str(m.group(1))
    m=re.search(r'OS=(.*)\sOX=',description)
    if m:
        OS=str(m.group(1))
    m=re.search(r'OX=(\w+)\s',description)
    if m:
        OX=str(m.group(1))
    m=re.search(r'GN=(\w*)\s',description)
    if m:
        GN=str(m.group(1))
    m=re.search(r'PE=(\d+)\s',description)
    if m:
        PE=str(m.group(1))
    m=re.search(r'SV=(\d+)',description)
    if m:
        SV=str(m.group(1))

    if 0:
        print("\ndb:(%s)\tUniprotID:(%s)\tEntryName:(%s)\tProteinName:(%s)" % (db, UniprotID, EntryName, ProteinName))
        print("OS:(%s)\tOX:(%s)\tGN:(%s)\tPE:(%s)\tSV:(%s)"%(OS,OX,GN,PE,SV))
    return db,UniprotID,EntryName,ProteinName,OS,OX,GN,PE,SV
