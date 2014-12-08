__author__ = 'j'
from Bio import SeqIO


vectors = SeqIO.parse('../../files/vectors-100.gb', 'genbank')

locationList = []
qualifierList = []

for seq_record in vectors:
    for feature in seq_record.features:
        if feature.type == 'rep_origin':
            locationList.append(feature.location)
            print feature.qualifiers.get("note")


#TODO: setup dictionary how often each term is used
#feature.qualifiers