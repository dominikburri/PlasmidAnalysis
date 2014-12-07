__author__ = 'Kevin'
from Bio import SeqIO

def terminator(records):
    seqList = []
    for record in records:
        if len(record.seq) > 1500:
            for feature in record.features:
                if feature.type == "CDS":
                    location = feature.location
                    seqList.append(feature.qualifiers['note'])
                    seqList.append(record.seq[location.start:location.end])
    return seqList

records = SeqIO.parse("C:\Users\Kevin\IdeaProjects\PlasmidAnalysis\\files\\vectors-100.gb", "genbank")
liste = terminator(records)

for eintrag in liste:
    print eintrag