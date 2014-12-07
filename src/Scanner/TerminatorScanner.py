__author__ = 'Kevin'
from Bio import SeqIO

class ResultObject:
    """
    An Object for storing the sequence and annotation of the feature
    """

    def __init__(self, sequence, annotation):
        self.annotation = annotation
        self.sequence = sequence

    def __str__(self):
        return "Sequence: "+str(self.sequence)+" Annotation: "+self.annotation

def getFeatureInfo(records, featureType):
    seqList = []
    for record in records:
        if len(record.seq) > 1500:
            for feature in record.features:
                if feature.type == featureType:
                    location = feature.location
                    try:
                        qualifier = seqList.append(feature.qualifiers['note'])
                    except(KeyError):
                        qualifier = "null"
                    sequence = seqList.append(record.seq[location.start:location.end])
                    seqList.append(ResultObject(sequence, qualifier))
    return seqList

records = SeqIO.parse("/Users/dominikburri/"
                      "PycharmProjects/Bioinformatik/PlasmidAnalysis/"
                      "files/vectors-100.gb", "genbank")

liste = getFeatureInfo(records, "terminator")

for entry in liste:
    print entry

