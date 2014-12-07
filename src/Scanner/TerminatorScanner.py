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
        return "Annotation: "+self.annotation+" Sequence: "+str(self.sequence)

def getFeatureInfo(records, featureType):
    for record in records:
        if len(record.seq) > 1500:
            for feature in record.features:
                if feature.type == featureType:
                    location = feature.location
                    try:
                        qualifier = feature.qualifiers['note'][0]
                    except(KeyError):
                        qualifier = "null"
                    sequence = record.seq[location.start:location.end]
                    yield ResultObject(sequence, qualifier)


records = SeqIO.parse("C:\Users\Kevin\IdeaProjects\PlasmidAnalysis\\files\\vectors-100.gb", "genbank")

liste = getFeatureInfo(records, "terminator")

for entry in liste:
    print entry

