__author__ = 'dominikburri'

from Bio import SeqIO

class ResultObject:
    """
    An Object for storing the sequence and annotation of the feature
    """
    sequence = ""
    annotation = []

    def __init__(self, sequence, annotation):
        self.annotation = annotation
        self.sequence = sequence

def generateList(feature_type):
    list_of_occurences = []
    for record in records:
        if len(record.seq) > 1500: # minimum for number of bases
            for feature in record.features:
                if feature.type == feature_type:
                    sequence_of_feature = record.seq[feature.location.start: feature.location.end]
                    annotation = feature.qualifiers['note']
                    list_of_occurences.append(ResultObject(sequence_of_feature, annotation))
    return list_of_occurences

records = SeqIO.parse("/Users/dominikburri/PycharmProjects/"
                      "Bioinformatik/PlasmidAnalysis/files/vectors-100.gb", "genbank")

#def clustering(list):


feature_type = 'promoter'
list = generateList(feature_type)
print("occurences of " + feature_type + ": " + str(len(list)))
for index in list:
    print index.sequence
    print index.annotation
