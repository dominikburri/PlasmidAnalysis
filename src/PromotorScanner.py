__author__ = 'dominikburri'

from Bio import SeqIO

class ResultObject:
    """
    An Object for storing the sequence and annotation of the feature
    """
    sequence = ""
    annotation = ""

    def __init__(self, sequence, annotation):
        self.annotation = annotation
        self.sequence = sequence
    def __str__(self):
        return str(self.sequence) + "; " + self.annotation

def generateList(feature_type, annotation_type):
    list_of_occurences = []
    for record in records:
        if len(record.seq) > 1500: # minimum for number of bases
            for feature in record.features:
                if feature.type == feature_type:
                    sequence_of_feature = record.seq[feature.location.start: feature.location.end]
                    try:
                        annotation = feature.qualifiers[annotation_type][0]
                    except KeyError:
                        annotation = "null"
                    result = ResultObject(sequence_of_feature, annotation)
                    list_of_occurences.append(result)
    return list_of_occurences

records = SeqIO.parse("/Users/dominikburri/PycharmProjects/"
                      "Bioinformatik/PlasmidAnalysis/files/vectors-100.gb", "genbank")

#def clustering(list):


feature_type = 'promoter'
annotation_type = 'note'
list = generateList(feature_type, annotation_type)
print("occurences of " + feature_type + ": " + str(len(list)))
for resultObject in list:
    print resultObject
