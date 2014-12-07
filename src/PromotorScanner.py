__author__ = 'dominikburri'

from Bio import SeqIO

def generateList(feature_type):
    list_of_occurences = []
    list_of_annotations = []
    for record in records:
        if len(record.seq) > 1500: # minimum for number of bases
            for feature in record.features:
                if feature.type == feature_type:
                    sequence_of_feature = record.seq[feature.location.start: feature.location.end]
                    list_of_occurences.append(sequence_of_feature)
                    list_of_annotations.append(feature.qualifiers['note'])
    return list_of_occurences, list_of_annotations

records = SeqIO.parse("/Users/dominikburri/PycharmProjects/"
                      "Bioinformatik/PlasmidAnalysis/files/vectors-100.gb", "genbank")

#def clustering(list):


feature_type = 'promoter'
lists = generateList(feature_type)
print("occurences of " + feature_type + ": " + str(len(lists[0])))
for index in range(len(lists[0])):
    print lists[0][index]
    print lists[1][index]
