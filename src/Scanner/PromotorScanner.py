__author__ = 'dominikburri'

from Bio import SeqIO
from bioservices import *

class ResultObject:
    """
    An Object for storing the sequence and annotation of the feature
    """

    def __init__(self, sequence, annotation):
        self.annotation = annotation
        self.sequence = sequence
    def __str__(self):
        return str(self.sequence) + "; " + self.annotation

def generateList(feature_type, annotation_type):
    """
    A generator
    :param feature_type:
    :param annotation_type:
    :return: a ResultObject with the desired sequence and annotation
    """
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
                    yield result


records = SeqIO.parse("/Users/dominikburri/PycharmProjects/"
                      "Bioinformatik/PlasmidAnalysis/files/vectors-100.gb", "genbank")

def clustering(list_of_sequences):
    m = MUSCLE(verbose=False)
    jobid = m.run(frmt="fasta", sequence=list_of_sequences, email="dominik.burri1@students.fhnw.ch")
    while m.getStatus(jobid) == u'RUNNING':
        print "Status: ", m.getStatus(jobid)

    result=m.getResult(jobid, "sequence")
    print "sequence:"
    print result
    result=m.getResult(jobid, "aln-fasta")
    print "aln-fasta:"
    print result
    result=m.getResult(jobid, "phylotree")
    print "phylotree:"
    print result
    result=m.getResult(jobid, "pim")
    print "pim:"
    print result
    result=m.getResult(jobid, "out")
    print "out:"
    print result

# make a list generator with the desired feature and its annotation
feature_type = 'promoter'
annotation_type = 'note' # key in dictionary
list_generator = generateList(feature_type, annotation_type)
occurences = 0
list_of_sequences = []
for resultObject in list_generator:
    print resultObject
    occurences += 1
    list_of_sequences.append(">test\n" + resultObject.sequence + "\n")
print("occurences of " + feature_type + ": " + str(occurences))

clustering(list_of_sequences)
