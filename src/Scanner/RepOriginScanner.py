__author__ = 'j'

# from Bio import SeqIO
# from bioservices import *


#TODO: setup dictionary how often each term is used
#PSSM nur machen wenn die Annotation identisch ist


from Bio import SeqIO
from bioservices import *

class ResultObject:
    """
    An Object for storing the sequence and annotation of the feature
    """

    def __init__(self, sequence, feature_type, annotation):
        self.annotation = annotation
        self.feature_type = feature_type
        self.sequence = sequence
    def __str__(self):
        return str(self.sequence) + "; " + str(self.feature_type) + "; " + str(self.annotation)

def generateList(feature_type, type_count):
    """
    A generator
    :param feature_type:
    :return: a ResultObject with the desired sequence and annotation
    """
    list_key = []
    for record in records:
        if len(record.seq) > 1500: # minimum for number of bases
            for feature in record.features:

                for key in feature.qualifiers:
                    if key not in list_key:
                        list_key.append(key)


                if feature.type == feature_type:
                    type_count[feature.type] += 1
                    sequence_of_feature = record.seq[feature.location.start: feature.location.end]
                    annotation = feature.qualifiers
                    feature_type = feature.type
                    result = ResultObject(sequence_of_feature, feature_type, annotation)
                    yield result
    #print type_count
    return



def clustering(list_of_sequences):
    """
    Compare the sequences to similarity. same sequences with similar annotations shall be clustered
    :param list_of_sequences:
    :return:
    """
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
records = SeqIO.parse("../../files/vectors.gb", "genbank")

feature_type = 'rep_origin'
jeremyFeatures = ['oriT', 'polyA_signal', 'rep_origin', 'primer_bind', 'rRNA', 'mRNA', 'tRNA']
type_count = dict((e, 0) for e in jeremyFeatures)


for features in jeremyFeatures:
    list_generator = generateList(features, type_count)
    occurences = 0
    list_of_sequences = []
    for resultObject in list_generator:
        print resultObject
        print resultObject.sequence.complement()
        occurences += 1
        list_of_sequences.append(">test\n" + resultObject.sequence + "\n")
    print("occurences of " + features + ": " + str(occurences))

# TODO: same sequences + annotations -> count occurences and prepare new list
