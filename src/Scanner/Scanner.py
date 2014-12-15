

from Bio import SeqIO
from bioservices import *
from Bio import motifs
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW


class ResultObject:
    """
    An Object for storing the sequence, feature type and annotation of the feature
    """
    def __init__(self, sequence, feature_type, annotation):
        self.occurences = 0
        self.annotation = annotation
        self.feature_type = feature_type
        self.sequence = sequence
    def __str__(self):
        return str(self.feature_type)+"; " + str(self.sequence)+ "; " + str(self.annotation)

    def setOccurences(self):
        self.occurences += 1
    def getOccurences(self):
        return self.occurences

def generateList(feature_type):
    """
    A generator
    :param feature_type:
    :return: a ResultObject with the desired sequence and annotation
    """
    for record in records:
        if len(record.seq) > 1500: # minimum for number of bases
            for feature in record.features:
                if feature.type == feature_type:
                    sequence_of_feature = record.seq[feature.location.start: feature.location.end]
                    annotation = feature.qualifiers
                    feature_type = feature.type
                    result = ResultObject(sequence_of_feature, feature_type, annotation)
                    yield result

def clustering(objects_of_sequences):
    """
    Compare the sequences to similarity. same sequences with similar annotations shall be clustered
    :param list_of_sequences:
    :return:
    """
    list_of_sequences = ""

    for k in objects_of_sequences:
        print k
        list_of_sequences += ">" + "identifier" +"\n"+str(k.sequence)+"\n\n"
        #list_of_sequences='>test_a\nAGAGAGAGAG\n\n>test_b\nAGAAAGAA\n\n>test_c\nAGAGGAGAG\n\n'

    m = MUSCLE(verbose=False)
    jobid = m.run(frmt="fasta", sequence=list_of_sequences, email="dominik.burri1@students.fhnw.ch")

    while m.getStatus(jobid) == u'RUNNING':
        print "Status: ", m.getStatus(jobid)

    result=m.getResult(jobid, "sequence")
    print "sequence:"
    print result


    result=m.getResult(jobid, "aln-fasta")
    print "aln-fasta:"
    sequencelist = result

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

    return sequencelist

def createPSSM(sequencelist):
    print "Start PSSM"

    #sequencelist = sequencelist.replace("-", ".")
    f = open('fastatmp', 'w')
    f.write(sequencelist)
    f.close()

    list = []

    for seq_record in SeqIO.parse("fastatmp", "fasta", IUPAC.unambiguous_dna):
        list.append(str(seq_record.seq))


    result_handle = NCBIWWW.qblast("blastn", "nt", list[0])
    save_file = open("my_blast.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

    #motifs.create(test, alphabet=Gapped(IUPAC.unambiguous_dna))
    m = motifs.create(list, alphabet=Gapped(IUPAC.unambiguous_dna))

    print "motif created"


    pwm = m.counts.normalize(pseudocounts=0.25)
    print "PWM done"
    pssm = pwm.log_odds()
    print "PSSM done"
    print pssm

    return pssm, list



def compare_sequences_and_annotations(generated_object):
    """
    same sequences + annotations -> count occurences and prepare new list
    :param generated_object:
    :return:
    """
    results = []

    results.append(generated_object.next())

    #print results

    for resultObject in generated_object:
        counter = 1
        foundMatch = False
        for result in results:
            if str(resultObject.sequence) == str(result.sequence) and str(resultObject.annotation)==str(result.annotation):
                foundMatch = True
                result.setOccurences()
            if len(results) == counter:
                if foundMatch == False:
                    results.append(resultObject)
            counter += 1

    return results




# - - - - start of skript - - - -
# - - - - - - - - - - - - - - - -
jeremyFeatures = ['oriT', 'polyA_signal', 'rep_origin', 'primer_bind', 'rRNA', 'mRNA', 'tRNA']
records = SeqIO.parse("../../files/vectors-100.gb", "genbank")

# make a list generator with the desired feature and its annotation
feature_type = 'rep_origin'
list_generator = generateList(feature_type)
# occurences = 0
# for resultObject in list_generator:
#     print resultObject
#     occurences += 1
# print("occurences of " + feature_type + ": " + str(occurences))

#  same sequences + annotations -> count occurences and prepare new list
list_of_identical_objects = compare_sequences_and_annotations(list_generator)
summe = 0
for i in list_of_identical_objects:
    #print i, i.getOccurences()
    #print i.sequence
    summe += i.getOccurences()
print len(list_of_identical_objects)
sequencelist = clustering(list_of_identical_objects)
pssm, list = createPSSM(sequencelist)

print summe



# TODO: PSSM only, if the annotations are the same

# TODO: Annotation statistics for every cluster via feature.qualifiers dictionary.

# TODO: Blast of typical sequence against nucleotide database.

#clustering(list_of_sequences)
