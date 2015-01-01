from Bio import SeqIO
from bioservices import *
from Bio import motifs
from Bio.Alphabet import IUPAC, Gapped
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

def generateList(feature_type, filePath):
    """
    A generator
    :param feature_type:
    :return: a ResultObject with the desired sequence and annotation
    """
    records = SeqIO.parse(filePath, "genbank")

    for record in records:
        if len(record.seq) > 1500: # minimum for number of bases
            for feature in record.features:
                if feature.type == feature_type:
                    sequence_of_feature = record.seq[feature.location.start: feature.location.end]
                    annotation = feature.qualifiers
                    importantAnnotation = dict((k,annotation[k]) for k in ('note','gene','bound_moiety','mobile_element_type','product')if k in annotation)
                    feature_type = feature.type
                    result = ResultObject(sequence_of_feature, feature_type, importantAnnotation)
                    yield result

def clustering(objects_of_sequences):
    """
    MUSCLE
    Compare the sequences to similarity. same sequences with similar annotations shall be clustered
    :param objects_of_sequences:
    :return:
    """
    print 'MUSCLE'
    list_of_sequences = ""

    if len(objects_of_sequences)<=1:
        return []

    for (i,k) in enumerate(objects_of_sequences):
        #print k
        list_of_sequences += ">" + "identifier" + str(i) +"\n"+str(k.sequence)+"\n\n"

    #print list_of_sequences
    m = MUSCLE(verbose=False)
    jobid = m.run(frmt="fasta", sequence=list_of_sequences, email="dominik.burri1@students.fhnw.ch")

    while m.getStatus(jobid) == u'RUNNING':
        print "Status: ", m.getStatus(jobid)

    result=m.getResult(jobid, "sequence")
    #print "sequence:"
    #print result
    sequencelist = result
    f = open('sequence_result.fasta', 'w')
    f.write(sequencelist)
    f.close()

    result=m.getResult(jobid, "aln-fasta")
    #print "aln-fasta:"

    #print result
    #print "Resulttypes: ", m.getResultTypes(jobid)

    result=m.getResult(jobid, "phylotree")
    #print "phylotree:"
    #print result
    result=m.getResult(jobid, "pim")
    #print "pim:"
    #print result
    pim_result = result
    f = open('pim_result.txt', 'w')
    f.write(pim_result)
    f.close()

    result=m.getResult(jobid, "out")
    #print "out:"
    #print result

    return sequencelist

def pim_evaluation(schwellenwert):
    '''
    Auswertung der Percent Identity Matrix
    :param schwellenwert: der Schwellenwert fuer die Erkennung von Matches
    :return: Liste mit aehnlichen Sequenzen (als Seq Object gespeichert), die jeweils in eine Liste gepackt sind
    '''

    identifier_list = []
    matches = []

    f = open('pim_result.txt', 'r')
    for i in range(6):
        f.readline()
    lines = 1
    while True:
        line = f.readline()
        if line == '':
            break
        words = line.split()
        words.pop(0) # deleting 1: etc
        name = words[0] # getting 'identifierXY'
        identifier_list.append(name)
        words.pop(0) # deleting 'identifierXY'
        if len(words) >= 1:
            index = 1
            for value in words:
                value = float(value)
                if index < lines:
                    if value > schwellenwert:
                        matches.append([name, index])
                        #print name, secondname
                else:
                    break
                index += 1
        lines += 1

    # get the correct index from the full identifier list
    # and set the name of the corresponding identifier
    names = []
    new_matches = []
    for match in matches:
        match[1] = identifier_list[match[1]-1]

        if match[1] in names:
            new_matches.append(match[0])
        else:
            names.append(match[1])

        print match # print the identifier names
        # if not match[0] in names:
        #     names.append(match[0])
        # if not match[1] in names:
        #     names.append(match[1])

    # TODO: get the multiple sequences that are similar
    # multiple_similar_sequences = []
    # for new_match in new_matches:
    #     multiple_similar_sequences.append(new_match)
    #     for match in matches:
    #         if new_match in match:
    #             for entry in match:
    #                 if not entry in multiple_similar_sequences:
    #                     multiple_similar_sequences.append(entry)
    #
    # print 'Multiple similar sequences: ' + str(multiple_similar_sequences)

    # TODO: get the unnessecary entries out
    # for match in matches:
    #     print match
    #     if match[0] in final_list:
    #         matches.remove(str(match))
    #     if match[1] in final_list:
    #         matches.remove(str(match))

    #matches.append(multiple_similar_sequences)

    print 'Matches: ' + str(matches)
    # get the sequence from the identifier name
    handle = open('sequence_result.fasta', 'r')
    for record in SeqIO.parse(handle, 'fasta', IUPAC.unambiguous_dna):
        for match in matches:
            for i in range(len(match)):
                if record.id == match[i]:
                    match[i] = record.seq

            #print match

    handle.close()
    return matches


def createPSSM(sequencelist):
    print "Start PSSM"

    if len(sequencelist)==0:
        return

    #sequencelist = sequencelist.replace("-", ".")
    f = open('fastatmp', 'w')
    f.write(sequencelist)
    f.close()

    list = []

    for seq_record in SeqIO.parse("fastatmp", "fasta", IUPAC.unambiguous_dna):
        list.append(str(seq_record.seq))

    #Blast typical sequence
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
    return pssm

def reduce_to_single_sequences(generated_object, feature):
    # TODO: implement in 'generateList'
    """
    same sequences + annotations -> count occurences and prepare new list
    :param generated_object:
    :return:
    """

    featureTypes = {
                    'oriT': ['gene', 'product'], 'polyA_signal': ['note'], 'rep_origin': ['note'],
                    'primer_bind': ['note'], 'rRNA': ['poduct'], 'mRNA': ['gene'], 'tRNA': ['product'],
                    'promoter': ['note'], "RBS": ['note', 'gene'], "-10_signal": ['note', 'gene'],
                    '-35_signal': ['note', 'gene'], 'terminator': ['note'], 'CDS': ['gene', 'product'],
                    'protein_bind': ['note', 'bound_moiety'], 'misc_binding': ['note', 'bound_moiety'],
                    'misc_recomb': ['note'], 'LTR': ['note'], 'misc_signal': ['note'], 'enhancer': ['note'],
                    'mobile_element': ['mobile_element_type', 'note'], 'sig_peptide': ['note']
                    }

    results = []
    try:
        results.append(generated_object.next())
    except (StopIteration):
        print "Warning: empty generator. ", feature, " not found!"
        return []
    #print results

    for resultObject in generated_object:
        counter = 1
        foundMatch = False
        for result in results:
            matchCounter = 0
            for key in featureTypes[feature]:
                if str(resultObject.sequence) == str(result.sequence) and resultObject.annotation.get(key)==result.annotation.get(key):
                    #print str(resultObject.annotation.get(key))
                    matchCounter += 1
            if matchCounter == len(featureTypes[feature]):
                result.setOccurences()
                foundMatch = True
            if len(results) == counter:
                if foundMatch == False:
                    results.append(resultObject)
            counter += 1

    return results




# - - - - start of skript - - - -
# - - - - - - - - - - - - - - - -
jeremyFeatures = ['oriT', 'polyA_signal', 'rep_origin', 'primer_bind', 'rRNA', 'mRNA', 'tRNA']
dominiks_list = ['promoter', 'RBS', '-10_signal', '-35_signal']
kevins_list = ['terminator', 'CDS']
alessandros_list = ['protein_bind', 'misc_binding', 'misc_recomb', 'LTR', 'misc_signal',
                    'enhancer', 'mobile_element', 'sig_peptide']


complete_list = jeremyFeatures + dominiks_list + kevins_list + alessandros_list

save_file_object = open("list_of_identical_objects.txt", "w")

# Schwellenwert fuer nahezu identische Sequenzen bei der percent identitiy matrix
schwellenwert = 90.0

for feature in jeremyFeatures:
    print 'Feature: ' + feature
    filePath = "../../files/vectors.gb"
    # make a list generator with the desired feature and its annotation
    list_generator = generateList(feature, filePath)

    # for eintrag in list_generator:
    #     for key in featureTypes[feature]:
    #
    #         print eintrag.annotation.get(key)


    #  same sequences + annotations -> count occurences and prepare new list
    list_of_identical_objects = reduce_to_single_sequences(list_generator, feature)
    summe = 0
    for object in list_of_identical_objects:
        summe += object.getOccurences()
        #Blast typical sequence
        save_file_object.write(str(object) + "\t" + str(object.getOccurences()) + "\n")
    print("Anzahl identischer objekte: \t" + str(len(list_of_identical_objects)))
    print("Summe aller Objekte: \t\t\t" + str(summe))
    # TODO: 'wichtige Annotation' Sequenzen in Liste speichern und MUSCLE uebergeben
    muscle_result = clustering(list_of_identical_objects)
    # PIM Auswertung: Sequenzen groesser Schwellenwert (bsp. 95%) rausspeichern. Rueckgabe: Liste von "fast identische Sequenzen"
    list_of_near_identical_sequences = pim_evaluation(schwellenwert)
    for sequences in list_of_near_identical_sequences:
        print sequences
    # TODO: neues MUSCLE
    # TODO: PSSM
    #createPSSM(sequencelist)
save_file_object.close()