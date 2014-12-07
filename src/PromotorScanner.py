__author__ = 'dominikburri'

from Bio import SeqIO

records = SeqIO.parse("/Users/dominikburri/PycharmProjects/"
                      "Bioinformatik/PlasmidAnalysis/files/vectors-100.gb", "genbank")

for record in records:
    print record.name
