#Make sure to install Biopython library using : pip install biopython
from Bio import Entrez
Entrez.email = "lyacine230@gmail.com"
handle = Entrez.esearch(db="nucleotide", term='"Bardaxima"[Organism] AND ("2004/07/24"[PDAT] : "2010/04/28"[PDAT])') 
record = Entrez.read(handle)
print(record["Count"])
