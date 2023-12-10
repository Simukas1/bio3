from Bio import AlignIO
from Bio import SeqIO
from Bio import SearchIO
from Bio.Align.Applications import MafftCommandline
from Bio.Blast import NCBIWWW
from Bio.SeqUtils import IUPACData
from os import path
import numpy as np

import matplotlib.pyplot as plt

def doBlastSearch(filePath):
	if not path.exists("rez.xml"):
		res = NCBIWWW.qblast(program="blastp", database="swissprot", sequence=getSequence(filePath).seq, perc_ident=80, entrez_query='mammals[Organism]')
		open("rez.xml", "w").write(res.read())
		res.close()
	return SearchIO.read("rez.xml", "blast-xml")

def getSequence(filePath):
	for seq in SeqIO.parse(filePath, "fasta"):
		return seq
	
def filter(blasRes):
	ans = []
	with open("data.fasta", "w") as data:
		for d in blasRes:
			if "Albumin" in d.description:
				ans.append(d)
				SeqIO.write(d.hsps[0].hit, data, "fasta")
	return ans

def doAlignment():
	mafft = MafftCommandline("mafft.bat", input = "data.fasta")
	out, _ = mafft()
	open("mafft_msa.fasta", "w").write(out)

def findSeq():
	alignment = AlignIO.read("mafft_msa.fasta", "fasta")
	length = len(alignment[0])
	mini = -1
	maxi = -1
	maxVal = -1e9
	minVal = 1e9
	values = []

	for i in range (length - 14):
		val = 0
		col = alignment[:, i:i+15]
		for j in range(15):
			row = col[:, j]
			aaCount = {}
			for l in IUPACData.protein_letters:
				aaCount[l] = 0
			for l in IUPACData.protein_letters:
				aaCount[l] = aaCount[l] + row.count(l)
			
			for v in aaCount.values():
				if v != 0:
					val += ((v / len(row)) * np.log2((v / len(row)) / 0.05))

		values.append(val)
			
		if val < minVal:
			minVal = val
			mini = i
		if val > maxVal:
			maxVal = val
			maxi = i

	print("Labiausiai skiriančio indeksas: ", end="")
	# pridedu vieną, nes čia indeksas eina nuo nulio, o Jalview eina indeksas nuo vieno.
	print(mini + 1)
	print("Labiausiai skiriančio seka: ", end="")
	print(alignment[0].seq[mini:mini+15])
	print("Labiausiai panašaus indeksas: ", end="")
	# pridedu vieną, nes čia indeksas eina nuo nulio, o Jalview eina indeksas nuo vieno.
	print(maxi + 1)
	print("Labiausiai panašaus seka: ", end="")
	print(alignment[0].seq[maxi:maxi+15])
	plt.plot(range(0, len(values)),values)
	plt.ylabel('information content')
	plt.xlabel('position')
	plt.show()

ans = filter(doBlastSearch("alb_human.fasta"))
doAlignment()
findSeq()
