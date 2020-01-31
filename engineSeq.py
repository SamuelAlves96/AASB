# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 17:09:29 2020

@author: faculdade
"""

from Bio import *
from collections import Counter
import pickle



class Seq:
    
    def __init__(self,sequencia):
        self.sequencia = sequencia.upper()
        self.alfabeto = ""
        self.propriedades = {}
        
            
        
    def valida(self):
        for c in self.sequencia:
            if c not in self.alfabeto:
                return False  
        return True

    def __len__(self):
        return len(self.sequencia)

    def __getitem__(self, n):
        return self.sequencia[n]


    def __str__(self):
        return self.sequencia
    
    def __repr__(self):
        return self.sequencia

    def printseq(self):
        print(self.sequencia)

    def get_propriedade(self, key):
        if key in self.propriedades.key():
            return self.propriedades[key]
        else:
            raise KeyError('Chave inexistente')
            
    def add_propriedade(self,key,value):
        self.propriedades[key]=value
        

    def ler_ficheiro(filename):
        '''
        Leitura de ficheiro de texto ou em formato FASTA
        '''
        ficheiro = []
        with open(filename, "r") as f:
            ls = f.readlines()
            for i in range (1,len(ls)):
                 s = ls[i].rstrip("\n")
                 if len(s)>0:
                     ficheiro.append(s)
        return ficheiro
        
        
        
    def download_NCBI(search_term):
        Entrez.email = 'diogomfbarbosa.93@gmail.com'
        handle = Entrez.esearch(db='nucleotide', term=search_term)
        seq_ids = Entrez.read(handle)['IdList']
        for seq_id in seq_ids:
            record = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            handle.close()
            filename = '%s.fasta' % seq_id
            print('Writing:{}'.format(filename))
            with open(filename, 'w') as f:
                f.write(record.read())
            record = SeqIO.read('%s.fasta' % seq_id, "fasta")
        return search_term
    
    
    def Blast(filename):
        record = SeqIO.read("%s.fasta" % filename, format="fasta")
        print(record)
        result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
        with open("blast.xml", "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()
        result_handle=open("blast.xml")
        blast_records=NCBIXML.read(result_handle)
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.05:
                    print("****Alignment****")
                    print("sequence:", alignment.title)
                    print('length:', alignment.length)
                    print('e value:', hsp.expect)
                    print(hsp.query[0:75] + '...')
                    print(hsp.match[0:75] + '...')
                    print(hsp.sbjct[0:75] + '...')
        return filename
    
    
    def Align(filename):
        clustexe=r'C:\Program Files (x86)\ClustalW2\clustalw2.exe'
        clustalw_cline = ClustalwCommandline(clustexe, infile="filename.fasta")
        clustalw_cline()
        return filename
    
    
    def Phylo(filename):     
        tree = Phylo.read("filename.dnd", "newick")
        Phylo.draw_ascii(tree)
        return filename    

        

class SeqDNA(Seq):
    
    def __init__(self, sequencia):
        super(SeqDNA,self).__init__(sequencia)
        self.tipo = "DNA"
        self.alfabeto ="ACTG"
        if not self.valida():
            raise ValueError(f"A sequencia {sequencia} é invalida")
        
        
    def testedna(self):
        
        print('Teste DNA')
    

    def transcricao(self):
        self.tipo='rna'
        return SeqRNA(self.sequencia.replace('T','U'))

    def compInverso(self):
        comp = ""
        for c in self.sequencia:
            if c == "A":
                comp = "T" + comp
            elif c == "T":
                comp = "A" + comp
            elif c == "C":
                comp = "G" + comp
            elif c == "G":
                comp = "C" + comp
            else:
                raise ValueError("Sequencia invalida!!!!")
        return SeqDNA(comp)

    def traduzSeq(self, iniPos=0):

        seqM = self.sequencia
        seqAA = ""
        for pos in range(iniPos, len(seqM)-2, 3):
            cod = seqM[pos:pos+3]
            aa= self.traduzCodao(cod)
            seqAA = seqAA + aa
        return SeqAA(seqAA)

    def orfs(self):  #open reading frames

        res = []
        res.append(self.traduzSeq(0))
        res.append(self.traduzSeq(1))
        res.append(self.traduzSeq(2))
        compinv = self.compInverso()
        res.append(compinv.traduzSeq(0))
        res.append(compinv.traduzSeq(1))
        res.append(compinv.traduzSeq(2))
        return res

    def traduzCodao(self, cod):
        tc = {"GCT": "A", "GCC": "A", "GCA": "A", "GCC": "A", "TGT": "C", "TGC": "C",
              "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TTT": "F", "TTC": "F",
              "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", "CAT": "H", "CAC": "H",
              "ATA": "I", "ATT": "I", "ATC": "I",
              "AAA": "K", "AAG": "K",
              "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
              "ATG": "M", "AAT": "N", "AAC": "N",
              "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
              "CAA": "Q", "CAG": "Q",
              "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
              "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
              "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
              "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
              "TGG": "W",
              "TAT": "Y", "TAC": "Y",
              "TAA": "_", "TAG": "_", "TGA": "_"}
        if cod in tc:
            aa = tc[cod]
        else:
            aa = "X"  # errors marked with X
        return aa


    # def frequencia(self, exp):
        # from shellSeq import shell
        # A = 0
        # C = 0
        # T = 0
        # G = 0
        # cont=0
        # for c in letra:
        #     if c.isalpha():
        #         cont+=1
        #     else:
        #         pass
        # print (f'Quantos {letra} tem: ', cont)
        
        #freqs = Counter(''.join(self.sequencia.splitlines()))
        #for symbol, count in freqs.most_common():
          #  print (symbol, count)
        
    def procuraPadraoER(seq, pad):
        from re import search
        mo = search(pad, seq)
        if (mo != None):
            return mo.span()[0]
        else:
            return -1

    def procuraTodasOcsER(seq, pad):
        from re import finditer
        mos = finditer(pad, seq)
        res = [ ]
        for x in mos:
            res.append(x.span()[0])
        return res
    
    
    
    
    
class SeqRNA(Seq):
    def __init__(self, sequencia):
        super(SeqRNA,self).__init__(sequencia)
        self.alfabeto ="ACUG"
        self.tipo = "RNA"
        if not self.valida():
            raise ValueError("A sequencia é invalida")
        













class SeqAA(Seq):
    def __init__(self, sequencia):
        super(SeqAA,self).__init__(sequencia)
        self.tipo = "AA"
        self.alfabeto ="ACDEFGHIKLMNPQRSTVWY_"
        if not self.valida():
            raise ValueError(f"A sequencia {sequencia} é invalida")





class DB:
    def __init__(self):
        self.dicionario = {}
        self.ficheiro_ser = "serializacao.dat"
        
        
    def adicionar(self,seq, key = None):
        if key is None:
            key = len(self.dicionario)
        self.dicionario[key]=seq
        print(f"Sequencia foi adicionada com a chave {key}")
        return key
        
    def __len__(self):
        return len(self.dicionario)
    
    def keys(self):
        return [k for k in self.dicionario.keys()]
    
    
    def get(self, key):
        return self.dicionario[key]
    
    
    def guardar(self,ficheiro):
        with open (ficheiro, 'w') as f:
            for key, sequencia in self.dicionario.items():
                f.write(f"{sequencia.seq} \n")
    
    
    def serializa(self):
        pickle.dump(self.dicionario, open( self.ficheiro_ser, "wb"))
    
    
    def deserializa(self):
       try:
           self.dicionario = pickle.load(open( self.ficheiro_ser, "rb"))
       except:
           pass
       
            




