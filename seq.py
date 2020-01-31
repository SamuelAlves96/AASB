# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 15:43:05 2019

@author: gil-a
"""

from abc import ABC

class Seq(ABC):
    def __init__(self,sequencia, tipo="dna"):
        """
        Cria uma nova sequencia
        argumentos:
            sequencia (str)
        """
        self.sequencia = sequencia.upper()
        self.alfabeto = "ACTG"
        self.tipo = tipo
        
    def __eq__(self,obj):
        if isinstance(obj,Seq):
            return self.sequencia == obj.sequencia
        else:
            return False
        
    def __str__(self):
        return self.sequencia
    
    def __getitem__(self, index):
        return self.sequencia[index]
        
    def valida(self):
        """
        return if the sequence is valid
        """
        for i in range(len(self.sequencia)):
            if self.sequencia[i] not in self.alfabeto:
                return False
        return True
    
    def validaER(self):
        import re
        if (self.tipo == "dna"):
            if re.search("[^ACTGactg]",self.sequencia) != None:
                return False
        elif (self.tipo == "rna"):
            if re.search("[^ACUGacug]",self.sequencia) != None:
                return False
        elif (self.tipo == "protein"):
            if re.search("[^ACDEFGHIKLMNPQRSTVWY_acdefghiklmnpqrstvwy]", self.sequencia) != None:
                return False
        return True
    
    
if __name__ == "__main__":
    sequencia = Seq('AcggTcGGG')
    print(sequencia.valida())
    s1 = Seq(sequencia)
    s1.printsequencia()

    if s1.validaER():
        print("Sequencia valida")
        # print("Complemento inverso:")
        # s1.compInverso().printseq()
        # print("Traducao: ")
        # s1.traduzSeq().printseq()
        # print("ORFs:")
        # for orf in s1.orfs():
        #     orf.printseq()
        # print("Maior proteina nas ORFs:")
        # s1.maiorProteinaORFs().printseq()
    else:
        print("Sequencia invalida")

