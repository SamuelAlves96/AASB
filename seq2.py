# -*- coding: utf-8 -*-

class Seq:
    
    def __init__(self,sequencia):
        self.sequencia = sequencia.upper()
        self.alfabeto = ""
        self.tipo = None
        self.propriedades = {}
        
        
    def valida(self):
        for c in self.sequencia:
            if c not in self.alfabeto:
                return False  
        return True

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, n):
        return self.sequencia[n]


    def __str__(self):
        return self.tipo + ":" + self.seq

    def printseq(self):
        print(self.seq)


    def get_propriedade(self, key):
        if key in self.propriedades.key():
            return self.propriedades[key]
        else:
            raise KeyError('chave inexistent')
            
    def add_propriedade(self,key,value):
        self.propriedades[key]=value
        
        

class SeqDNA(Seq):
    
    def __init__(self, sequencia):
        super.__init__(sequencia)
        self.tipo = "DNA"
        self.alfabeto ="ACTG"
        if not self.valida():
            raise ValueError("A sequencia invalida")
        
        


class SeqRNA(Seq):
    def __init__(self, sequencia):
        super.__init__(sequencia)
        self.alfabeto ="ACUG"
        
        
        


class SeqAA(Seq):
    def __init__(self, sequencia):
        super.__init__(sequencia)
        self.alfabeto ="ACTG"