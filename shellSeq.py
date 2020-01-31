# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 17:14:50 2020

@author: faculdade
"""

from cmd import *
from engineSeq import *
from collections import  OrderedDict
from MatrixNum import MatrixNum
from ClustHier import ClustHier
from MySeq import MySeq
from AlignSeq import AlignSeq
from SubstMatrix import SubstMatrix
from upgma import *

db = DB()
db.deserializa()

class shell(Cmd):
    global intro
    intro = '\nInterpretador de comandos de análise de sequências biológicas. \n \nPara inserir a sequência, escreva inserir. \nPara ler uma sequência inserida num ficheiro, escreva ler. \nPara obter o tamanho do dicionário, escreva len. \nPara obter as keys, escreva keys. \nPara obter a sequência de determinada key, escreva get. \nPara realizar alinhamentos múltiplos, escreve align. \nPara fazer o download de um ficheiro, escreva download. \nPara fazer o blast, escreva blast. \nPara calcular a frequência de um símbolo, escreva frequencia. \nPara procurar ocorrências de um padrão, escreva padrao. \nPara construir uma árvore filogenética, escreva phylo. \nPara construir uma árvore usando o algoritmo UPGMA, escrever upgma.'
    print(intro)
    prompt = 'Seq> ' 


    def menu(self):
        print(intro)
     
    def nova_seq(self,sequencia):
        dna = SeqDNA(sequencia)
        db.adicionar(dna)
        
    def do_len(self,arg):
        print(len(db))
        self.menu()
        
    def do_keys(self,arg):
        print(db.keys())
        self.menu()
        
    def do_get(self,arg):
        i= input('Qual o índice? ')
        print(db.get(int(i)))
        self.menu()
        
    def do_sair(self, arg):
 
        """Sair do programa: sair"""
        db.serializa()
        print('Obrigado!')
        global janela  # pois pretendo atribuir um valor a um identificador global
        if janela is not None:
                    del janela  # invoca o metodo destruidor de instancia __del__()
        return True

    def do_inserir(self,arg):
        lista=[]
        global sequencia
        sequencia = input("Sequencia: ")
        # try:
        s1 = SeqDNA(sequencia)
        # db.adicionar(s1)
        s1.printseq()
        # SeqDNA.testedna(self)
        print("Transcricao:")
        transcricao = s1.transcricao()
        print(transcricao)
        print("Traducao:")
        traducao = s1.traduzSeq()
        print(traducao)
        # except ValueError:
        print ("ORFs:")
        orf=s1.orfs()
        listorf=[]
        for x in orf:
            listorf.append(x)
            print(x)
        op = input('Quer guardar as ORFs? ')
        if op == 'sim':
            ficheiro = input('Indique o nome do ficheiro: ')
            file_csv = ficheiro
            with open(file_csv,'w') as f:
                f.write(f"{listorf}{';'}")
                
        elif op == 'nao':
            self.menu()
            return
        # print('Erro na sequência inserida!')
        lista.extend((s1, transcricao, traducao, orf))
        db.adicionar(lista)
        self.menu()

    def do_ler(self, arg):
        fasta = Seq.ler_ficheiro(arg)
        lista_fasta=[]
        for f in fasta:
            s2 = SeqDNA(f)
            s2.printseq()
            lista_fasta.append(s2)
        db.adicionar(lista_fasta)
        self.menu()
    
    def do_frequencia(self,arg):
        # seq = input("Indice:")
        seq = db.get(int(input("Indice: ")))
        pad = input("Simbolo: ")
        todos_res = SeqDNA.procuraTodasOcsER(str(seq), pad)
        if len(todos_res) > 0:
            print("Simbolo descoberto nas posicoes: ", todos_res)
            print(pad,len(todos_res))
        else:
            print("Simbolo nao encontrado")
        self.menu()
    
    def do_padrao(self,arg):
        seq = db.get(int(input("Indice: ")))
        pad = input("Padrao: ")
        todos_res = SeqDNA.procuraTodasOcsER(str(seq), pad)
        if len(todos_res) > 0:
            print("Padrao descoberto nas posicoes: ", todos_res)
            print(pad,len(todos_res))
        else:
            print("Padrao nao encontrado")
        self.menu()
    
    def do_upgma(self,arg):
        i = int(input('Quantas sequências quer? '))
        listaupgma = []
        j=int(input('Qual o valor do mismatch? '))
        l=int(input('Qual o valor do gap? '))
        for k in range(0,i):
            a = db.get(int(input("Indice: ")))
            listaupgma.append(a)
        sm=SubstMatrix()
        sm.createFromMatchPars(i, j, "ACTG")
        alseq = AlignSeq(sm,l)
        up = UPGMA(str(listaupgma), alseq)
        arv = up.run()
        arv.printtree()
        self.menu()
    
        
    def do_download(self, arg):
        Seq.download_NCBI(arg)
        self.menu()
        
    def do_blast(self, arg):
        Seq.Blast(arg)
        self.menu()
    
    def do_align(self, arg):
        Seq.Align(arg)
        self.menu()
    
    def do_phylo(self, arg):
        Seq.Phylo(arg)
        self.menu()


    
if __name__ == '__main__':
    janela = None


    sh = shell()

    sh.cmdloop()