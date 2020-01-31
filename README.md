# Trabalho prático da Unidade Curricular Algoritmos para Análise de Sequências Biológicas

Este trabalho consiste em um programa de análise de sequências de DNA desenvolvidos pelo grupo 5.

## Explicação do programa

O trabalho está dividido em dois ficheiros (engineSeq e shellSeq).

O engineSeq é composto por cinco classes:

* ```Seq```, em que definimos a função valida, dos índices, de leitura de um ficheiro, do download de um ficheiro fasta do NCBI, do blast, do alinhamento e da árvore filogenética.

* ```SeqDNA``` valida a sequência inserida (se é ou não de DNA). Definimos funções para a transcrição, complemento inverso, tradução, ORFs e para a procura de um padrão/frequência de símbolo(s) na sequência.

* ```SeqRNA``` valida a sequência inserida (se é ou não de RNA).

* ```SeqAA``` valida a sequência inserida (se é ou não de Aminoácidos).

* ```DB```, em que definimos a função de adicionar chaves às sequências, das keys (para nos fornecer informação de quantas chaves existem), do get (para nos dar a informação da sequência com determinada chave), do guardar (para guardar o ficheiro), do serializa (para guardar o ficheiro no pickle, em formato binário), e o deserializa (para ir buscar o ficheiro ao pickle e o retirar de formato binário).
