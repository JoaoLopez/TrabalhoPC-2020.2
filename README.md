# Diferenças finitas PVC-1D

## Trabalho de implementação da disciplina Programação Científica - Universidade Federal Fluminense

## Alunos

- Pedro Paulo
- Gabriel Henrique
- Nicholas Quintella
- Gabriel Duarte
- Joao Lopez

<p align="left">Este experimento tem como objetivo desenvolver um algoritmo que resolva numericamente Problemas de Valor de Contorno para Equações Diferenciais pelo Método das Diferenças Finitas.</p>

Dependências que precisam ser instaladas para executar o projeto:
- sympy: pip install sympy
- numpy: pip install numpy
- latex2sympy: pip install latex2sympy

Para executar o método das diferenças finitas PVC: 

* python3 PVC1D.py

Consultar arquivo latexcheatsheet.pdf para um guia rápido da representação utilizada neste projeto.  

## Caso de teste

f(x) em LaTeX : (0.5*(\sqrt{1+x})) + \arctan{x} + ((\sin{x})^{2}) \n

a: -1

b: 2

n: 9

p: 1

q: 2

y(a): 0

dy/dx (b): -0.111

