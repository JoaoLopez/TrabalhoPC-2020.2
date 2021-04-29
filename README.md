# Diferenças finitas PVC-1D

## Trabalho de implementação da disciplina Programação Científica - Universidade Federal Fluminense

## Alunos

- Pedro Paulo
- Gabriel Henrique
- Nicholas Quintella
- Gabriel Duarte
- Joao Lopez

<p align="left">Este experimento tem como objetivo desenvolver um algoritmo que resolva numericamente Problemas de Valor de Contorno para Equações Diferenciais pelo Método das Diferenças Finitas.</p>

### Dependências que precisam ser instaladas para executar o projeto:
- sympy: pip install sympy
- numpy: pip install numpy
- latex2sympy: pip install latex2sympy
- antlr4-python3-runtime: pip install antlr4-python3-runtime

### Para executar o método das diferenças finitas PVC: 

* python3 PVC1D.py

Consultar arquivo latexcheatsheet.pdf para um guia rápido da representação utilizada neste projeto.  

## Casos de teste
- teste 1
  f(x) em LaTeX : (0.5*(\sqrt{1+x})) + \arctan{x} + ((\sin{x})^{2}) \n

  a: -1

  b: 2

  n: 9

  p: 1

  q: 2

  y(a): 0

  dy/dx (b): -0.111
- teste 2 (exemplo mostrado em  : https://www.ufrgs.br/reamat/CalculoNumerico/livro-sci/pdvdc-metodo_de_diferencas_finitas.html )
  f(x) em LaTeX : \exp{-x}

  a: 0.5

  b: 1.5

  n: 10

  p: -1

  q: 1

  y(a): 1

  dy/dx (b): 2
