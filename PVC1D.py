
from sympy import *
from latex2sympy import *
from numpy.linalg import inv
import numpy as np

FFX=latex2sympy(r"\frac {1+ \sqrt {\a}} {\b}",ConfigL2S())
FF = parsing.sympy_parser.parse_expr("(1*(1*(1*a)**0.5 + 1*1))/((1*b))")


#TESTADO OK
def discretiza_dominio1D(a,b,n):
    '''Dado um dominio delimitado pelo intervalo (a,b)
discretiza esse domínio, dividindo ele em n pontos equidistantes
ou seja, dado um dominio [a,b] , e um valor n retorna os n pontos
desse domínio cuja distancia entre eles, é (a-b)/n ; incluindo a e b
visão gráfica em :
https://drive.google.com/file/d/1bzIc99DaNZOIob1cHPAraBw7MMCe2ECK/view?usp=sharing
'''
    if n<=1:
        return None
    # tem q ter pelo menos 2 pontos
    x=[]
    h=(b-a)/(n)
    for i in range(0,n):
        x.append(a+(i*h))
    x.append(b)
    return x


def calcula_f(funcao_em_lateX,x_i):
    """"retorna f(x_i) = f_i ; dados o ponto x_i e a string da representação de f em LaTeX"""
    pass

def get_vetor_f(funcao_f,y0,yn,pontos_dominio,p,h,n):
    """"retorna o vetor f que contém os valores de f nos pontos do dominio discretizado"""
    # primeiro valor é f_1 + y0*p/h²
    resp=[calcula_f(funcao_f,pontos_dominio[1])+(p*y0/(h**2))]
    #itera sobre os pontos xk do dominio , calcula f(xk) e guarda no vetor
    for k in range(2,n):
        xk=pontos_dominio[k]
        resp.append(calcula_f(funcao_f,xk))
    return resp


def monta_matriz_D(n,p,q,h):
    """monta a matriz do PVC aproximado por diferenças finitas, dado o n, que vale...
    [ q+(2p/h²)     -p/h²     0      0       ...      0           0   ]    1
    [ -p/h²      q+(2p/h²)   -p/h²   0       ...      0           0   ]    2
    [   0        -p/h²    q+(2p/h²) -p/h²    ...      0           0   ]    3
    [   0         0         -p/h²  q+(2p/h²) ...      0           0   ]    4
    [   .         .          .       .      .         .           .   ]    .
    [   :         :          :       :        .       :           :   ]    :
    [   0         0          0       0       ...   q+(2p/h²)    -p/h  ]   n-2
    [   0         0          0       0       ...   -p/h       q+(p/h²)]   n-1
    """
    matriz = []
    for i in range(n-1):
        coluna = []
        
        for j in range(n-1):
            if(i == j):
                coluna.append(2*p/(h*h) + q)
            elif(abs(i - j) == 1):
                coluna.append(-p/(h*h))
            else:
                coluna.append(0)
        
        matriz.append(coluna)

    matriz[-1][-1] = p/(h*h) + q
    
    return matriz

def PVC_1D(n,p,q,h,funcao_f,y0,yn,pontos_dominio):
    """ Assumindo que a matriz do PVC aprox. por dif. finitas é dada por D
       o vetor f com os valores de f nos pontos do dominio discretizado dado por F
       se quer descobrir o vetor Y, onde vale a equação matricial a seguir...
       DY = F
       aplicando propriedades de alg. linear temos que
       (D-¹)DY = (D-¹)F
        ONDE (D-¹) é a inversa de D, e então
       Y=(D-¹)F
       E ASSIM TEREMOS A SOLUÇÃO
       """
    D_inv=inv(np.matrix(np.array(monta_matriz_D(n,p,q,h))))  #calcula a inversa
    F=np.transpose(get_vetor_f(funcao_f,y0,yn,pontos_dominio,p,h,n))  #obtem o vetor vertical F
    #retorna o produto matricial
    return D_inv.dot(F)

while(True):
    a = float(input("a: "))
    b = float(input("b: "))
    n = int(input("n: "))
    p = float(input("p: "))
    q = float(input("q: "))
    h = (b-a)/n

    print("DISCRETIZAÇÃO:")    
    print(discretiza_dominio1D(a,b,n))
    print("\n")

    matriz = monta_matriz_D(n, p, q, h)
    print("MATRIZ:")
    for linha in matriz:
        print(linha)
    print("\n\n")

latex2sympy()
