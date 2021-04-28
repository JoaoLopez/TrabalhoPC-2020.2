import math
from sympy import *
from sympy.parsing.latex import parse_latex
from latex2sympy import *
from numpy.linalg import inv
import numpy as np
#from process_latex import *

#FFX=latex2sympy(r"\frac {1+ \sqrt {\x}} {2}",ConfigL2S())
#FF = parsing.sympy_parser.parse_expr("(1*(1*(1*a)**0.5 + 1*1))/((1*b))")


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

#DEPRECIADO
def calcula_f_OLD(funcao_em_lateX,x_i):
    """"[DEPRECIADO] retorna f(x_i) = f_i ; dados o ponto x_i e a string da representação de f em LaTeX"""
    # exemplo : 0.5* ( √(x+1) ) → em latex → \frac {\sqrt{1+x}} {2}
    FLX_no_ponto=funcao_em_lateX.replace("x",str(x_i))
    FFX = latex2sympy(FLX_no_ponto, ConfigL2S())
    valor = parsing.sympy_parser.parse_expr(str(FFX))
    return valor
    #pass

# TESTADO OK
def calcula_f(funcao_em_lateX,x_i):
    """"retorna f(x_i) = f_i ; dados o ponto x_i e a string da representação de f em LaTeX"""
    # exemplo : arctg(x) + sen²(x) + 0.5* ( √(x+1) ) → em latex → \arctan{x} + (\sin{x})^{2}) + \frac {\sqrt{1+x}} {2}
    return parse_latex(funcao_em_lateX).evalf(subs=dict(x=x_i))

#VERIFICAR PARTE TEÓRICA (SE ESTÁ PEGANDO OS VALORES QUE DEVERIA SEGUNDO O SIDE)
def get_vetor_f(funcao_f,y0,pontos_dominio,p,h,n):
    """"retorna o vetor f que contém os valores de f nos pontos do dominio discretizado"""
    # primeiro valor é f_1 + y0*p/h²
    resp=[calcula_f(funcao_f,pontos_dominio[1])+(p*y0/(h**2))]
    #itera sobre os pontos xk do dominio , calcula f(xk) e guarda no vetor
    for k in range(2,n):
        xk=pontos_dominio[k]
        resp.append(calcula_f(funcao_f,xk))
    return resp

#OK
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

#TESTADO OK
def PVC_1D(D,F):#,n,p,q,h,funcao_f,y0,yn,pontos_dominio):
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
    D_inv=inv(np.matrix(np.array(D)))#monta_matriz_D(n,p,q,h))))  #calcula a inversa
    F_t=np.transpose(F)#get_vetor_f(funcao_f,y0,yn,pontos_dominio,p,h,n))  #obtem o vetor vertical F
    #retorna o produto matricial
    return D_inv.dot(F_t)

#testado ok
def imp_mat(m):
    for linha in m:
        print(linha)

#testado OK
def retorna_f(funcao_em_lateX):
    """"retorna f(x) ; dada a string da representação de f em LaTeX"""
    # exemplo : 0.5* ( √(x+1) ) → em latex → \frac {\sqrt{1+x}} {2}
    return parse_latex(funcao_em_lateX)

#DEPRECIADO
def retorna_f_old(funcao_em_lateX):
    """"retorna f(x) ; dada a string da representação de f em LaTeX"""
    # exemplo : 0.5* ( √(x+1) ) → em latex → \frac {\sqrt{1+x}} {2}
    #FLX_no_ponto=
    FFX = latex2sympy(funcao_em_lateX, ConfigL2S())
    valor = parsing.sympy_parser.parse_expr(str(FFX))
    return valor
    #pass

#TEstado ok
def verifica_condicoes_problema(a,b,n,f_de_x):
    """Verifica se a entrada viola as condições do problema
    → a NÃO PODE SER IGUAL A b
    → n deve ser maior q 2
    → f deve fazer sentido em a e b
    → deve existir d²f/dx² (derivada segunda) para a função f
    """
    if a == b:
        print("[ERRO] Intervalo Vazio! a=b")
        exit(1)
    if n <= 2 :
        print("[ERRO] intervalo deve ser dividido pelo menos em 3 partes! n deve ser maior q 2!")
        exit(1)
    try:
        fA = calcula_f(f_de_x, a)
        if str(fA) == "zoo" or str(fA).find("*I") >= 0:
            print("[ERRO] a função f(x) =", f_de_x, 'não faz sentido no ponto (', a, ')')
            exit(1)
#        print(str(calcula_f(f_de_x,a)))
    except Exception as e:
        print("[ERRO] a função f(x) =",f_de_x,'não faz sentido no ponto (',a,')',e)
        exit(1)
    try:
        fA = calcula_f(f_de_x, b)
        if str(fA) == "zoo" or str(fA).find("*I") >= 0:
            print("[ERRO] a função f(x) =", f_de_x, 'não faz sentido no ponto (', b, ')')
            exit(1)
    except Exception as e:
        print("[ERRO] a função f(x) =", f_de_x, 'não faz sentido no ponto (', b, ')', e)
        exit(1)
    # CONDIÇÃO PRO PVC1D : Se existe a derivada segunda de f então existe a derivada quarta de y
    # logo vale o problema
    f = retorna_f(f_de_x)
    x = symbols("x")
    d2f_dx2 = None
    try:
        d2f_dx2 = diff(diff(f, x), x)
        if (not d2f_dx2) or (str(d2f_dx2) == 'zoo'):
            print("[ERRO] não existe d²f/dx² (derivada segunda) para a função f(x)=", f_de_x, "violando as condições do problema")
            exit(1)
        else:
            print("""[OK!] → """, a, """= a ≠ b=""", b, """
        → n > 2
        → existe f(a) e f(b), f(a)=""", calcula_f(f_de_x, a), """ f(b)=""", calcula_f(f_de_x, b), """
        → Existe d²f/dx² ; d²f/dx²=""", d2f_dx2,
    """ 
    → f(x) = """, f)
    except Exception as e:
        print("[ERRO] não existe d²f/dx² (derivada segunda) para a função f(x)=", f_de_x,
              "violando as condições do problema", e)
        exit(1)

def input_num(msg):
    t=input(msg)
    if t.lower()=="pi":
        return math.pi
    elif t.lower()=="e" or t.lower()=="exp":
        return math.e
    else :
        return float(t)


while(True):
    f_de_X = input("Insira a função f(x) em LaTeX: ")
    # \frac {( \sqrt{x-1})^3 } {x}

    inicio_dominio = input_num("Início de domínio (a): ")
    final_dominio = input_num("Final de domínio (b): ")
    qtd_sub_intervalos = int(input("Qtd. de sub-intervalos (n): "))
    p_de_X = input_num("Insira a função p(x) em LaTeX: ")
    q_de_X = input_num("Insira a função q(x) em LaTeX: ")
    y_0 = input_num("Valor de y(a): ")
    dy_dx_b = input("Valor de dy/dx (b): ")  ##ESTE VALOR É DADO PELO PROBLEMA ???
    h = (final_dominio - inicio_dominio) / qtd_sub_intervalos
    fun_F = retorna_f(f_de_X)
    fun_P = retorna_f(p_de_X)
    fun_Q = retorna_f(q_de_X)

    print("""Sistema em questão :
           { -""", fun_P, """dy/dx + """, fun_Q, """y =""", fun_F, """
           { y(""", inicio_dominio, """)=""", y_0, """
           { dy/dx (""", final_dominio, """)=""", dy_dx_b)

    print('\n Verificando condições do problema:')
    verifica_condicoes_problema(inicio_dominio, final_dominio, qtd_sub_intervalos, f_de_X)

    print("\nDISCRETIZAÇÃO:")
    dominio = discretiza_dominio1D(inicio_dominio, final_dominio, qtd_sub_intervalos)
    print(dominio)
    print("\n")

    print("CALCULANDO VETOR F:")
    vetor_F = get_vetor_f(f_de_X, y_0, dominio, p_de_X, h, qtd_sub_intervalos) # usar y_0 ao invés de a
    #for x_i in X:
    #    Y.append(calcula_f(f_de_X,x_i))
    print(vetor_F)
    print("\n")

    matriz = monta_matriz_D(qtd_sub_intervalos, p_de_X, q_de_X, h)
    print("MATRIZ:")
    imp_mat(matriz)
#    for linha in matriz:
#        print(linha)
    print("\n\n")

    print("SOLUCIONANDO PROBLEMA:")
    vetor_Y = PVC_1D(matriz, vetor_F)
    print('VETOR Y=', vetor_Y)

    print("-" * 50, '\n')


#
#latex2sympy()
