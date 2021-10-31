#!/bin/python3
import matplotlib.pyplot as plt
import numpy as np
def standart_values():
    N1 = 20 # кол-во промежутков на пластину
    N2 = 20
    N3 = 20
    t_end = 160 # окончание по времени

    lambda1 = 48 # теплопроводности матереала А
    lambda2 = 384 # теплопроводности матереала B
    lambda3 = 50

    ro1 = 7800 # плотность матереала А
    ro2 = 8800 # плотность матереала B
    ro3 = 6000

    c1 = 460 # теплоемкость матереала А
    c2 = 381 # теплоемкость матереала В
    c3 = 300

    tl = 42.5 # температура на границе х = 0
    tr = 22 # температура на границе х = L 
    t0 = 20 # температура начальная

    l = 0.3 # толщина пластина 1m = 0.1
    arr_values=[N1,N2,N3,t_end,lambda1,lambda2,lambda3,
            ro1,ro2,ro3,c1,c2,c3,tl,tr,t0,l]
    return  arr_values

def calculation(arr_values):
    N1,N2,N3,t_end,lambda1,lambda2,lambda3, ro1,ro2,ro3,c1,c2,c3,tl,tr,t0,l = arr_values
    N1 = int(N1)
    N2 = int(N2)
    N3 = int(N3)
    # общее число узов:
    N = N1 + N2 + N3 + 1
    # N = int(N)
    # расчетный шаг сетки по пространственной координате:
    h = l / (N - 1)

    # определяем коэфициент температуропроводности:
    a1 = lambda1/(ro1 * c1)
    a2 = lambda2/(ro2 * c2)
    a3 = lambda3/(ro3 * c3)
    # определяем расчетный шаг сетки по времени:
    tau = t_end / 100.0

    # определяем поле температуры в начальный момент времени:
    t = [None] * N
    for i in range(N):
        t[i] = t0

    # проводим интегрирование нестационарного уравнения теплопроводности
    time = 0
    alpha = [None] * N
    beta = [None] * N
    while time < t_end:
        time = time + tau

        # определяем начальные прогоночные коэфициэнты на основе левого граничного условия
        alpha[0] = 0.0
        beta[0] = tl
        
        # цикл с параметром для определения прогоночных коэфициентов по формуле 8 в первой части пластины
        for i in range(1, N1):
            # ai, bi, ci, fi коэфициент канонического представления СЛАУ с трехдиагональной матрицей
            ai = lambda1 / pow(h, 2)
            bi = 2.0 * lambda1 / pow(h, 2) + (ro1 * c1 / tau)
            ci = lambda1 / pow(h, 2)
            fi = (-ro1) * c1 * t[i] / tau

            # alpha[i], beta[i] – прогоночные коэффициенты
            alpha[i] = ai / (bi - ci * alpha[i-1])
            beta[i] = (ci * beta[i-1]-fi) / (bi -ci * alpha[i-1])

        #  определяем прогоночные коэффициенты на границе раздела первой и второй
        # частей, используем соотношения (28)

        alpha[N1]=2.0*a1*a2*tau*lambda2/(2.0*a1*a2*tau*(lambda2+lambda1 
        *(1-alpha[N1-1]))+pow(h,2)*(a1*lambda2+a2*lambda1))
        beta[N1]=(2.0*a1*a2*tau*lambda1*beta[N1-1]+pow(h,2)*(a1*lambda2+a2 
        *lambda1)*t[N1])/(2.0*a1*a2*tau*(lambda2+lambda1 
        *(1-alpha[N1-1]))+pow(h,2)*(a1*lambda2+a2*lambda1))
        
        # цикл с параметром для определения прогоночных коэффициентов по
        # формуле (8) во второй части пластины
        for i in range(N1+1,N-1):
            # {ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с
            # трехдиагональной матрицей
            ai = lambda2/pow(h,2);
            bi = 2.0*lambda2/pow(h,2)+ro2*c2/tau;
            ci = lambda2/pow(h,2);
            fi = -ro2*c2*t[i]/tau;
            # alpha[i], beta[i] – прогоночные коэффициенты
            alpha[i] = ai/(bi-ci*alpha[i-1]);
            beta[i] = (ci*beta[i-1]-fi)/(bi-ci*alpha[i-1]);

        #  определяем прогоночные коэффициенты на границе раздела второй ит третей
        # частей, используем соотношения (28)
        alpha[N1*2]=2.0*a2*a3*tau*lambda3/(2.0*a2*a3*tau*(lambda3+lambda2 
        *(1-alpha[N1*2-1]))+pow(h,2)*(a2*lambda3+a3*lambda2))
        beta[N1*2]=(2.0*a2*a3*tau*lambda2*beta[N1*2-1]+pow(h,2)*(a2*lambda3+a3 
        *lambda2)*t[N1*2])/(2.0*a2*a3*tau*(lambda3+lambda2 
        *(1-alpha[N1*2-1]))+pow(h,2)*(a2*lambda3+a3*lambda2))
        
        # цикл с параметром для определения прогоночных коэффициентов по
        # формуле (8) во третей части пластины
        for i in range(N1*2+1,N-1):
            # {ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с
            # трехдиагональной матрицей
            ai = lambda3/pow(h,2);
            bi = 2.0*lambda3/pow(h,2)+ro3*c3/tau;
            ci = lambda3/pow(h,2);
            fi = -ro3*c3*t[i]/tau;
            # alpha[i], beta[i] – прогоночные коэффициенты
            alpha[i] = ai/(bi-ci*alpha[i-1]);
            beta[i] = (ci*beta[i-1]-fi)/(bi-ci*alpha[i-1]);




        # определяем значение температуры на правой границе
        t[N-1] = tr

        # используя соотношение (7) определяем неизвестное поле
        # температуры
        for i in range(N-2,-1,-1):
            t[i] = alpha[i] * t[i + 1] + beta[i]


    plt.plot(t)
    plt.show()

def overheat():
    N = 10000
    t_end = 160 # окончание по времени

    lambda1 = 48 # теплопроводности матереала А
    ro1 = 7800 # плотность матереала А
    c1 = 460 # теплоемкость матереала А
    tl = 1900 # температура на границе х = 0
    tr = 5000 # температура на границе х = L 
    t0 = 1900 # температура начальная

    l = 0.1 # толщина пластина 1 m = 0.1
    a1 = lambda1/(ro1 * c1)
    h = l / (N - 1)
    tau = t_end / 100.0

    # определяем поле температуры в начальный момент времени:
    t = [None] * N
    for i in range(N):
        t[i] = t0
    time = 0
    alpha = [None] * N
    beta = [None] * N
    while time < t_end:
        time = time + tau

        # определяем начальные прогоночные коэфициэнты на основе левого граничного условия
        alpha[0] = 0.0
        beta[0] = tl
        for i in range(1, N):
            # ai, bi, ci, fi коэфициент канонического представления СЛАУ с трехдиагональной матрицей
            ai = lambda1 / pow(h, 2)
            bi = 2.0 * lambda1 / pow(h, 2) + (ro1 * c1 / tau)
            ci = lambda1 / pow(h, 2)
            fi = (-ro1) * c1 * t[i] / tau

            # alpha[i],  прогоночные коэффициенты
            alpha[i] = ai / (bi - ci * alpha[i-1])
            beta[i] = (ci*beta[i-1]-fi)/(bi-ci*alpha[i-1]);

        t[N-1] = tr
        for i in range(N-2,-1,-1):

            t[i] = alpha[i] * t[i + 1] + beta[i]
    for i in range(N-1,-1,-1):
        if t[i] > 3000:
            t[i] = 3000
        # else:
            # t
    plt.plot(t)
    plt.show()



if __name__ == '__main__':
    print("take input from file [yes/y] or from program [no/n]?")
    input_ = input()
    # соблюдайте пробелы в файле!
    if input_ == 'y' or input_ =='yes':
        # print("eeee")
        arr_values = []
        with open('input.txt','r') as file:
            # for line in file:
            arr_all_text = file.read().split("\n")
            for i in arr_all_text:
                arr_values.append(float(i.split(" ")[2]))
        calculation(arr_values)
    elif input_ == 'n' or input_ == 'no':
        # print("nnnn")
        calculation(standart_values())
    overheat()
    