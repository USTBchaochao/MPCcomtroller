import time
import numpy
from math import *
from cvxopt import matrix , solvers

def MPCController(State , Speed , Ts , Np , Nc , dumax , dumin , dsmax , dsmin , dmax , dmin , ROW , L , Route):
    start = time.time()
    Ns = numpy.shape(State)[0]#状态维数
    Nu = numpy.shape(Speed)[0]#控制量维数

    M = numpy.zeros((2*(Nu*(Nc+1)+1) ,(Nu*(Nc+1)+1)))#边界条件系数矩阵
    for i in range(Nu*(Nc+1)+1):
        for j in range(Nu*(Nc+1)+1):
            if i ==j:
                M[i,j] = 1
                M[i+(Nu*(Nc+1)+1) , j] = -1

    N = numpy.zeros((2*(Nu*(Nc+1)+1) , 1))
    for i in range(Nu*(Nc+1)):
        if i%2 == 0:
            N[i,0] = dumax
        else:
            N[i,0] = dsmax
        N[i+(Nu*(Nc+1)+1),0] = N[i,0]
    N[(Nu*(Nc+1)+1)-1 , 0] = dmax
    N[2*(Nu*(Nc+1)+1)-1 , 0] = dmin

    DSpeed = numpy.array(([0],[0]))

    PosX = State[0,0]
    PosY = State[1,0]
    Heading = State[2,0]
    speed = Speed[0,0]
    steer = Speed[1,0]
    stateref = Route[i]
    Xref = stateref[0]
    Yref = stateref[1]
    Headingref = stateref[2]
    Stateref = numpy.array(([Xref],[Yref],[Headingref]))
    StateError = State - Stateref

    Ki = numpy.vstack((StateError , DSpeed))

    a = numpy.array(([1 ,0 , -Ts*speed*sin(Heading)]
                     ,[0 , 1 , Ts*speed*cos(Heading)]
                     ,[0 , 0 , 1]))
    b = numpy.array(([Ts*cos(Heading) , 0]
                    ,[Ts*sin(Heading) , 0]
                     ,[Ts*tan(Heading)/L , Ts*speed/(L*cos(steer)*cos(steer))]))
    A1 = numpy.hstack((a , b))
    A2 = numpy.hstack((numpy.zeros((Nu , Ns)) , numpy.eye(Nu)))
    A = numpy.vstack((A1 , A2))

    B = numpy.vstack((b , numpy.eye(Nu)))

    C = numpy.hstack((numpy.eye(Ns) , numpy.zeros((Ns , Nu))))

    PHI = C.dot(A)
    for i in range(2,Np+1):
        phi = (C.dot(A)).dot(A)
        PHI = numpy.vstack((PHI , phi))
        phi = phi.dot(A)

    list = []
    for i in range(1,Np+1):
        for j in range(1,Nc+1+1):
            if i>j:
                theta = C
                k = i-j
                for n in range(k):
                    theta = theta.dot(A)
                theta = theta.dot(B)
            elif i == j:
                theta = C.dot(B)
            else:
                theta = numpy.zeros((Ns , Nu))
            list.append(theta)
    THETA = numpy.array((list))
    THETAline1 = numpy.hstack((THETA[0] , THETA[1] , THETA[2] , THETA[3]))
    THETAline2 = numpy.hstack((THETA[4] , THETA[5] , THETA[6] , THETA[7]))
    THETAline3 = numpy.hstack((THETA[8] , THETA[9] , THETA[10] , THETA[11]))
    THETAline4 = numpy.hstack((THETA[12] , THETA[13] , THETA[14] , THETA[15]))
    THETA = numpy.vstack((THETAline1 , THETAline2 , THETAline3 , THETAline4))

    list = []
    for i in range(1,Np+1):
        for j in range(1,Np+1):
            if i == j:
                q = numpy.eye(Ns)
            else:
                q = numpy.zeros((Ns , Ns))
            list.append(q)
    Q = numpy.array(list)
    Qline1 = numpy.hstack((Q[0] , Q[1] , Q[2] , Q[3]))
    Qline2 = numpy.hstack((Q[4] , Q[5] , Q[6] , Q[7]))
    Qline3 = numpy.hstack((Q[8] , Q[9] , Q[10] , Q[11]))
    Qline4 = numpy.hstack((Q[12] , Q[13] , Q[14] , Q[15]))
    Q = numpy.vstack((Qline1 , Qline2 , Qline3 , Qline4))

    list = []
    for i in range(1, Nc+1+1):
        for j in range(1 , Nc+1+1):
            if i == j:
                r = numpy.eye(Nu)
            else:
                r = numpy.zeros((Nu , Nu))
            list.append(r)
    R = numpy.array(list)
    Rline1 = numpy.hstack((R[0] , R[1] , R[2] , R[3]))
    Rline2 = numpy.hstack((R[4] , R[5] , R[6] , R[7]))
    Rline3 = numpy.hstack((R[8] , R[9] , R[10] , R[11]))
    Rline4 = numpy.hstack((R[12] , R[13] , R[14] , R[15]))
    R = numpy.vstack((Rline1 , Rline2 , Rline3 , Rline4))

    E = PHI.dot(Ki)
    H1 = numpy.hstack((numpy.transpose(THETA).dot(Q).dot(THETA)+R , numpy.zeros((Nu*(Nc+1) , 1))))
    # H2 = numpy.hstack((numpy.zeros((1 , (Nu*(Nc+1)))) , ROW))
    H2 = numpy.zeros((1 , (Nu*(Nc+1)+1)))
    H2[0,(Nu*(Nc+1))] = ROW
    H = numpy.vstack((H1 , H2))

    G = numpy.hstack((2*numpy.transpose(E).dot(Q).dot(THETA) , numpy.zeros((1,1))))

    H = numpy.transpose(H)
    G = numpy.transpose(G)
    # print(numpy.shape(G))
    H = matrix(H)
    G = matrix(G)
    M = matrix(M)
    N = matrix(N)
    result = solvers.qp(H , G , M , N , kktsolver=None)
    end = time.time()
    print('time cost:' , end-start)
    x0 = numpy.array(result['x'])
    dspeed = x0[0][0]
    dsteer = x0[1][0]
    DSpeed[0][0] = dspeed
    DSpeed[1][0] = dsteer

    return  dspeed , dsteer
