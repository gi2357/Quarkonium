
"""Module with functions that calculate the energy of states of quarkonium systems and plot their radial probability density functions"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp
import math


def myfunc1(u, r, l, m, E, a, b):
    return (u[1], u[0]*((l*(l+1)/(r**2))-2*m*(E-b*r+(4/3)*a/r)))


def solve_it1(r, l, m, E, a, b):

    """Solves the radial Schrodinger equation and returns the solution and its derivative in one array"""

    solve_me=np.array([0.,1.])
    return sp.odeint(myfunc1, solve_me, r, args=( l, m, E, a, b, ), )


def myfunc2(u, r, l, m, E, c1, d, c2):
    return (u[1], u[0]*((l*(l+1)/(r**2))-2*m*(E-(c1*(r**d)-c2))))


def solve_it2(r, l, m, E, c1, d, c2):

    """Solves the radial Schrodinger equation and returns the solution and its derivative in one array"""

    solve_me=np.array([0.,1.])
    return sp.odeint(myfunc2, solve_me, r, args=( l, m, E, c1, d, c2, ), )


def myfunc3(u, r, l, m, E, a, b):
    return (u[1], u[0]*((l*(l+1)/(r**2))-2*m*(E+b*(r**(-1/2))-a*(r**(1/2)))))


def solve_it3(r, l, m, E, a, b):

    """Solves the radial Schrodinger equation and returns the solution and its derivative in one array"""

    solve_me=np.array([0.,1.])
    return sp.odeint(myfunc3, solve_me, r, args=( l, m, E, a, b, ), )


def find_number_of_nodes_and_turning_points(solArray):

    """Finds the positions and the number of nodes and turning points of the solution and returns them in a tuple"""

    length1=len(solArray[:,0])-1
    nArray=np.zeros(length1+1)  
    tpArray=np.zeros(length1+1)
    for i in range(length1):
        if math.copysign(1,solArray[i,0])!=math.copysign(1,solArray[i+1,0]):  #compare the signs of the ith and the (i+1)th element 
            nArray[i]=1                                                       #store the positions of nodes
        if math.copysign(1,solArray[i,1])!=math.copysign(1,solArray[i+1,1]):  #compare the signs of the ith and the (i+1)th element  
            tpArray[i]=1                                                      #store the positions of turning points
    nodes=sum(nArray)  #find number of nodes
    t_points=sum(tpArray)  #find number of turning points
    return nodes, t_points, nArray, tpArray
    


def find_cutoff_point(solArray, nodes, t_points, nArray, tpArray):

    """Finds the position of the last node or turning point of the solution"""

    length2=len(nArray)-1
    ntpArray=np.zeros(length2+1) 
    for j in range(length2+1):
        ntpArray[j]=nArray[j]+tpArray[j]  #put the positions of nodes and turning points in a single array
    for k in range(length2):
        if ntpArray[length2-k]!=0:        #iterate from the end of the array to find
            p=length2-k                   #the position of the last node or turning point
            break
    return p


def normalise_solution(solArray, p, delta_r):

    """Calculates the normalisation constant for the solution and returns an array with the normalised solution"""

    const=sum(solArray[0:p,0]**2)*delta_r        #calculate normalisation constant
    normSolArray=(1/const)*(solArray[0:p,0]**2)  #find normalised solution of the radial probability density function
    return normSolArray, const


def plot_graph(p, normSolArray, delta_r):

    """Plots a graph of the radial probability density function against radius"""

    y=normSolArray
    x=np.arange(1.0,p+1)*delta_r
    plt.figure()
    plt.plot(x,y, color='k')
    plt.xlabel('Radius (GeV^-1)')
    plt.ylabel('Radial probability density')
    plt.show()


def calculate_b():
    l=0
    E=0.368
    mass=1.35/2
    a=0.40
    b1=0.05
    b3=0.5
    b2=(b1+b3)/2.0
    delta_r=0.001
    r=np.arange(1.0,25000.0)*delta_r
    for k in range(20):
        answer1=solve_it1(r,l,mass,E,a,b1)
        answer2=solve_it1(r,l,mass,E,a,b2)
        answer3=solve_it1(r,l,mass,E,a,b3)
        data1=find_number_of_nodes_and_turning_points(answer1)
        data2=find_number_of_nodes_and_turning_points(answer2)
        data3=find_number_of_nodes_and_turning_points(answer3)
        if data1[1]!=data2[1] and data1[0]!=data2[0]:    #iterate towards solution by comparing the number of
            b3=b2                                        #of nodes and turning points of the test values of b
        elif data3[0]!=data2[0] and data3[1]!=data2[1]:  
            b1=b2                                        #data[0]=nodes and data[1]=t_points
        else:
            pass
        b2=(b1+b3)/2.0
    print "b=", b2


def charmonium_energy_and_wavefunction(n,l):
    if n==1:
        if l==0:
            E1=0.1
            E3=0.5
        if l==1:
            E1=0.5
            E3=0.8
        if l==2:
            E1=0.8
            E3=1.1
    if n==2:
        if l==0:
            E1=0.5
            E3=0.9
        if l==1:
            E1=0.9
            E3=1.2
    E2=(E1+E3)/2.0
    mass=1.35/2
    a=0.40
    b=0.145
    delta_r=0.001
    r=np.arange(1.0,25000.0)*delta_r
    for k in range(20):
        answer1=solve_it1(r,l,mass,E1,a,b)
        answer2=solve_it1(r,l,mass,E2,a,b)
        answer3=solve_it1(r,l,mass,E3,a,b)
        data1=find_number_of_nodes_and_turning_points(answer1)
        data2=find_number_of_nodes_and_turning_points(answer2)
        data3=find_number_of_nodes_and_turning_points(answer3)
        if data1[1]!=data2[1] and data1[0]!=data2[0]:    #iterate towards solution by comparing the number of
            E3=E2                                        #of nodes and turning points of the test energies
        elif data3[0]!=data2[0] and data3[1]!=data2[1]:  
            E1=E2                                        #data[0]=nodes and data[1]=t_points
        else:
            pass
        E2=(E1+E3)/2.0
    print "The Energy of the n=",n,", l=",l, "state is", E2, "GeV"
    point=find_cutoff_point(answer2, data2[0], data2[1], data2[2], data2[3])
    normAnswer=normalise_solution(answer2, point, delta_r)[0]
    plot_graph(point, normAnswer, delta_r)


def bottomonium_energy_and_wavefunction(n,l):
    if n==1:
        if l==0:
            E1=0
            E3=0.3
        if l==1:
            E1=0.3
            E3=0.6
        if l==2:
            E1=0.6
            E3=1
    if n==2:
        if l==0:
            E1=0.3
            E3=0.6
        if l==1:
            E1=0.6
            E3=1
    if n==3:
        if l==0:
            E1=0.6
            E3=0.9
        if l==1:
            E1=0.9
            E3=1.2
    if n==4:
        E1=0.9
        E3=1.2
    E2=(E1+E3)/2.0
    mass=4.70/2
    a=0.28
    b=0.145
    delta_r=0.001
    r=np.arange(1.0,25000.0)*delta_r
    for k in range(20):
        answer1=solve_it1(r,l,mass,E1,a,b)
        answer2=solve_it1(r,l,mass,E2,a,b)
        answer3=solve_it1(r,l,mass,E3,a,b)
        data1=find_number_of_nodes_and_turning_points(answer1)
        data2=find_number_of_nodes_and_turning_points(answer2)
        data3=find_number_of_nodes_and_turning_points(answer3)
        if data1[1]!=data2[1] and data1[0]!=data2[0]:    #iterate towards solution by comparing the number of
            E3=E2                                        #of nodes and turning points of the test energies
        elif data3[0]!=data2[0] and data3[1]!=data2[1]:  
            E1=E2                                        #data[0]=nodes and data[1]=t_points
        else:
            pass
        E2=(E1+E3)/2.0
    print "The Energy of the n=",n,", l=",l, "state is", E2, "GeV"
    point=find_cutoff_point(answer2, data2[0], data2[1], data2[2], data2[3])
    normAnswer=normalise_solution(answer2, point, delta_r)[0]
    plot_graph(point, normAnswer, delta_r)


def hyperfine_splitting_of_charmonium(n):
    if n==1:
        E1=0
        E3=0.4
    if n==2:
        E1=0.4
        E3=0.8
    E2=(E1+E3)/2.0
    mass=1.35/2
    a=0.40
    b=0.145
    delta_r=0.001
    r=np.arange(1.0,25000.0)*delta_r
    for k in range(20):
        answer1=solve_it1(r,0,mass,E1,a,b)
        answer2=solve_it1(r,0,mass,E2,a,b)
        answer3=solve_it1(r,0,mass,E3,a,b)
        data1=find_number_of_nodes_and_turning_points(answer1)
        data2=find_number_of_nodes_and_turning_points(answer2)
        data3=find_number_of_nodes_and_turning_points(answer3)
        if data1[1]!=data2[1] and data1[0]!=data2[0]:    #iterate towards solution by comparing the number of
            E3=E2                                        #of nodes and turning points of the test energies
        elif data3[0]!=data2[0] and data3[1]!=data2[1]:  
            E1=E2                                        #data[0]=nodes and data[1]=t_points
        else:
            pass
        E2=(E1+E3)/2.0
    point=find_cutoff_point(answer2, data2[0], data2[1], data2[2], data2[3])
    const=normalise_solution(answer2, point, delta_r)[1]
    hypSplit=(8.0/9)*(a/((2*mass)**2))*(2/const**2)
    print "The hyperfine splitting of the n=", n, ",l= 0 state is", hypSplit, "GeV"


def hyperfine_splitting_of_bottomonium(n):
    if n==1:
        E1=0
        E3=0.3
    if n==2:
        E1=0.3
        E3=0.6
    E2=(E1+E3)/2.0
    mass=4.7/2
    a=0.28
    b=0.145
    delta_r=0.001
    r=np.arange(1.0,25000.0)*delta_r
    for k in range(20):
        answer1=solve_it1(r,0,mass,E1,a,b)
        answer2=solve_it1(r,0,mass,E2,a,b)
        answer3=solve_it1(r,0,mass,E3,a,b)
        data1=find_number_of_nodes_and_turning_points(answer1)
        data2=find_number_of_nodes_and_turning_points(answer2)
        data3=find_number_of_nodes_and_turning_points(answer3)
        if data1[1]!=data2[1] and data1[0]!=data2[0]:    #iterate towards solution by comparing the number of
            E3=E2                                        #of nodes and turning points of the test energies
        elif data3[0]!=data2[0] and data3[1]!=data2[1]:  
            E1=E2                                        #data[0]=nodes and data[1]=t_points
        else:
            pass
        E2=(E1+E3)/2.0
    point=find_cutoff_point(answer2, data2[0], data2[1], data2[2], data2[3])
    const=normalise_solution(answer2, point, delta_r)[1]
    hypSplit=(4.0/9)*(a/((2*mass)**2))*(1/(2*const**2))
    print "The hyperfine splitting of the n=", n, ",l= 0 state is", hypSplit, "GeV"


def charmonium_energy_and_wavefunction_potential_2(n,l):
    if n==1:
        if l==0:
            E1=0.1
            E3=0.5
        if l==1:
            E1=0.5
            E3=0.8
        if l==2:
            E1=0.8
            E3=1.2
    if n==2:
        if l==0:
            E1=0.5
            E3=1
        if l==1:
            E1=1
            E3=1.2
    E2=(E1+E3)/2.0
    mass=1.35/2
    c1=6.87
    d=0.1
    c2=7.3
    delta_r=0.001
    r=np.arange(1.0,25000.0)*delta_r
    for k in range(20):
        answer1=solve_it2(r,l,mass,E1,c1,d,c2)
        answer2=solve_it2(r,l,mass,E2,c1,d,c2)
        answer3=solve_it2(r,l,mass,E3,c1,d,c2)
        data1=find_number_of_nodes_and_turning_points(answer1)
        data2=find_number_of_nodes_and_turning_points(answer2)
        data3=find_number_of_nodes_and_turning_points(answer3)
        if data1[1]!=data2[1] and data1[0]!=data2[0]:    #iterate towards solution by comparing the number of
            E3=E2                                        #of nodes and turning points of the test energies
        elif data3[0]!=data2[0] and data3[1]!=data2[1]:  
            E1=E2                                        #data[0]=nodes and data[1]=t_points
        else:
            pass
        E2=(E1+E3)/2.0
    print "The Energy of the n=",n,", l=",l, "state is", E2, "GeV"
    point=find_cutoff_point(answer2, data2[0], data2[1], data2[2], data2[3])
    normAnswer=normalise_solution(answer2, point, delta_r)[0]
    plot_graph(point, normAnswer, delta_r)


def charmonium_energy_and_wavefunction_potential_3(n,l):
    if n==1:
        if l==0:
            E1=0
            E3=0.3
        if l==1:
            E1=0.3
            E3=0.5
        if l==2:
            E1=0.45
            E3=0.6
    if n==2:
        if l==0:
            E1=0.3
            E3=0.7
        if l==1:
            E1=0.45
            E3=0.9
    E2=(E1+E3)/2.0
    mass=1.35/2
    a=0.51
    b=0.924
    delta_r=0.001
    r=np.arange(1.0,60000.0)*delta_r
    for k in range(20):
        answer1=solve_it3(r,l,mass,E1,a,b)
        answer2=solve_it3(r,l,mass,E2,a,b)
        answer3=solve_it3(r,l,mass,E3,a,b)
        data1=find_number_of_nodes_and_turning_points(answer1)
        data2=find_number_of_nodes_and_turning_points(answer2)
        data3=find_number_of_nodes_and_turning_points(answer3)
        if data1[1]!=data2[1] and data1[0]!=data2[0]:    #iterate towards solution by comparing the number of
            E3=E2                                        #of nodes and turning points of the test energies
        elif data3[0]!=data2[0] and data3[1]!=data2[1]:  
            E1=E2                                        #data[0]=nodes and data[1]=t_points
        else:
            pass
        E2=(E1+E3)/2.0
    print "The Energy of the n=",n,", l=",l, "state is", E2, "GeV"
    point=find_cutoff_point(answer2, data2[0], data2[1], data2[2], data2[3])
    normAnswer=normalise_solution(answer2, point, delta_r)[0]
    plot_graph(point, normAnswer, delta_r)
    
