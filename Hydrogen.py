
"""Module with functions that calculate the energy of states of the
Hydrogen atom and plot their radial probability density functions"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp
import math



def myfunc(u, r, l, m, E, a):
    return (u[1], u[0]*((l*(l+1)/r**2)-2*m*(E+(a/r))))



def solve_it(r, l, m, E, a):

    """Solves the radial Schrodinger equation and returns the solution and its derivative in one array"""

    solve_me=np.array([0.,1.])
    return sp.odeint(myfunc, solve_me, r, args=( l, m, E, a, ),)



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
    return normSolArray



def plot_graph(p, normSolArray):

    """Plots a graph of the radial probability density function against radius"""

    y=normSolArray
    x=np.arange(1.0,p+1)*delta_r
    plt.figure()
    plt.plot(x,y)
    plt.xlabel('Radius (MeV^-1)')
    plt.ylabel('Radial probability density')
    plt.title('Radial probability density function vs. radius') 
    plt.show()



if __name__ == "__main__" :
    n=2
    l=1
    E1=-10**(-5)
    E3=-10**(-7)
    E2=(E1+E3)/2.0
    m=20
    delta_r=0.01
    r=np.arange(1.0,700000.0)*delta_r
    mass=0.510998928
    a=1.0/137.0
    for k in range(m):
        answer1=solve_it(r,0,mass,E1,a)
        answer2=solve_it(r,0,mass,E2,a)
        answer3=solve_it(r,0,mass,E3,a)
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
    print "The Energy of the n=2, l=1 state is", E2, "MeV"
    finalAnswer=solve_it(r,l,mass,E2,a)
    finalData=find_number_of_nodes_and_turning_points(finalAnswer)
    point=find_cutoff_point(finalAnswer, finalData[0], finalData[1], finalData[2], finalData[3])
    normAnswer=normalise_solution(finalAnswer, point, delta_r)
    plot_graph(point, normAnswer)
