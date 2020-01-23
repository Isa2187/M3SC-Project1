import numpy as np
import scipy as sp
import math as ma
import sys
import time




def Dijkst(ist,isp,wei):
    # Dijkstra algorithm for shortest path in a graph
    #    ist: index of starting node
    #    isp: index of stopping node
    #    wei: weight matrix

    # exception handling (start = stop)
    if (ist == isp):
        shpath = [ist]
        return shpath

    # initialization
    N         =  len(wei)
    Inf       =  sys.maxint
    UnVisited =  np.ones(N,int)
    cost      =  np.ones(N)*1.e6
    par       = -np.ones(N,int)*Inf

    # set the source point and get its (unvisited) neighbors
    jj            = ist
    cost[jj]      = 0
    UnVisited[jj] = 0
    tmp           = UnVisited*wei[jj,:]
    ineigh        = np.array(tmp.nonzero()).flatten()
    L             = np.array(UnVisited.nonzero()).flatten().size

    # start Dijkstra algorithm
    while (L != 0):
        # step 1: update cost of unvisited neighbors,
        #         compare and (maybe) update
        for k in ineigh:
            newcost = cost[jj] + wei[jj,k]
            if ( newcost < cost[k] ):
                cost[k] = newcost
                par[k]  = jj

        # step 2: determine minimum-cost point among UnVisited
        #         vertices and make this point the new point
        icnsdr     = np.array(UnVisited.nonzero()).flatten()
        cmin,icmin = cost[icnsdr].min(0),cost[icnsdr].argmin(0)
        jj         = icnsdr[icmin]

        # step 3: update "visited"-status and determine neighbors of new point
        UnVisited[jj] = 0
        tmp           = UnVisited*wei[jj,:]
        ineigh        = np.array(tmp.nonzero()).flatten()
        L             = np.array(UnVisited.nonzero()).flatten().size

    # determine the shortest path
    shpath = [isp]
    while par[isp] != ist:
        shpath.append(par[isp])
        isp = par[isp]
    shpath.append(ist)

    return shpath[::-1]




def calcWei(RX,RY,RA,RB,RV):
    # calculate the weight matrix between the points

    n    = len(RX)
    wei = np.zeros((n,n),dtype=float)
    m    = len(RA)
    for i in range(m):
        xa = RX[RA[i]-1]
        ya = RY[RA[i]-1]
        xb = RX[RB[i]-1]
        yb = RY[RB[i]-1]
        dd = ma.sqrt((xb-xa)**2 + (yb-ya)**2)
        tt = dd/RV[i]
        wei[RA[i]-1,RB[i]-1] = tt
    return wei





#Define a function where we can model the traffic flow easily 
#for different parameter values
def trafficFlow(ist,iend,wei0,minutes,numInject,freqInject,e):
    
    #ist = index of the starting node
    #iend = index of the end node
    #wei0 = weight matrix at start of the traffic flow
    #minutes = number of iterations we simulate
    #numInject = number of cars to inject at the end of
    #            each iteration we are adding cars at
    #freqInject = number of times we inject new cars into 
    #             the network
    #e = parameter used when updating the weight matrix
    
    #Initialise the weight matrix as the original one, but
    #note that wei will be updated at the end of each
    #iteration after the cars have moved through the network
    wei = wei0.copy()
        
    #Initialise the matrix carFreq which stores the number of
    #cars at each node at the end of each iteration once
    #the cars have moved and we have added more cars to the
    #network if necessary. 
    
    #Note: the first row i.e. carFreq[0,:] is a row of zeros
    #because I start off the simulation with no cars in the
    #network.
    
    #Note: numNodes is the number of nodes in the network, and
    #can be determined by looking at the size of the weight
    #matrix
    numNodes = int(np.shape(wei0)[0])
    carFreq = np.zeros((minutes,numNodes),dtype=int)
    
    #Initialise an empty list which will update and store the
    #edges that we use as we simulate the traffic flow
    usedEdges = []
    
    #Complete the iterations using a for loop, iterating
    #over t
    for t in range(minutes):
            
        #After the first iteration is complete i.e. once we
        #have cars in the network, update the current row of
        #zeros by the previous row so we know how many cars are
        #at each node at the start of the current iteration
        
        #Note: we start with no cars in the network, so at 
        #the start of the first iteration with python-index
        #t=0 we have that carFreq[0,:] is a row of zeros
        if t > 0:
            carFreq[t,:] = carFreq[t-1,:]
        
        #Make a copy of how many cars are at each node at the
        #start of the current iteration so we can refer back
        #to it when computing how many cars (if any) need to 
        #move from each node
        tempVec = carFreq[t,:].copy()
        
        #Loop through each of the nodes and determine the next
        #node that the cars should move to using the path 
        #given by Dijkstra's algorithm, where we use the 
        #current weight matrix to find the current optimal
        #Dijkstra path
        for i in range(numNodes):
        
            #Only apply DIjkstra's algorithm to nodes where 
            #there are currently cars using the indexes of 
            #these nodes and the current weight matrix, 
            #because otherwise we will compute a path but
            #won't have any cars to move
            if tempVec[i] != 0:
                
                #Store the Dijkstra path in the list shpath so
                #that we can obtain the next node to move the 
                #cars to later on
                shpath = Dijkst(i,iend,wei)
                
                #Assign numCars to be the current number of  
                #cars at the node we are currently interested 
                #in, so we can refer to it later when 
                #computing how many cars to move
                numCars = tempVec[i]
                
                #Move cars through the network if we are not 
                #at the end of the path i.e. not at node with
                #python-index iend, so that we do not lose
                #cars incorrectly
                if i != iend:
                    
                    #If we aren't at the end node, compute 70%
                    #of the number of cars currently there to 
                    #the nearest integer so these are the cars 
                    #to be moved. Note that we only compute  
                    #this value once and we add it to the next
                    #node and subtract it from the current 
                    #node in order to conserve the number of 
                    #cars in the network, so there is no net 
                    #gain nor loss of cars in the network
                    numMove = int(round(numCars*0.7))
                                        
                    #Obtain the python-index of the next node 
                    #for the cars to move to from the Dijkstra
                    #path
                    nextNode = shpath[1]
                                        
                    #Add the cars to be moved to the number of 
                    #cars at the next node
                    carFreq[t,nextNode] = carFreq[t,nextNode] + numMove
                                        
                    #Remove these cars from the node that they
                    #were at when the iteration started
                    carFreq[t,i] = carFreq[t,i] - numMove
                    
                    #Add the edge used by the cars moving from
                    #node i+1 (python-index i) to node 
                    #nextNode+1 (python-index nextNode) to the 
                    #list of used edges using append
                    usedEdges.append([i+1,nextNode+1])
                
                #If we are at the end node, then compute 40% 
                #of the cars that are currently there to the 
                #nearest integer and remove them from the  
                #network by subtracting them from the end node
                else:                    
                    numMove = int(round(numCars*0.4))
                    carFreq[t,i] = carFreq[t,i] - numMove
                
        #If we haven't completed the first freqInject
        #iterations, add numInject cars into the network at  
        #the starting node with python-index ist in 
        #conjunction with the function argument
        if t <= (freqInject-1):
            carFreq[t,ist] = carFreq[t,ist] + numInject
        
        
        
        #Update the weight matrix wei for the elements that 
        #are non-zero, since we only update the weights 
        #corresponding to places where there are actually 
        #links between the nodes. We do this using the formula 
        #given to us, so we add the original weight to e 
        #(argument of the function) multiplied by half of the 
        #sum of the new number of cars which are there at the 
        #end of the current iteration
        for i in range(numNodes):
            for j in range(numNodes):
                if wei[i,j] != 0:
                    wei[i,j] = wei0[i,j] + e*(carFreq[t,i] + carFreq[t,j])/2
    
    #Note: at this point, we have completed the simulation of
    #the model
    
    #Find the maximum number of cars at each node by finding
    #the maximum value in each column, since each column
    #contains the number of cars at each node at the end of 
    #each iteration
    maximumLoad = np.amax(carFreq, axis = 0)
    
    #Find the 5 most congested nodes from the array 
    #maximumLoad, and list them from the most congested in 
    #descending order in the array mostCongested.
    #Note: I need to add 1 to the indexes since python-indexes 
    #start from 0 and not 1
    indexCongested = np.argsort(maximumLoad)[-5:][::-1]
    mostCongested = indexCongested + 1
        
    #Find how many cars are at each of the 5 most congested 
    #nodes using the indexes of these nodes, and store these
    #in the array loadCongested
    loadCongested = maximumLoad[indexCongested]
        
    #Find the unvisited nodes by finding where the maximum 
    #load of cars over all iterations was 0 in maximumLoad,
    #meaning there were never any cars there at any point 
    #in the flow model. Again, I need to add 1 to remove the
    #effect of python-indexes 
    unvisitedNodes = np.where(maximumLoad == 0)[0] + 1
  
    #Make a list of the used edges called usedUnique, which
    #removes the repeats of the used edges in usedEdges
    usedUnique = []
    [usedUnique.append(edge) for edge in usedEdges if edge not in usedUnique]
    #Using the edge list, which is a list of lists called
    #edgeList and will be computed from the files RomeVertices
    #and RomeEdges before using this function trafficFlow, 
    #take out the edges which have been used so that we can 
    #obtain the unused edges. To do this, set unusedUnique to  
    #be edgeList, and then remove any edges that have been 
    #used from it to obtain the unused edges in the simulation
    unusedUnique = edgeList[:]    
    for edge in usedUnique:
        if edge in unusedUnique:
            unusedUnique.remove(edge)
    
    #If e = 0, then there is one flow path for the traffic so
    #I want this path to be returned by the function along 
    #with the other bits of useful information computed by the
    #function
    if e == 0:
    
        #Compute the Dijkstra path
        flowPattern = Dijkst(ist,iend,wei0)
        
        #Add 1 to each of the python-indexes in flowPattern so 
        #we obtain the nodes of the Dijkstra path and not 
        #their pyhton-indexes
        for i in range(len(flowPattern)):
            flowPattern[i] = flowPattern[i] + 1  
        
        #Return all of this data
        return carFreq,maximumLoad,mostCongested,loadCongested,unvisitedNodes,usedUnique,unusedUnique,flowPattern
    
    #If e != 0, just return the relevant information and no 
    #path since there wouldn't be a unique path taken by all
    #of the cars injected into the network as there would be 
    #so many different paths taken by different sets of cars
    else:
        return carFreq,maximumLoad,mostCongested,loadCongested,unvisitedNodes,usedUnique,unusedUnique





if __name__ == '__main__':
    import numpy as np
    import scipy as sp
    import math as ma
    import sys
    import time
    import csv    

    #Initialise weight matrix using the given data files 
    RomeX = np.empty(0,dtype=float)
    RomeY = np.empty(0,dtype=float)
    with open('RomeVertices','r') as file:
        AAA = csv.reader(file)
        for row in AAA:
            RomeX = np.concatenate((RomeX,[float(row[1])]))
            RomeY = np.concatenate((RomeY,[float(row[2])]))
        file.close()
        
        RomeA = np.empty(0,dtype=int)
        RomeB = np.empty(0,dtype=int)
        RomeV = np.empty(0,dtype=float)
        with open('RomeEdges','r') as file:
            AAA = csv.reader(file)
            for row in AAA:
                RomeA = np.concatenate((RomeA,[int(row[0])]))
                RomeB = np.concatenate((RomeB,[int(row[1])]))
                RomeV = np.concatenate((RomeV,[float(row[2])]))
            file.close()    
            
    #Make a copy of the original weight matrix using the 
    #provided function calcWei and the data from RomeVertices
    #and RomeEdges, and call this original weight matrix 
    #weiOriginal
    weiOriginal = calcWei(RomeX,RomeY,RomeA,RomeB,RomeV)      
    
    
    
    #Using RomeA and RomeB which came from the file RomeEdges,
    #I make an edge list of all of the edges which can be
    #referred to when the traffic flow is modelled below 
    #within the function.
    #Note: I will make a list of lists, so each element of the
    #list is in fact a list consisting of an edge in the
    #network
    edgeList1 = []
    for i in range(np.size(RomeA)):
        edge = []
        edge.append(RomeA[i])    
        edge.append(RomeB[i])
        edgeList1.append(edge)
    
    #Ensure that the edge list has no repeated edges in it for
    #when it is referred to when the function trafficFlow is
    #called upon below 
    edgeList = []
    [edgeList.append(x) for x in edgeList1 if x not in
                                                     edgeList]
    #Make a copy of the original weight matrix to use in the
    #first 2 models
    wei0 = weiOriginal.copy()
    
    #Model 1: start at node 13 (python-index 12) and end at
    #node 52 (python-index 51) with e = 0.01 and no blocked
    #nodes
    carFreq1,maximumLoad1,mostCongested1,loadCongested1,unvisitedNodes1,usedUnique1,unusedUnique1 = trafficFlow(12,51,wei0,200,20,180,0.01)
    
    #Model 2: same as Model 1, but for e = 0 instead
    carFreq2,maximumLoad2,mostCongested2,loadCongested2,unvisitedNodes2,usedUnique2,unusedUnique2,flowPattern1 = trafficFlow(12,51,wei0,200,20,180,0)
    
    #Get the weight matrix for Model 3 where node 30 is
    #blocked, so that any edges linked to node 30
    #(python-index 29) are blocked off so we do not use any of
    #the edges connected to node 30, and call it wei1
    wei1 = weiOriginal.copy()
    wei1[29,:] = 0
    wei1[:,29] = 0
    
    #Model 3: same as Model 1, but with node 30 blocked
    carFreq3,maximumLoad3,mostCongested3,loadCongested3,unvisitedNodes3,usedUnique3,unusedUnique3 = trafficFlow(12,51,wei1,200,20,180,0.01)