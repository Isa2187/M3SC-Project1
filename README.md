# M3SC-Project1
First project of Scientific Computation (M3SC) module taken in 3rd year. (Grade = 100%)

All  code was done in Python, and details of the task can be found in the folder Project_Files.

The task involved modelling Rome as a network of nodes, and simulating the flow of cars through the city. At each minute, new cars would be injected into the network at a starting node, and start making their way through the network along the shortest route, which is computed using Dijkstra's algorithm. To account for many cars taking the same route, the weight matrix representing the edges between adjacent nodes was updated to account for this, meaning some cars would take different routes to the finishing node than others.

To account for real-life scenarios, we investigated how changing certain parameters would affect the traffic flow. These included blocking off one node to simulate an accident, which would block off all incoming and outgoing edges from this road.
