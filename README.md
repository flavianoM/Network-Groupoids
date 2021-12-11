Created in 2017

Author: Flaviano Morone

# Network-Groupoids
Groupoid Algorithm to find balanced network partitions 

To compile the source code use the command

 gcc -o get_groupoid get_groupoid.c -lm -O3

The program requires as input the file containing the network.

To run the executable use the command:

./get_groupoid NETWORK_FILENAME

where - NETWORK_FILENAME - is the file containing the network.

The data structure in the file containing the network must be formatted as an adjacency matrix A_ij with 3 admissibile values: A_ij = -1,0,1.  

At output, the program prints at screen the node ID and its "color ID". Nodes with same color belong to the same balanced partition. 



For more information, visit https://www.pnas.org/content/pnas/117/15/8306.full.pdf
