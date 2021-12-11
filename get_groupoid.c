//Created on October 4th 2017
//Author: Flaviano Morone

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#define MAX_DEGREE 1000
#define NAME_LENGTH 100

char line[MAX_DEGREE];

int ***A;
char **gene_name;

typedef enum{NO, YES} no_yes;
typedef enum{DOWN = 1, UP = 2} down_up;

typedef struct{
	int color;
	int in_degree;
}varGene;

varGene *gene;      

int init_gene(int N) {
	int i, col;
	col = 2;
	for(i = 1; i <= N; i++) {
		if(gene[i].in_degree == 0){
			gene[i].color = col;
			col++;
		}
		else 
			gene[i].color = 1;
	}
	return (col - 1);
}


int make_network(const char *network, int num_link_types) {
	int i, j, k, node, N;
	char *start;
    FILE *list;
	
    list = fopen(network, "r");
    N = 0;
    while( fgets(line, MAX_DEGREE, list) != NULL){
        N++;
    }
    fclose(list);

	A = (int ***)calloc(N + 1, sizeof(int **));
	for(i = 0; i <= N; i++)
		A[i] = (int **)calloc(N + 1, sizeof(int *)); 
	for(i = 0; i <= N; i++)
		for(j = 0; j <= N; j++)
			A[i][j] = (int *)calloc(num_link_types + 1, sizeof(int));

    gene = (varGene *)calloc(N + 1, sizeof(varGene));

	gene_name = (char **)calloc(N+1, sizeof(char *));
	for(i = 0;i <= N; i++)
		gene_name[i] = (char *)calloc(NAME_LENGTH, sizeof(char));

	list = fopen(network, "r");
    i = 1;
    while( fgets(line, MAX_DEGREE, list) != NULL){
        start = line;
        j = 1;
        while( sscanf(start, "%d%n", &node, &k) == 1) {
			if(node == -1){
				A[i][j][DOWN] = 1;
				gene[i].in_degree += 1;
			}
			if(node == 1){
				A[i][j][UP] = 1;
				gene[i].in_degree += 1;
			}
			start += k;
			j++;
		}
		i++;
	}
    fclose(list);

	for(i = 1;i <= N; i++){
		sprintf(gene_name[i], "node_%d has color", i);
	}

	return N;
}

void free_all(int N){
	int i, j;

	for(i = 0; i <= N; i++){
		for(j = 0; j <= N; j++){
			free(A[i][j]);
		}
	}
	for(i = 0; i <= N; i++){
		free(A[i]);
	}
	free(A);
	for(i = 0;i <= N; i++){
		free(gene_name[i]);
	}
	free(gene_name);
	free(gene);
}

no_yes have_same_colors(int **vj, int **vi, int num_colors, int num_link_types) {
	int n, t, diff;
	
	diff = 0; 
	for(n = 1; n <= num_colors; n++) {
		diff += abs(vj[n][DOWN] - vi[n][DOWN]);
		diff += abs(vj[n][UP] - vi[n][UP]);
	}
	return (diff > 0 ? NO : YES);
}

int get_new_partition(int N, int num_colors, int num_link_types, int initial_colors) {
	int i, j, new_num_colors, cnt;
	int ***table;
	
	table = (int ***)calloc(N + 1, sizeof(int **));
	for(i = 0; i <= N; i++)
		table[i] = (int **)calloc(num_colors + 1, sizeof(int *));

	for(i = 0; i <= N; i++)
		for(j = 0; j <= num_colors; j++)
			table[i][j] = (int *)calloc(num_link_types + 1, sizeof(int));
			
	for(i = 1; i <= N; i++) {
		for(j = 1; j <= N; j++) {
			
			if(A[i][j][DOWN] == 1)
				table[i][gene[j].color][DOWN]++;
			if(A[i][j][UP] == 1)
				table[i][gene[j].color][UP]++;
		}
	}	
	new_num_colors = initial_colors;
	for(i = 1; i <= N; i++) {
		if(gene[i].in_degree > 0) {
			cnt = 1;
			for(j = 1; j < i; j++) {
				if( have_same_colors(table[j], table[i], num_colors, num_link_types) ){
					gene[i].color = gene[j].color;
					break;
				}
				cnt++;
			}
			if(cnt == i) {
				new_num_colors++;
				gene[i].color = new_num_colors;
			}
		}
	}		
	for(i = 0; i <= N; i++)
		for(j = 0; j <= num_colors; j++)
			free(table[i][j]);
	for(i = 0; i <= N; i++)
		free(table[i]);
	free(table);
	return new_num_colors;
}


int main(int argc, char *argv[]){
	int i, j, N, cnt, num_colors, new_num_colors, initial_colors; 
	int num_link_types = 2; // Number of LINK equivalence classes (e.g. up/down, +1/-1, ...)
	 const char *network;

	network = argv[1];

    N = make_network(network, num_link_types);

	// Arabinose 
    //gene 1: CRP
	//gene 2: araC
	//gene 3: araBAD
	// A[2][1][UP] = 1;     This means that gene2 (araC) is UP-regulated by gene1 (CRP), so the link is 1-->2. 
	
	
	initial_colors = init_gene(N);	
	num_colors = initial_colors;
	new_num_colors = get_new_partition(N, num_colors, num_link_types, initial_colors);
	
	while(new_num_colors != num_colors) {
		num_colors = new_num_colors;
		new_num_colors = get_new_partition(N, num_colors, num_link_types, initial_colors);
	}
	cnt = 0;
	for(i = 1; i <= N; i++){
		if(gene[i].color == 1) {
			cnt = 1;
			break;
		}
	}
	if(cnt == 0){
		for(i = 1; i <= N; i++)
			gene[i].color -= 1;
	}
	printf("\n\n\t\t Groupoids of %s\n\n", network);
	for(i = 1; i <= N; i++)
		printf("%s ---> color # %d\n", gene_name[i], gene[i].color);
	printf("\n\n");	
	

    free_all(N);	
	return 0;
}















