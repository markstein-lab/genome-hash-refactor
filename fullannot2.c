//Copyright (c) World Internet Productions, LLC, 1999
//Copyright (c) in silico Labs, LLC, 2001
//modified by in silico Labs,LLC 10/4/01 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <malloc.h>
#include "group_defs2.h"




void plusmsg(char message[], int i, int pos, int posend);
void minusmsg(char message[], int i, int pos, int posend);

void makerevcompl(char *from, char *to, int num);

typedef   struct ann_index {
	int annot_number;
	int start;
	int stop;
	char direction;
	char *gene_name;
   } ann_index;
ann_index *map;
int map_size;
int gene_names;
int margin;

void fullannot (int geneindex, int pos, int posend, char *gene_filter, 
				char *results) {
	
	//geneindex is a condensation of the annotations - one entry per gene
	//pos is the first position in the cluster
	//posend is the last position in the cluster
	//gene_filter is no longer used
	//results is the character string that describes where the cluster was found
	
	int i, j;  
	int lfilter;
	char message[3000], mess1[3000], extramessage[10000];
	
	
	// get length of gene filter
	lfilter = strlen(gene_filter);
	//initialize output strings to empty
	//message[0] = mess1[0] = extramessage[0] = 0;
	strcpy(message, ""); strcpy(mess1, ""); strcpy(extramessage, "");
    
	for (i = geneindex; 
	i < map_size && map[i].annot_number == map[geneindex].annot_number;
	i++) {
		if (map[i].start > posend) break; 
		if (map[i].start <= pos && map[i].stop >= posend) {
			//it is in gene at map[i]
			i++;
			break;
		}
	}
	if (i == map_size ||
		map[i].annot_number != map[geneindex].annot_number) {
		//Cluster is after all genes
		if (map[i-1].direction == '+') {
			plusmsg(message, i-1, pos, posend);
		}
		else {
			minusmsg(message, i-1, pos, posend);
		}
		sprintf(results, "\t%s", message);
		return;
	}
	if (i == geneindex || posend > map[i-1].stop) {
		// Normal case: cluster is between two genes
		if (i > geneindex) {
			if (map[i-1].direction == '+') plusmsg(mess1, i-1, pos, posend);
			else minusmsg (mess1, i-1, pos, posend);
		} else mess1[0] = 0;  //null string
		if (map[i].direction == '+') plusmsg(message, i, pos, posend);
		else minusmsg(message, i, pos, posend);
		if (strlen(mess1)) {
			if (gene_names) {
				if (strlen(message)) {
					sprintf(results, "%s\n%s", mess1, message);
				}
				else
					sprintf(results, "%s", mess1);
			}
			else {
				sprintf(results, "\t%s,\n\t%s", mess1, message);
			}
		}
		else {
			if (gene_names) {
				sprintf(results, "%s", message);
			}
			else {
			    sprintf(results, "\t%s", message);
			}
		}
	}
	else {
		// previous gene overlaps cluster
		if (map[i-1].direction == '+') plusmsg(message, i-1, pos, posend);
		else minusmsg(message, i-1, pos, posend);
		if (gene_names) {
			sprintf(results, "%s", message);
		}
		else {
			sprintf(results, "\t%s", message);
		}
	}
	for (j = i+1; map[j].annot_number == map[i].annot_number; j++) {
		if(map[j].start > map[i-1].stop) break;
		if (map[j].start > posend && map[j-1].start <= pos ||
			map[j].start <= pos && map[j].stop >= posend ||
			map[j].start <= posend && map[j].stop > posend)
		{
			fullannot(j-1, pos, posend, gene_filter, extramessage);
			strcat(results, "\n");
			strcat(results, extramessage);
		}
	}
}

void plusmsg(char message[], int i, int pos, int posend) {
	if (posend < map[i].start) {
		if (gene_names) {
			if (margin == 0 || map[i].start - pos <= margin) {
				sprintf(message, "%s", map[i].gene_name);
			}
		}
		else {
			sprintf(message, "%d bp upstr of %s(+)",
				map[i].start - posend, map[i].gene_name);
		}
	}
	else if (pos < map[i].start) {
		if (gene_names) {
			sprintf(message, "%s", map[i].gene_name);
		}
		else {
			sprintf(message, "%d bp upstr of %s(+) and overlaps",
				map[i].start - pos, map[i].gene_name);
		}
	}
	else if (posend <= map[i].stop) {
		if (gene_names) {
			sprintf(message, "%s", map[i].gene_name);
		}
		else {
			sprintf(message, "%d bp inside and totally within %s(+)",
				pos - map[i].start, map[i].gene_name);
		}
	} 
	else if (pos <= map[i].stop) {
		if (gene_names) {
			sprintf(message, "%s", map[i].gene_name);
		}
		else {
			sprintf(message, "%d bp from end of %s(+) and overlaps", 
				map[i].stop - pos, map[i].gene_name); 
		}
	}
	else {
		if (gene_names) {
			if(margin == 0 || posend - map[i].stop <= margin) {
				sprintf(message, "%s", map[i].gene_name);
			}
		}
		else {
			sprintf(message, "%d bp dnstr of %s(+)",
				pos - map[i].stop, map[i].gene_name);
		}
	}
}

void minusmsg(char message[], int i, int pos, int posend) {
	if (posend < map[i].start) {										
		if (gene_names) {
			if (margin == 0 || map[i].start - pos <= margin) {
				sprintf(message, "%s", map[i].gene_name);
			}
		}
		else {
			sprintf(message, "%d bp dnstr of %s(-)",
				map[i].start - posend, map[i].gene_name);
		}
    }
	else if (pos < map[i].start) {
		if (gene_names) {
			sprintf(message, "%s", map[i].gene_name);
		}
		else {
			sprintf(message, "%d bp from end of %s(-)and overlaps",
				posend - map[i].start, map[i].gene_name);
		}
	}
	else if (posend <= map[i].stop) {
		if (gene_names) {
			sprintf(message, "%s", map[i].gene_name);
		}
		else {
			sprintf(message, "%d bp inside and totally within %s(-)",
				map[i].stop - posend, map[i].gene_name);
		}
	} 
	else if (pos <= map[i].stop) {
		if (gene_names) {
			sprintf(message, "%s", map[i].gene_name);
		}
		else {
			sprintf(message, "%d bp upstr of %s(-) and overlaps", 
				posend - map[i].stop, map[i].gene_name); 
		}
	}
	else {
		if (gene_names) {
			if (margin == 0 || posend - map[i].stop <= margin) {
				sprintf(message, "%s", map[i].gene_name);
			}
		}
		else {
			sprintf(message, "%d bp upstr of %s(-)",
				pos - map[i].stop, map[i].gene_name);
		}
	}
}
