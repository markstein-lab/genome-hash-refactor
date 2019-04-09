//copyright (c) in silico Labs, LLC 2004

#include "group_defs2.h"
typedef struct ann_index {
		int annot_number;
		int start;
		int stop;
		char direction;
		char *gene_name;
} ann_index;
ann_index *map;
char** annot_name;
char** chromo_name;
unsigned int* arm;
int annotation_search(unsigned int chr, int index);
int findchrm(unsigned int);

int findchrm(unsigned int x)
{
		int i;
		for (i = 1; arm[i-1] < arm[i]; i++) {
				if ( x < arm[i]) break;
		}

		return i-1;
}

int locategene(unsigned int where, int neargene)
{
		int j, index, i;

		//find contig in which argument lies
		j = findchrm(where);

		index = annotation_search(j, where - arm[j]);
		//index now points to first gene on appropriate chromosome
		for (i = index; map[i].annot_number == j; i++) {
				if (where >= arm[j] + map[i].start - neargene &&
					where <= arm[j] + map[i].stop + neargene) {
						return (arm[j] + map[i].stop + neargene - where + 1);
				}
		}
		return 0;
}
