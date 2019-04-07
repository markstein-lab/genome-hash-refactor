//    copyright (c) in silico Labs, LLC 2005
#include <stdlib.h>
#include <stdio.h>
//#include <malloc.h>
#include "group_defs2.h"

char *DNAstring;


#define FOLD   hash += ((hash>>11) + (hash>>20));
char lookup[4] = {'A', 'C', 'G', 'T'};
char lookrev[4] = {'T', 'G', 'C', 'A'};



htPtr newTable(int size) {
	htPtr htptr;
	if (size & (size - 1)) {
		printf("Table size for newTable must be a power of 2; size given was %d\n",
			    size);
	}
	htptr = calloc(1, sizeof(hashTable));
	htptr->tablesize = size;
	htptr->mask = size - 1;
	htptr->htable = calloc(size, sizeof(hashEntry));
	htptr->nextavail = 1;		//cannot start overflow at 0
	return htptr;
}


unsigned int hashSubsequence(char *vector, int len) {
	int result = 0;
	int i;
	for (i = 0; i < len; i++) {
		switch (vector[i]) {
		case 'A':
  		case 'a':
			result <<= 2;
			break;
		case 'C':
  		case 'c':
			result = 1 + (result<< 2);
			break;
		case 'G':
  		case 'g':
			result = 2 + (result<< 2);
			break;
		case 'T':
  		case 't':
			result = 3 + (result<< 2);
			break;
		default:
			result += 0;
//			printf("Error in position %d, character is %c\n",
//				i, vector[i]);
			errorflag = TRUE;
		}
	}
	return result;
}

//Create a registry of patterns 
hashPtr hashRegister(htPtr table, char *vector, int patternlen) {
	unsigned int hash, h;
	int j;
	char  *copy, *t;
	errorflag = FALSE;
	hash = hashSubsequence(vector, patternlen);
	FOLD;
	if (errorflag == TRUE) return NULL;
	h = hash & table->mask;
	while(table->htable[h].key != 0) {
		if (hash != table->htable[h].hash) goto next;
		if (patternlen != table->htable[h].length) goto next;
		for (j = 0; j < patternlen; j++) {  //does sequence match hash entry?
			if (vector[j] != table->htable[h].key[j]) goto next;
		}		//we found a match

		return NULL;    //don't re-enter item already in registry
next:
		if (table->htable[h].overflow != 0) {//if there is an overflow chain,
			h = table->htable[h].overflow;	//follow it!
			continue;
		}
		for (j = table->nextavail; j < table->tablesize; j++) {
			//search for an overflow slot;
			if (table->htable[j].key == 0) break;	//empty slot found
		}
		if (j == table->tablesize) {
			//we cannot find more space in this hashtable
			//search the next hashtable, or create a new one, and search it
			if (table->next == NULL)
				table->next = newTable(2*table->tablesize);	//next table twice as big
			return hashRegister(table->next, vector, patternlen);
		}
		table->htable[h].overflow = j;		//indicate overflow position
		table->nextavail = j+1;
		h = j;
		break;
	}
	// a new position has been found. Make it point to current substring, 
	//and set up chain of pointers to instances of this sequence
	copy = calloc(1,patternlen+1);
	t = copy;
	for (j = 0; j < patternlen; j++) {
	   *t++ = *vector++;   //make a copy of the pattern
        }
	*t = 0;                //zero-terminate copied string
	table->htable[h].key = copy;
	table->htable[h].length = patternlen;
	table->htable[h].hash = hash;
	return &table->htable[h];
}

hashPtr hashFindAndListFast(htPtr table, unsigned int where, int patternlen, 
			    int* pos, char* *pattern, unsigned int rawhash) {
	unsigned int hash, h;
	unsigned int j, jptr, jj;

	errorflag = FALSE;
	hash = rawhash;
	FOLD;
	if (errorflag == TRUE) return NULL;
	h = hash & table->mask;
	while(table->htable[h].key != 0) {
		if (hash != table->htable[h].hash) goto next;
		if (patternlen != table->htable[h].length) goto next;
		for (j = 0; j < (unsigned)patternlen; j++) {  //does sequence match hash entry?
			jptr = where+j;
			if (fetchDNA(jptr) != table->htable[h].key[j]) goto next;
		}		//we found a match
		// wait a minute -- check for Ns in the genome
		if ((where+j)/32 < extable[2*exindex]) goto goodcompare;  //all is ok
		while(where/32 < extable[2*exindex]) exindex++;
		if (where/32 == extable[2*exindex])	{
			for (jj = 0; jj < (unsigned)patternlen; jj++)	{
				if (extable[2*exindex+1] & (1<<(31-((where+jj)&31)))) goto next;
				if (((where + jj + 1)&31) == 0) { //finished a word of exceptions
					//maybe no more exception words overlap the pattern
					if (((where + j - 1)/32) < extable[2*exindex+2]) goto goodcompare;
					//skip pattern characters which do not overlap next
					//exception word
					while ((where + jj + 1)/32 < extable[2*exindex +2]) {
						jj += 32;
					}
					exindex++;  //continue looking with next word of exceptions
				}
			}
		}
goodcompare:		
		*pos = where;
		*pattern = table->htable[h].key;
		return &table->htable[h];	//return pointer to hashEntry
next:
		if (table->htable[h].overflow != 0) {
			//if there is an overflow chain,
			h = table->htable[h].overflow;	//follow it!
			continue;
		}
		
		//we cannot find pattern in this hashtable
		//search the next hashtable, return NULL
		
		if (table->next == NULL) {
			return NULL;
		}
		
		return hashFindAndListFast(table->next, where, patternlen, 
			pos, pattern, rawhash);
		
	}
	return NULL;
}

hashPtr hashQuery(htPtr table, unsigned int where, int patternlen) {
	unsigned int hash, h;
	char vector[500];
	unsigned int j, jptr;

	errorflag = FALSE;
	for (j = 0; j < (unsigned)patternlen; j++) {	//translate to ACGT char. codes
		jptr = where + j;
		vector[j] = fetchDNA(jptr);
	}
	vector[j] = 0;	 //end string
	hash = hashSubsequence(vector, patternlen);
	FOLD;
	if (errorflag == TRUE) return NULL;
	h = hash & table->mask;
	while(table->htable[h].key != 0) {
		if (hash != table->htable[h].hash) goto next;
		if (patternlen != table->htable[h].length) goto next;
		for (j = 0; j < (unsigned)patternlen; j++) {  //does sequence match hash entry?
			if (vector[j] != table->htable[h].key[j]) goto next;
		}		//we found a match 
		return &table->htable[h];	//return pointer to hashEntry
next:
		if (table->htable[h].overflow != 0) {
			//if there is an overflow chain,
			h = table->htable[h].overflow;	//follow it!
			continue;
		}
		
		
		//It is not in this hash table
		//search the hashtable extension (table->next) if it exists
		//or return NULL
		if (table->next == NULL) {
			return NULL;
		}
		return hashQuery(table->next, where, patternlen);
		
	}
	return NULL;
}

