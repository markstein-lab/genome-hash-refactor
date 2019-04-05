//Copyright (c) World Internet Productions, LLC, 1999
//Copyright (c) in silico Labs, LLC, 2006
#ifndef group_defs
#define group_defs

extern char lookup[];
extern char lookrev[];
#define fetchDNA(X) lookup[(DNAstring[(X)/4]>> (6 - 2*((X)&3)))&3]
#define fetchRevDNA(X) lookrev[(DNAstring[(X)/4]>> (6 - 2*((X)&3)))&3]

typedef struct chain chain;
struct chain{
	char* position;
	chain *next;
	chain *previous;
};
typedef struct {
	char *key;
	int length;
	unsigned int hash;
	chain *first;
	chain *last;
	int overflow;
	char flags;
	char palindrome;
	char logical_name;
} hashEntry, *hashPtr;
typedef struct ht hashTable, *htPtr;
struct ht {
	unsigned int mask;
	int tablesize;
	hashEntry *htable;
	htPtr next;
	int nextavail;
};
extern int errorflag;
extern int extablesize;
extern unsigned int* extable;
extern int exindex;
extern htPtr mainHashTable;
#define TRUE    1
#define FALSE   0
unsigned int hashSubsequence(char *vector, int len);
htPtr newTable(int size);
hashPtr hashRegister(htPtr table, char *vector, int patternlen);
hashPtr hashQuery(htPtr table, unsigned int where, int patternlen);
hashPtr hashFindAndListFast(htPtr table, unsigned int where, int patternlen, int *pos, char **pat, unsigned int rawhash);
#endif				//group_defs

