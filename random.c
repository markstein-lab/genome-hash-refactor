//Copyright (c) World Internet Productions, LLC, 1999, 2000, 2001
//Copyright (c) in silico Labs, LLC. 2004, 2005, 2006

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "group_defs2.h"
#include "timer.h"
#include "newbool.h"

#define BORDERWIDTH 10
void readgenome();
void readinput(char line[]);
int ncontigs;
int extablesize;
int exindex;
unsigned int* extable;
unsigned int DNALen;
int errorflag;
htPtr mainHashTable;
int flank[2 * BORDERWIDTH][4];
char *DNAstring;
char *motif[50];
void makerevcompl(char *from, char *to, int num);
void recursiveRegister(char *line, char* target[], int *patnum); 
void insertClusters(unsigned int **av, char ***pv, char line[]);
void hsort2(unsigned int x[], unsigned int y[], int n);
int do_the_search(int start, int stop, unsigned int *av, char **pv, int slide_size,
				   unsigned int slide[], unsigned int patsize[]);
FILE *save;
FILE *insource;
FILE *index_file;
FILE *archive;
int map_size;
char pattern_id;
typedef   struct ann_index {
	int annot_number;
	int start;
	int stop;
	char direction;
	char *gene_name;
   } ann_index;
ann_index *map;
int annotation_search(int chr, int index);
void getLine(char line[],FILE *f);
void fullannot(int, int, int, char *, char *);
int readContigs(FILE *qnames);
char** annot_name;
char** chromo_name;
unsigned int* arm;
unsigned int search_start, search_stop;
int remote_input;
char* booleansyntax(char message[]);
int locategene(unsigned int where, int neargene);
int gene_names;
int margin;
int findchrm(unsigned int x) ;
void readgenenames(FILE *allnames, int syn_count);
int* gene_number;
char  **gene_name;
void printallsynonyms(char* resultstring, int syn_count);
int xdata[500];
int ydata[500];

int main (int argc, char *argv[])
{
   char result[100];
   FILE *stream;
   FILE *qnames;  
   FILE *allnames;
   char line[500], comparestring[200], hitlist[500], decnumber[20];
   char resultstring[100000], signature[250];
   int i, j, k, ix, ichr;
   int winLen, linelength;
   int slide_size, group_count, maxpatsize;
   unsigned int ref;
   int patnum, minGrpSize, groupsize, oldpatnum;
   char* target[50000];
   unsigned int endcluster;
   char temp[500];
   char boolean_expression[100], gene_filter[100];
   unsigned int *addrvec, *av, *ax, *ay, *az, *aq;
   char **patvec, **pv, **px, **py, **pz, **pq;
   unsigned int cp;
   hashPtr hpt;
   //void tally(unsigned int where, int patLen);
   unsigned int rawhash, slide[50], patsize[50];
   int annotation_index;
   int full_results;
   int counter;
   int num_contigs;
   int repeat;
   int numPatterns;
   char count_only;
   TIMEVAL time1, time2;
   int counts[27][3];
   char longsig[500];
   char region_name[50];
   unsigned int temploc, temploc2;
   char *clusterDNA;
   char* remoteinput[14];
   int remotectr;
   int syn_count;
 
   int dist;
   int varsfound;

   double instances, p, q, sigma;
   double p0, p1, p2, pmore;

   count_only = ' ';
   result[0] = 0;
   pattern_id = 'A';
   full_results = 0;
   remotectr = 0;
   insource = fopen("enhancerinput.txt", "r");
   if (insource) {
	   remote_input = 1;
	   for (i = 0; i < 14; i++) {
		   readinput(line);
		   remoteinput[i] = malloc(strlen(line) + 2);
		   strcpy(remoteinput[i], line);
	   }														  
   }
   else	{
	   remote_input = 0;
       insource = stdin;
   }
   for (i = 1; i < argc; i++) {
      strcpy(result, argv[i]);
      break;
   }
   if (result[0] == 0) strcpy(result, "tempres");
   save = fopen(result, "w");
   if (save == NULL) {
      printf("Could not open output file %d\n", result);
      exit(99);
   }
   archive = 0;
   stream = NULL;
   errorflag = 0;
   slide_size = 0;
   counter = 0;                      //point to start of DNAstring
   ix = 0;
   printf("Reading in the genome\n");
   qnames = fopen("datasize.txt", "r");
   if (qnames == NULL) {
	   printf("No datasize.txt file in this directory -- quitting\n");
	   exit(0);
   }
   fscanf(qnames, "%d %d %d %d %d", &DNALen, &map_size, &ncontigs, &extablesize,
	   &syn_count);
printf("syn_count = %d\n", syn_count);
//  printf("genome size = %u; number of genes = %d, number of contigs = %d\n",
//	   DNALen, map_size, ncontigs);
   fclose(qnames);
   DNAstring = malloc(DNALen/4 + 8);
   if (DNAstring == 0) {
	   //No space available!
	   printf("Not enough space for genome.\n");
	   exit(0);
   }
   arm = malloc((ncontigs+2) * sizeof(int*));
   annot_name = malloc ((ncontigs+2)* sizeof(char*));
   chromo_name = malloc ((ncontigs+2)* sizeof(char*));
   extable = malloc ((extablesize + 1) * 2 * sizeof (int));
   map = malloc (map_size * sizeof(ann_index));
   if (map == 0) {
	   //No space available!
	   printf("Not enough space for the gene map!\n");
	   exit(0);
   }
   
  
   qnames = fopen("xcontigs.txt", "r");
   k = 0;
   k = readContigs(qnames); 
   if (k != ncontigs) printf("wrong number of contings, wanted %d, found %d\n", ncontigs, ix);
   fclose(qnames);
   num_contigs = k;						     //ix is the number of chromosomes
   arm[k] = DNALen;                //end of genome
   arm[k+1] = 0;
   chromo_name[k] = malloc(20);
   strcpy(chromo_name[k], "end_of_genome");
   chromo_name[k+1] = chromo_name[k];
   annot_name[k] = chromo_name[k];
   annot_name[k+1] = annot_name[k];
   qnames = fopen("exceptions.txt", "r");
   for (i = 0; i < extablesize; i++) {
	   getLine(line, qnames);
	   sscanf(line, "%d %8.8ullx", &extable[2*i], &extable[2*i+1]);
	   if (feof(qnames)) {
		   printf("Exception table too short\n");
		   exit(0);
	   }
   }
   extable[2*i] = (DNALen/32) + 1;    //dummy last entry beyond the genome
   extable[2*i + 1] = 0xFFFFFFFF;
   fclose (qnames);
   getTime(&time1);
   readgenome(); //read in the genome itself
   getTime(&time2);
   printf("Time to read genome was %g seconds\n", getDiffMillisecs(&time1, &time2)/1000.0);
/*  to check that N's separate the contigs
   for (i = 0; i < ix; i++) {
	   //relocate items in arm to point into the DNA vector
	   //check that inter-contig points have an 'N'
	   if (*arm[i] != 'N') {
		   printf("Contig %d starts with %d (%c), should have been an 'N'\n",
			   i, *arm[i], *arm[i]);
	   }
   }
*/
   printf("DNALen = %u\n", DNALen);
 
   index_file = fopen("xchrdata.txt", "r");
   if(index_file == NULL) {
      printf("Annotation index file not found\n");
   }
   else {
	   for (i = 0; i < map_size; i++) {
		   if (feof(index_file)) {
			   printf("premature end of gene info file!\n");
			   exit (0);
		   }
		   for (j = 0; j < 10; j++) line[j] = 0;
		   getLine(line, index_file);
		   j = sscanf(line, "%d%d%d %c %498[^\n\r]",&map[i].annot_number, &map[i].start,
			   &map[i].stop, &map[i].direction, temp);
		   map[i].gene_name = malloc(strlen(temp) + 1);
		   strcpy(map[i].gene_name, temp);
	   }
   }
   if (i != map_size) printf("chrdata file did not have right number of records.\n");
   fclose(index_file);


   allnames = fopen("xallnames.txt", "r");

   if (allnames == NULL) {
	   syn_count = 0;
	   printf("File of gene-names is missing\n");
	   exit(1);
   }
   gene_number = malloc(syn_count * sizeof(int));
   gene_name = malloc(syn_count * sizeof(char**));


   if(syn_count) readgenenames(allnames, syn_count);


   mainHashTable = newTable(4096);

   patvec = malloc(10000000*4);
   addrvec = malloc(10000000*4);
   
   patnum = 0;
   oldpatnum = 0;
   av = addrvec;
   pv = patvec;

search_limiters:
   printf("Enter the chromosome name where to localize search, or hit enter\n");
   readinput(region_name);
   if (strlen(region_name) == 0) {
      //scan the entire genome
      search_start = 0;
      search_stop = (DNALen/4);
   }
   else {
      for (i = 0; i < ncontigs; i++) {
         if (strcmp(chromo_name[i], region_name)) continue;
         printf("Enter the starting position on chromosome:");
		 readinput(line);
         sscanf(line,"%d", &temploc);
         search_start = (arm[i]+temploc)/4;
		 j = i+1;
		 while(strcmp(chromo_name[i], chromo_name[j]) == 0) j++;
		 if (arm[i] + temploc > arm[j]) {
			 printf("starting location is beyond chromosome end\n");
			 printf("Chromosome %s has %d base pairs\n", region_name, arm[j]-arm[i]);
			 printf("re-enter search limit parameters\n");
			 goto search_limiters;
		 }
         printf("Enter the stopping position on chromosome:");
		 readinput(line);
         sscanf(line, "%d", &temploc2);
		 printf("Temploc2 = %d\n", temploc2);
		 if (temploc2 == 0) temploc2 = arm[j] - arm[i] - 1;
		 printf("Temploc2 = %d\n", temploc2);
         search_stop = (arm[i]+temploc2)/4;
         if (arm[i] + temploc2 > arm[j] ||
			 temploc2 < temploc) {
            printf("stopping location is beyond chromosome end\n");
			printf("or stopping location preceeds start point\n");
			printf("Chromosome %s has %d base pairs\n", region_name, arm[j]-arm[i]);
            printf("re-enter search limit parameters\n");
            goto search_limiters;
        }
        goto read_pattern;
      }
      printf("Chromosome name not found. Respecify search limits\n");
      goto search_limiters;
   }

read_pattern: 
for (repeat = 0; repeat < 20; repeat++) {
  if(repeat) {
    numPatterns = pattern_id - 'A';
    pattern_id = 'A';
    free (mainHashTable->htable);
    free (mainHashTable);
    mainHashTable = newTable(2048 * 8);

    patnum = 0;
    oldpatnum = 0;
    av = addrvec;
    pv = patvec;
    for (i = 0; i < numPatterns; i++) {
       for (j = 0; j < strlen(motif[i]); j++) {
          xdata[j] = Random();
          ydata[j] = motif[i][j];
       }
       hsort2(xdata, ydata, strlen(motif[i]));
       for (j = 0; j < strlen(motif[i]); j++) motif[i][j] = ydata[j];
       printf ("new pattern %c = %s\n", pattern_id, motif[i]);
       fprintf (save, "new pattern %c = %s\n", pattern_id, motif[i]);
       recursiveRegister(motif[i], target, &patnum);
       printf("Pattern generates %d instances\n", patnum - oldpatnum);
       fprintf(save, "Pattern generates %d instances\n", patnum - oldpatnum);
       oldpatnum = patnum;
       pattern_id++;
    }
 } else {
read_next_pattern:
   if (remote_input) {
	   if (remoteinput[remotectr][2] == 0) {
		   linelength = 0;
		   line[0] = 0;
	   }
	   else {
		   strcpy (line, remoteinput[remotectr]+3);
		   linelength = strlen(line);
		   printf("Linelength = %d, line is %s\n", linelength, remoteinput[remotectr]+3);
		   fprintf(save, "General pattern %c:= %s\n", pattern_id, line);
		   remotectr++;
	   }
   }
   else {
	   printf ("Type in the pattern to find (limit 500 bp)\n%c:= ", 
		   pattern_id);
	   readinput(line);
	   linelength = strlen(line);
	   fprintf(save, "General pattern %c:= %s\n", pattern_id, line);
   }
   if (linelength) {
	   if (line[0] == '"') {  //specified a file name
		   insertClusters(&av, &pv, line);
		   pattern_id++;
 		   goto read_pattern;
	   }
  
      for (i = 0; i < linelength; i++) line[i] &= 0xDF; //force upper case
      motif[pattern_id-'A'] = malloc(linelength+1);
      for (i = 0; i <= linelength; i++) motif[pattern_id-'A'][i] = line[i];
      recursiveRegister(line, target, &patnum);
      printf("Pattern generated %d sequences\n", patnum - oldpatnum);
	  if (remote_input && patnum == oldpatnum) {
		  fprintf(save, "ERROR - pattern %s contains inappropriate characters, no search patters generated\n", line);
		  exit(0);
	  }
      oldpatnum = patnum;
      pattern_id++;
      for (ix = 0; ix < slide_size; ix++) {
		  if (linelength == (int)patsize[ix]) goto read_next_pattern;
      }
      patsize[ix] = linelength;
      slide_size++;

      goto read_next_pattern;
   }
  }
  // patLen = strlen(target[0]);    //in this version, all patterns have same len.
   

   printf("Patience: constructing list of occurrances of target in DNAstring\n");
   getTime(&time1);
 
   //Construct the hash table

   maxpatsize = 0;
   rawhash = 0;
   for (i = 0; i < slide_size; i++) {
	   if (patsize[i] > 15) slide[i] = 0xFFFFFFFF;
	   else slide[i] = (1<<(2*patsize[i])) - 1;
   }
   exindex = 0;
   ix = do_the_search(search_start, search_stop, av, pv, slide_size, slide, patsize);
   av += ix;
   pv += ix;

   getTime(&time2);
   printf("Time to scan genome was %g seconds\n", getDiffMillisecs(&time1, &time2)/1000.0);
   hsort2(addrvec, (unsigned int *)patvec, av-addrvec);
   getTime(&time1);
   printf("Time to sort the hits was %g seconds\n", getDiffMillisecs(&time2, &time1)/1000.0);

   printf("Instances found = %d\n", av-addrvec);
   fprintf(save, "Instances found = %d\n", av-addrvec);

   instances = av - addrvec;
   p = instances/DNALen;
   q = 1.0 - p;
   printf("#motifs found = %g, prob of finding a motif = %g\n", instances, p);

   ref = 0;
   comparestring[0] = 0;    //not using a start pos input right now!
   boolean_expression[0] = 0;

next_grouping:
   count_only = 'C';
   if (repeat) goto skipParameters; 
   printf("\n==========================================================\n");
   printf("==========================================================\n");
   fprintf(save,
	  "==========================================================\n");
   if (remote_input) {
	   sscanf(remoteinput[12]+13, "%d", &winLen);
	   sscanf(remoteinput[11]+13, "%d", &minGrpSize);
	   strcpy(boolean_expression, remoteinput[13]+19);
	   printf("Read boolean expression #%s#\n", boolean_expression);
       if (booleansyntax(boolean_expression) == 0) {
		   printf("Ill formed boolean expression, run terminated \n");
		   exit(0);
	   }
	   goto start_clustering;
   }
   printf("Do you want to keep expanded results? (Y, N. enter==N, C=count only)\n");
   readinput(comparestring);
   full_results =
      ((comparestring[0] & 0xDF) == 'Y')? 1 : 0;
   count_only = comparestring[0] & 0xDF;
   comparestring[0] = 0;

read_boolean:
   printf("Enter a boolean expression or just hit enter:\n");  
   readinput(boolean_expression);
   if (booleansyntax(boolean_expression) == 0) {
	   printf("Ill formed boolean expression, re-enter \n");
	   goto read_boolean;
   }
   fprintf(save, "Boolean expression is %s\n", boolean_expression);
   
 
/**
   // check for filter	ll 

   printf("Enter gene filter characters (case sensitive) or just hit enter:\n");
   comparestring[0] = 0;
   readinput(comparestring);
   strcpy(gene_filter, comparestring);
   fprintf(save, "Controlling gene filter %s\n",gene_filter);
**/

read_winLen: 
   printf("What is the window size? (enter 0 to exit program): ");
   printf("Also enter a second integer n for limiting searches to only genes and\n");
   printf("   n bp upstream and downstrem of a gene\n");
   printf("   Add a third parameter (contents immaterial) to print out all\n");
   printf("   synonyms of genes near hits\n");
   comparestring[0]=0;
   margin = 0;
   readinput(comparestring);
   varsfound = sscanf (comparestring, "%d %d %s", &winLen, &margin, line);
   printf("%s\n", comparestring);
   gene_names = 0;
   if (varsfound <= 0 || winLen == 0) exit(0);
   if (varsfound >= 2) {								
	   fprintf(save, "Search within %d bp of any gene\n", margin);
	   printf("Search within %d bp of any gene\n", margin);
	   if (varsfound > 2) gene_names = 1;
   }
   else	{
	   if ((unsigned int) winLen > DNALen) {
		   printf("Window size inconsistant with DNA string size\n");
		   printf("Please re-enter\n");
		   goto read_winLen;
	   }	   
   }
   fprintf(save, "Window size is %d\n", winLen);
   
   printf("Enter the minimum group size to report: ");
   scanf("%d", &minGrpSize);	
   getchar();           //remove character that ended previous line
   printf("%d\n", minGrpSize);
   fprintf(save, "Minimum group size is %d\n", minGrpSize);
 

start_clustering:
skipParameters:
   if (archive) remove("enhancerres.txt");
   archive = fopen("enhancerres.txt", "w");

   p0 = pow(q, winLen - 1.0);
   p1 = (winLen - 1.0) * p * p0/q;
   p2 = (winLen - 2.0) * p / q / 2.0 * p1;
   pmore = 1.0 - p0 - p1 ;
   sigma = sqrt(instances * pmore * (1.0 - pmore));
   printf("prob clusters of 3 or more is %7.5f\n",
      pmore);
   printf("Expected number of clusters = %8.5f, std = %10.5f\n",
      instances * pmore, sigma);
   fprintf(save, "prob clusters of 3 or more is %7.5f\n",
      pmore);
   fprintf(save, "Expected number of clusters = %8.5f, std = %10.5f\n",
      instances * pmore, sigma);

   //for (i = 0; i < 2 * BORDERWIDTH; i++) for(j = 0; j < 4; j++) 
   //   flank[i][j] = 0;
   group_count = 0;
   for (ax = addrvec, px = patvec; ax < av; ax++, px++) {
	   aq = ax;
	   pq = px;
	   if (*aq - *addrvec < ref) continue;


	   groupsize = 1;
	   if (ax + 1 < av || minGrpSize == 1) {
		   if (*(ax+1) < *ax + winLen || minGrpSize == 1) {
			   //We may have a window!
			   for (ay = ax+1, py = px+1; ay < av; ay++, py++){
				   if (*ay < *aq + strlen(*pq)) continue;
				   if (*ay >= *ax + winLen) {
					   break;    //stop when ay is outside the window
				   }
				   groupsize++;	
				   aq = ay;
				   pq = py;
			   }  //close for ay loop
			   
			   //Now the total window runs is [ax, ay).
			   if (groupsize < minGrpSize) continue;
			   if (margin) { //is hit near enough to a gene?
				   dist = locategene(*ax, margin);
				   if (dist == 0) continue;
			   }
			   
			   
			   //Find chromosome
			   ichr = findchrm(*ax);
			   
			   //Do not allow cluster to cross a contig
			   if (*(ay-1) >= arm[ichr+1]) continue; 

			   //	printf("pos = %d, arm start at %d,annot index = %d\n", *ax, arm[ichr],ichr);
			   //	fprintf(save,"pos = %d, arm start at %d, annot index = %d\n", *ax,arm[ichr], ichr);
			   //Now ichr indexes into annotation name (annot_name) array
			   //as well as annotation contig arm start array (arm).
			   
			   for (j = 0; j < 26; j++) counts[j][0] = counts[j][1] = counts[j][2]= 0;
			   hitlist[0] = 0;  //start with empty hit list;
			   for (az = ax, pz = px, i=0; az < ay; az++, pz++){		
				   if (az != ax && *az < *(az-1) + strlen( *(pz-1))) continue;
				   hpt = hashQuery(mainHashTable, *az, strlen(*pz));
				   endcluster = *az+strlen(*pz)-1;
				   if (hpt) {
					   signature[i] = hpt->logical_name;
					   longsig[2*i] = signature[i];
					   i++;
					   counts[signature[i-1] - 'A'][0]++;
					   if (hpt->palindrome) {
						   longsig[2*i-1] = '*';
						   counts[signature[i-1] - 'A'][1]++;
						   counts[signature[i-1] - 'A'][2]++;
					   }
					   else if (hpt->flags) {
						   longsig[2*i - 1] = '+';
						   counts[signature[i-1] - 'A'][1]++;
					   }
					   else {
						   longsig[2*i - 1] = '-';
						   counts[signature[i-1] - 'A'][2]++;
					   }
					   sprintf(decnumber, "%c%d,", longsig[2*i - 1], *az - *ax);
					   strcat(hitlist, decnumber);

				   }
			   } 
			   signature[i] = 0;
			   longsig[2*i] = 0;
			   hitlist[strlen(hitlist)-1] = 0;  //remove last character
			   if (0 == newbool(counts, 26, boolean_expression)) {
				   continue;
			   }
			   if (count_only != 'C') printf("%s %s pos %d (length %d) in %s:%s\n", 
				   longsig, hitlist, (*ax) - arm[ichr], 
				   endcluster - (*ax) + 1, chromo_name[ichr], annot_name[ichr]);
			   if (count_only != 'C' && gene_names == 0) fprintf(save, "%s  %s pos %d (length %d) in %s:%s\n", 
				   longsig, hitlist, (*ax) - arm[ichr], 
				   endcluster - (*ax) + 1, chromo_name[ichr], annot_name[ichr]);
			   fprintf(archive, "%x %x %s\n", 
				   *ax, *(ay-1) + strlen(*(py-1)) - 1, longsig);
			   group_count++;
			   annotation_index = 
				   annotation_search(ichr, 
				   *ax - arm[ichr]/* + 1*/);
			   	//printf("annotation file %s (%d), starting index is %d\n",
			        //			annot_name[ichr], ichr, annotation_index);
			   
			   
			   if (annotation_index>= 0) {
				   //	printf("annotation index = %d\n", annotation_index);
				   //	fprintf( save, "annotation index = %d\n", annotation_index);
				   fullannot(annotation_index, 
					   *ax - arm[ichr] /*+ 1 *//*- map[annotation_index].start*/,
					   endcluster - arm[ichr]/* + 1*/ /*-map[annotation_index].start*/,
					   gene_filter, resultstring); 
				   if (count_only != 'C') {
					   printf("%s\n\n", resultstring);
					   if (gene_names) {
						   if (strlen(resultstring)) {
							 printallsynonyms(resultstring, syn_count);
							  // fprintf(save, "%s\n", resultstring);
						   }
					   }
					   else {
						   fprintf(save, "%s\n\n", resultstring);
					   }
				   }
			   }
			   else {
				  // full_results = 0;
			   }
			   if (full_results && gene_names == 0) {
				   fprintf(save, "          Window start: Annot %s, pos %d.", 
					   annot_name[ichr], *ax - arm[ichr] - ref + 1);
				   fprintf(save, "\n");
			   }
			   
			   if (full_results) {
				   //fetch entire cluster
				   clusterDNA = malloc(endcluster - (*ax) + 2);
				   for (cp = *ax, i = 0; cp <= endcluster; cp++, i++) {
					   clusterDNA[i] = fetchDNA(cp) | 0x20;  //force lower case
				   }
				   j = i;
				   clusterDNA[i] = 0;
				   for (az = ax, pz = px, aq = NULL; az < ay; az++, pz++) {
					   if (aq == NULL || *az >= *aq + strlen(*pq)) {
						   aq = az;
						   pq = pz;
						   for (cp = *az, i = 0; i < (int)strlen(*pz); cp++, i++) {
							   clusterDNA[cp - *ax] &= 0xDF;  //make upper case
						   }
					   }
				   }
				   
				   for (i = 0; i < j; i+=10) {
					   for (k = 0; k < 10; k++) {
						   line[k] = clusterDNA[i+k];
						   if (line[k] == 0) break;
					   }
					   line[k] = 0;

					   printf("%s ", line);
					   if (gene_names == 0) fprintf(save, "%s ", line);
					   
					   if (i%60 == 50) {
						   printf("\n");
						   if (gene_names == 0)fprintf(save, "\n");
					   }
				   }
				   printf("\n");
				   if (gene_names == 0) fprintf(save, "\n");
				   free(clusterDNA);
				   		   

			   }
			   ax = aq;
			   px = pq;
		   }
       }
   }
   
   printf("repeat = %d, Clusters found: %d\n", repeat, group_count);
   fprintf(save, "repeat = %d, Clusters found: %d\n", repeat, group_count);
   printf("expected %10.1f clusters, difference in std. deviations = %g\n",
      instances * pmore, (group_count - instances * pmore) / sigma);   
   fprintf(save, "expected %10.1f clusters, difference in std. deviations = %g\n",
      instances * pmore, (group_count - instances * pmore) / sigma);   
   fprintf(save, "======================================\n");
   fclose(archive);
}
   exit(0);

   printf("to save result, type name (without blanks) on the following line, otherwise press enter\n");
   readinput(line);
   if (strlen(line)) {
	   archive = fopen(line, "r");
	   if (archive) { // named file exists -- erase it before trying rename operation
		   fclose (archive);
		   remove(line);
	   }
	   rename("enhancerres.txt", line);
	   archive = 0;  //avoid trying to delete enhancerres.txt, which no longer exists
   }

   if(remote_input) return(0);
   
   goto next_grouping;
}

void makerevcompl(char *from, char *to, int num) {
   int i;
   for (i = 0; i < num; i++) {
      switch(from[i]) {
	 case 'A':
	 case 'a':
	     to[num - i - 1] = 'T';
	     break;
         case 'C':
	 case 'c':
	     to[num - i - 1] = 'G';
	     break;
         case 'G':
	 case 'g':
	     to[num - i - 1] = 'C';
	     break;
         case 'T':
	 case 't':
	     to[num - i - 1] = 'A';
	     break;
         default:
	     to[num - i - 1] = '0';
      }
      to[num] = 0;
   }
}

/***
void tally(unsigned int where, int patLen) {
	unsigned int cp, jptr, jmod;
	int i, j;
	hashPtr hpt;
	//Count left and right flanking nucleotides:
	hpt = hashQuery(mainHashTable, where, patLen);
	if (hpt && hpt->flags) {  //left-to-right direct reading
		cp = where - BORDERWIDTH;
		for (i = 0; i < 2 * BORDERWIDTH; i++) {
			if ( cp < DNALen) {
				jptr = cp / 4;
				jmod = cp & 3;
				j = (DNAstring[jptr]>>(6-2*jmod))&3;
				flank[i][j]++;
			}

			if (i != (BORDERWIDTH - 1)) cp++;
			else cp = where + patLen;
		}
	} else if(hpt) {
		cp = where + patLen + BORDERWIDTH - 1;
		for (i = 0; i < 2 * BORDERWIDTH; i++) {
			if (cp >= 0 && cp < DNALen) {
				jptr = cp / 4;
				jmod = cp & 3;
				j = 3 - ((DNAstring[jptr]>>(6-2*jmod))&3); //complement coding
				flank[i][j]++;
			}

			if (i != (BORDERWIDTH - 1)) cp--;
			else cp = where -1;
		}
	}
}
*/
/*
int findSpecific(char v[]) {
   int i, n, j;
   n = 0;
   while(v[n]) {
      switch (v[n]) {
	 case 'A': case 'C': case 'G': case 'T':
	    n++;
	    goto findmore;
	    break;
         default:
	    goto alldone;
      }
findmore:  continue;
   }
alldone:

   for (i = 0; i < (int)(DNALen - n); i++) {
      j = 0;
      while (v[j]) {
	 if (v[j] != DNAstring[i+j]) break;
	 j++;
      }
      if (j == n) return i;
   }
   return 0;
}
*/

void recursiveRegister(char *line, char* target[], int *index) {
      int i, n = strlen(line);
      int patnum = *index;
      char u[500];
      char *s, *t;
      int patLen;
      hashPtr hpt, hptrev;

      for (i = 0, t = line; i < n; i++) u[i] = *t++;    //make copy of input

      u[i] = 0;
      for (i = 0; i < n; i++) {
	 switch (u[i]) {
	    case 'B':            //C or G or T
	       u[i] = 'C';
	       recursiveRegister(u, target, index);
	       u[i] = 'G';
	       recursiveRegister(u, target, index);
	       u[i] = 'T';
	       recursiveRegister(u, target, index);
	       return;
	    case 'D':            //A or G or T
	       u[i] = 'A';
	       recursiveRegister(u, target, index);
	       u[i] = 'G';
	       recursiveRegister(u, target, index);
	       u[i] = 'T';
	       recursiveRegister(u, target, index);
	       return;
	    case 'H':            //A or C or T
	       u[i] = 'A';
	       recursiveRegister(u, target, index);
	       u[i] = 'C';
	       recursiveRegister(u, target, index);
	       u[i] = 'T';
	       recursiveRegister(u, target, index);
	       return;
	    case 'K':            //G or T
	       u[i] = 'G';
	       recursiveRegister(u, target, index);
	       u[i] = 'T';
	       recursiveRegister(u, target, index);
	       return;
	    case 'M':            //A or C
	       u[i] = 'A';
	       recursiveRegister(u, target, index);
	       u[i] = 'C';
	       recursiveRegister(u, target, index);
	       return;
	    case 'N':            //A or C or G or T
	       u[i] = 'A';
	       recursiveRegister(u, target, index);
	       u[i] = 'C';
	       recursiveRegister(u, target, index);
	       u[i] = 'G';
	       recursiveRegister(u, target, index);
	       u[i] = 'T';
	       recursiveRegister(u, target, index);
	       return;
	    case 'R':            //A or G
	       u[i] = 'A';
	       recursiveRegister(u, target, index);
	       u[i] = 'G';
	       recursiveRegister(u, target, index);
	       return;
	    case 'S':            //C or G
	       u[i] = 'C';
	       recursiveRegister(u, target, index);
	       u[i] = 'G';
	       recursiveRegister(u, target, index);
	       return;
            case 'U':
	       u[i] = 'T';       //allow U instead of T
	       recursiveRegister(u, target, index);
	       return;
	    case 'V':            //A or C or G
	       u[i] = 'A';
	       recursiveRegister(u, target, index);
	       u[i] = 'C';
	       recursiveRegister(u, target, index);
	       u[i] = 'G';
	       recursiveRegister(u, target, index);
	       return;
	    case 'W':            //A or T
	       u[i] = 'A';
	       recursiveRegister(u, target, index);
	       u[i] = 'T';
	       recursiveRegister(u, target, index);
	       return;
	    case 'Y':            //C or T
	       u[i] = 'C';
	       recursiveRegister(u, target, index);
	       u[i] = 'T';
	       recursiveRegister(u, target, index);
	       return;
	    default:
	       break;
	 }
      }

      target[patnum] = malloc(n+1);
      s = target[patnum]; t = line;
      while (*s++ = *t++);
      patLen = n;
      hpt = hashRegister(mainHashTable, target[patnum], patLen);
      if (hpt == NULL) return;    //We've already seen this pattern
      hpt->logical_name = pattern_id;
//      fprintf(save,
//	      "Target pattern is %s, of length %d\n", target[patnum],
//patLen);
      (*index)++;
      hpt->flags = 1;
	  hpt->palindrome = 0;
      target[patnum+1] = malloc(strlen(line)+1);
      makerevcompl(line, target[patnum+1], patLen);
      hptrev = hashRegister(mainHashTable, target[patnum+1], patLen);
      if (hptrev != NULL) {
		  hptrev->logical_name = pattern_id;
		  //         fprintf(save, "RevComp pattern is %s\n", target[patnum+1]);
		  (*index)++;
		  //printf("Register pointer = %x\n", hptrev);
		  hptrev->flags = 0;
		  hptrev->palindrome = 0;
      }
	  else hpt->palindrome = 1;
      return;
}

void readinput(char line[]) {

   int i;
   if (remote_input == 0) for (i=0; i<500; i++) {
      line[i] = getchar();
      if (line[i] == '\n' || line[i] == EOF) break;
   }
   else for (i = 0; i<500; i++) {
	   line[i] = getc(insource);
	   if (line[i] == '\n' || line[i] == EOF) break;
   }
   line[i] = 0;
}
void getLine(char line[], FILE *f) {
   int i;
   for (i=0; i<500; i++) {
      line[i] = getc(f);
      if (line[i] == '\n' || line[i] == EOF) break;
   }
   line[i] = 0;
}

int annotation_search(int chr, int index) {
   int i;
   
   for (i = 0; i < map_size; i++) {
      if (map[i].annot_number > chr)  break;
	  if (map[i].annot_number < chr) continue;
      if (map[i].start > index) continue;
      return i;
   }
   return -1;
}

void readgenome() {
	FILE* qnames;
	int i, counter;
	
	qnames = fopen("genome.txt", "rb");
	if (qnames == NULL) {
		printf("Genome does not exist -- quitting\n");
		exit(0);
	}
	
	for ( i = 0; i < (int)(DNALen/4+8); i++) {
		DNAstring[i] = fgetc(qnames);
		if(feof(qnames)) break;
	}
	
	counter = i;
	printf("compact genome took %u bytes (%u nucleotides)\n", counter, 4*counter);
	fclose (qnames);
	printf("expect %u nucleotides, used %u nucleotides for genome\n", DNALen+1, 4*i);
		
}

void insertClusters(unsigned int **av, char ***pv, char line[]) {
	int i;
	unsigned int j, *k;
	unsigned int start, stop;
	char *tempstring;
	hashPtr hpt;
	FILE *tempfile;

	k = *av;
	i = strlen(line);
	line[i-1] = 0;
//	printf("Line is: #%s#\n", line+1);
	tempfile = fopen (line+1, "r");
	if (tempfile == 0) {
	   printf("File %s cannot be opened; please respecify pattern\n", line);
	   return ;
	 }
	while(feof(tempfile) == 0) {
		getLine(line, tempfile);
		if (line[0] == 0) break;
		sscanf(line, "%x %x", &start, &stop);
		tempstring = malloc(stop - start + 2);
		if (stop/4 < search_start || start/4 > search_stop) {
			//cluster lies outside region of search
			continue;
		}
		for (j = start; j <= stop; j++) {
			tempstring[j-start] = fetchDNA((j));  //put cluster into string
		}
		tempstring[j-start] = 0;
		hpt = hashRegister(mainHashTable, tempstring, stop-start+1);
		**av = start;
		if (hpt) {
			hpt->logical_name = pattern_id;
			hpt->flags = 1;
			hpt->palindrome = 0;
		}
		else {
			hpt = hashQuery(mainHashTable, start, stop-start+1);
		}
		**pv = hpt->key;
/***
		**pv = malloc(stop - start + 2);
		for (j = start; j <= stop; j++) {
			(**pv)[j-start] = fetchDNA((j));  //put cluster into string
		}
		(**pv)[j-start] = 0;   //finish string
**/
		(*av)++;
		(*pv)++;
	}
	printf("File contained %d items\n", *av - k);
	fclose(tempfile);
}
int readContigs(FILE *qnames) {
	char line[1000], temp[100], temp2[100];
	int k = 0;
	int i, isize;
	while (!feof(qnames)) {
	   getLine(line, qnames);   //read a line of contig info
	   //Format of a contig line is name of contig file, chr. no, contig no, start of contig,
	   //size of contig.  We only use the start of contig for now.
	   if (line[0] == 0) break;
	   sscanf(line, "%s%s%d%d%d", temp, temp2, &i, &arm[k], &isize);
	   annot_name[k] = malloc(strlen(temp)+1);
	   chromo_name[k] = malloc(strlen(temp2)+1);
	   strcpy(annot_name[k], temp);
	   strcpy(chromo_name[k], temp2);
	   k++;
   }
	return k;}

int do_the_search(int start, int stop, unsigned int *av, char **pv, int slide_size,
				   unsigned int slide[], unsigned int patsize[]) {
	int i, ii, ix;
	unsigned char c;
	unsigned int rawhash = 0;
	hashPtr hpt;
	unsigned int *avold;

	avold = av;
   for (i = search_start; i < (int)search_stop; i++) {
	   c = DNAstring[i];
	   for (ii = 0; ii < 4; ii++) {
		   rawhash = (rawhash<<2) + (c>>6);
		   c<<=2;
		   if (i < 2) continue;
		   		   
		   for (ix = 0; ix < slide_size; ix++) {
			   hpt = hashFindAndListFast(mainHashTable, 4*i + ii - patsize[ix]+1,
				   patsize[ix], 
				   av, pv, rawhash & slide[ix]);
			   if (hpt) {
				   //found an instance!
				   //printf("i = %d, ptr = %x;, key = %s; addrptr = %x; vecptr = %x\n",
				   //          i, hpt, hpt->key, av, pv);
				   av++;
				   pv++;
			   }
			   //if ((i & 0xffff) == 0) printf("Char %d\n", i);
		   }
	   }
   }
   return (av - avold);

}

void printallsynonyms(char* resultstring, int syn_count){
	int i, j;
	fprintf(save, "%s\n", resultstring);
	for (i = 0; i < syn_count; i++) {
                if (strlen(gene_name[i]) < 3) continue;
		if ( strstr(resultstring, gene_name[i]) ) {
//                      fprintf(save, "index in xallnames = %d, gene number = %d\n", i, gene_number[i]);
			for (j = i+1; j < syn_count; j++) {
				if (gene_number[i] != gene_number[j]) break;
				fprintf(save, "%s\n", gene_name[j]);
			}
		}
	}
	return;
}

void readgenenames(FILE *allnames, int syn_count){

	int i, n;
	char line[500], line2[1000];

	for (i = 0; i < syn_count; i++) {
		getLine(line2, allnames);
		sscanf(line2, "%[^,],%d", line, &gene_number[i]);
		n = strlen(line);
		gene_name[i] = malloc(n+1);
		strcpy(gene_name[i], line);
	}
}

static int seed = 1;
int Random() {
   seed = seed * 3125;
   return(seed);
}
