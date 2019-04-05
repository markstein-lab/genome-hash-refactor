
//Copyright (c) World Internet Productions, LLC, 1999
//Copyright (c) in silico Labs, LLC 2006

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
unsigned char * searchstr ( unsigned char *, unsigned char *);
void getLine (char* line, FILE *stream);
unsigned char * skipwhitesp(unsigned char * string);
void add(unsigned char);
#define wordsize sizeof(unsigned int)
int j, k, excount;
unsigned int counter, oldcounter;
unsigned char genome;
unsigned int exceptions;
FILE *dnadata;
FILE *exceptionfile;
FILE *allnames;
typedef struct synonym synonym;
struct synonym{
	char *name;
	int  mainname;
} ;
int total_genes;
int max_genes;
unsigned int fileSize;
void add_name(char *, int);
char* readlimits(char* ptr, int*, int*);

   FILE *master;

synonym *synonyms;

char* readlimits(char* ptr, int* start, int* stop) {
	char *memory; 
	memory = ptr;
	if (*ptr == '<') ptr++;    //move past "<" if present;
	sscanf(ptr, "%d..", start);
	ptr = searchstr(ptr, "..");
	if (ptr == 0) {
//		printf("%s", memory);
		*stop = 0;
		return ptr;
	}
	if (*ptr == '>') ptr++;    //move past ">" if present
	sscanf(ptr, "%d", stop);
        ptr = skipwhitesp(ptr);    //avoid white space
        while (*ptr >= '0' && *ptr <= '9') ptr++;  //skip over number
	return ptr;
}

int end_of_sequence(char *line, FILE *stream, char direction);

int main( void )
{

   FILE *stream;
   FILE *names;
   FILE *data;
   FILE *UTRs;
   char line[500], line2[500], head[500];
   char annot_name[20];
   char junk[100];
   char gene_name[500];
   char direction;
   char temp_direction;
   char *ptr, *ptr2, *ptr3;
   char chrm_name[500];
   char full_name[500];
   char first_name[500];
   char file_name[500];
  
   int i, z;
   unsigned int origin;
   int start, stop, size;
   int annot_id;
   int number_of_genes;

   int number_of_mrnas, number_of_cds;
   int mrna[1000], cds[1000], sizeUTR[1000]; // ends of mrna, cds, respectively
   int flag;
   
   stream = NULL;
   dnadata = fopen("genome.txt", "wb");
   UTRs = fopen("utrinfo.txt", "w");
   master = fopen("master.txt", "r");
   exceptionfile = fopen("exceptions.txt", "w");
   if (master == NULL) {
	   printf("file master.txt, containing chromosome names, is missing!\n");
	   exit(101);
   }
   annot_id = 0;
   origin = 0;
   number_of_genes = 0;
   excount = 0;
   counter = 0;
   oldcounter = 0;
   names = fopen("xcontigs.txt", "w");
   data = fopen("xchrdata.txt", "w");
   j = 0;	//number of bits in	genome char. buffer
   k = 0;   //number of bits in exception word buffer
   genome = 0; //genome buffer
   fileSize = 0;  //bytes written to genome.txt
   max_genes = 100000;
   synonyms = calloc (max_genes, sizeof(synonym));
   total_genes = 0;
   while (1) {
	   line[0] = 0;    //initialize line to NULL string
	   fscanf(master, "%s", line);
	   if (strlen(line) == 0) break;  //no more chromosomes
	   strcpy(chrm_name, line);
	   strcpy(full_name, line);	  //chromosome name for the record
	   strcat(chrm_name, ".gbk");
	   strcpy(file_name, "chr");
	   strcat(file_name, chrm_name);  //construct name chrxx.gbk
	   stream = fopen(file_name, "r");
	   if (stream == NULL) {
		   printf("attempt to open file %s produced %d\n", file_name, stream);
		   exit(100);
	   }
	   else printf("Opened input file %s\n", file_name);
	   
	   if (dnadata == NULL) {
		   printf("Attempt to open file %s for writing failed\n",
			   "genome.txt");
	   }
	   
	   printf ("Chromosome %s starts at position %u\n", full_name, origin);
	   strcat(line, ".txt");
	   strcpy(file_name, "contigs");
	   strcat(file_name, line);
//	   names = fopen(file_name,"w");
	   
	   strcpy(file_name, "chrdata");
	   strcat(file_name, line);
//	   data = fopen(file_name, "w");
next_annotation:
	   getLine(head, stream);
	   if (strstr(head,"LOCUS") != head) {
		   if (feof(stream)) goto finish;
		   goto next_annotation;
	   }
find_definition:
	   getLine(line, stream);
	   if (strstr(line, "DEFINITION") != line) {
		   if (feof(stream)) goto finish;
		   goto find_definition;
	   }
	   if (strstr(line, "alternate assembly")) {
	       goto next_annotation;
	   }
   	   
find_annot_id:
	   //get name of the annotation file.
	   getLine(line, stream);
	   if (strstr(line,"VERSION") != line) goto find_annot_id;
	   sscanf(line, "VERSION %s", annot_name);
	   
find_size:
	   getLine(line, stream);
	   if (!strstr(line, "source")) goto find_size;
           readlimits(line, &start, &size);
	   //sscanf(line, "%s%d..%d", junk, &start, &size);
	   
	   fprintf(names, "%s %s %d %u %u\n", annot_name, full_name, 
		   annot_id, origin, size);
#ifdef DEBUG
	   printf("%s %s %s, size = %u\n", head, line, annot_name, size);
#endif
	   
find_features:
	   //This section finds all the genes in the annotation file, writes 
           //out the start, stop, and direction of all the genes, together with
           //the annotation file index. Finally, when it finds the DNA sequence 
           //represented by the annotation file, the DNA sequence is written to 
           //a file whose name is the name of the annotation sequence.
	   
	   getLine(line, stream);
	   sscanf(line, "%s", junk);
	   if (strcmp(junk, "gene") == 0) {
		   //process a gene
		   ptr = searchstr(line, " gene ");
		   if (ptr == 0) goto find_features;  // false alarm?
		   ptr = skipwhitesp(ptr);    //find start of gene limits
		   if (*ptr == 'c') {
			   //reverse complement -- indicate in table
			   temp_direction = '-';  //may not be a gene
			   ptr = searchstr(ptr, "(");
		   }
		   else {
			   temp_direction = '+';  //may not be a gene
		   }
		   readlimits(ptr, &start, &stop);
		   if (stop == 0) goto find_features;   //badly formed data -- ignore it
//		   sscanf(ptr, "%d..%d", &start, &stop);
		   
                   //looks like a new gene entry has been found
                   //keep the new direction
                   direction = temp_direction;

		   //Now go to next line to get the gene name
		   getLine(line, stream);
		   ptr = searchstr(line, "gene=\"");
		   //If gene name field is found, enter it into list of names
		   if (ptr) {
			   ptr3 = ptr;  //remember start of gene name
			   add_name(ptr, number_of_genes);
			   getLine(line2, stream);
			   ptr = searchstr(line2, "locus_tag=\"");
		   }
		   else {
			   //if gene does not have a gene name, it might
			   //have a locus_tag identifier, which we use instead
			   ptr = searchstr(line, "locus_tag=\"");
			   if (ptr == 0) goto find_features;  //false alarm?
			   ptr3 = ptr;  //remember start of gene name
		   }
           if (ptr) {
			   add_name(ptr, number_of_genes);  //enter name onto lisT
		   }
		   ptr2 = gene_name;
		   while(*ptr3 != '"') {
			   *ptr2++ = *ptr3++;
		   }
		   *ptr2 = 0;
		   fprintf(data, "%d %u %u %c %s\n", 
			   annot_id, start, stop, direction, gene_name);
                   strcpy(first_name, gene_name);
		   //Now look to see if there are more aliases for our gene
		   
		   getLine(line, stream);
		   ptr = searchstr(line, "synonym=\"");
		   if (ptr == 0) { //look for alternate form: synonym part of a note field
			   ptr = searchstr(line, "synonym");
			   if (ptr) {
				   if (*ptr == 's') ptr++;
				   if (*ptr == ':') ptr++; else ptr = 0;
			   }
		   }
		   if (ptr) {
			   //extract all the gene names and insert into name list
			   do {
				   while (*ptr <= ' ') {
					   if (*ptr) ptr++;
					   else {
						   for (i = 0; i < 500; i++) line[i] = 0;
						   getLine(line, stream);
						   ptr = line;
					   }
				   }  
				   sscanf(ptr, "%497[^\r\n,;\"]", gene_name); 
				   add_name(gene_name, number_of_genes);
				   ptr += strlen(gene_name);
			   } while (*ptr++ != '"'); 
		   } 
		   
		   number_of_genes++;
           number_of_mrnas = 0;
           number_of_cds = 0;
	   }
	   else if (strcmp(junk, "mRNA") == 0) {
		   if (number_of_mrnas >=0) {
               mrna[number_of_mrnas] = end_of_sequence(line, stream, direction);
               number_of_mrnas++;  //increase number of mRNAs found; 
		   }
	   }
	   else if (strcmp(junk, "CDS") == 0) {
		   if (number_of_cds >= 0) {
               cds[number_of_cds] = end_of_sequence(line, stream, direction);
               number_of_cds++;
               if (number_of_cds == number_of_mrnas) {
				   //we are ready to look at the sizes of the UTRs.
				   for (i = 0; i < number_of_mrnas; i++) {
					   //examine each splice variant
					   sizeUTR[i] = direction == '+' ? mrna[i] - cds[i] 
						   : cds[i] - mrna[i];
					   if (sizeUTR[i] <= 0) continue; //do nothing for no UTR
                                           if (sizeUTR[i] < 0) {
                                               printf("warning - negative UTR size\n");
                                           }
					   flag = 0;
					   for (z = 0; z < i; z++) {
						   if (sizeUTR[i] == sizeUTR[z] ) flag = 1;
					   }
					   if (flag) continue; //not a unique UTR
					   //we have a unique UTR: output it
					   if (direction == '+') {
						   fprintf(UTRs, "%d %u %d %c %s %s\n", 
							   number_of_genes -1, cds[i] + 1 + origin, sizeUTR[i], direction, full_name, first_name);
					   }
					   else {
						   fprintf(UTRs, "%d %u %d %c %s %s\n", 
							   number_of_genes -1, cds[i] - 1 + origin, sizeUTR[i], direction, full_name, first_name);
					   }
				   }
				   number_of_cds = -1;
				   number_of_mrnas = -1;				   
               }						
		   }
	   }
	   else if (strcmp(junk, "ORIGIN") == 0) {
		   //process DNA string;
		   //start by inserting an N into the genome
		   add('N');
		   while (1) {
			   getLine(line, stream);
			   if (line) for (i = 0; line[i]; i++) if (line[i] != ' ') break;
			   if (line[i] == '/' || feof(stream) )
				   break;  //found end of data for annot. file
			   for (i = 0; line[i]; i++) {
				   if (line[i] >= 'A') {
					   add(line[i]);
				   }
			   }
		   }
		   
                   if (counter/4 != fileSize) {
                      printf("Counter is %u, origin = %u, number of bytes should be %u, number of bytes written is %u\n", counter, origin+size+1, counter/4, fileSize);
                   }
		   if (counter - oldcounter != (unsigned int)size +1) {
			   printf("%s %s %s, size = %u\n", 
				   head, line, annot_name, size);
			   printf("size should be %u, counted %u\n",
				   size, counter-oldcounter-1);
		   }
		   oldcounter = counter;
		   goto finish_features;
	   }
	   goto find_features;
	   
finish_features:
	   annot_id++;
	   origin = origin + size + 1;
	   goto next_annotation;
	   
finish:
	   
	   printf("found %d annotation files for chromosome %s\n", annot_id, full_name);
	   add ('N');    //end chromosome with an 'N'
	   oldcounter = counter;    //update origin of next item to read
	   origin = origin + 1;    //account for the extra 'N'

//	   fclose(stream);
   }
   //Fill out exception word with N's
   while(k) (add ('N'));
   add(0); //make sure dnastring looks like a string, i.e. ends in a 0.
//   fputc(0, dnadata);   
   fclose(dnadata);
   fclose(names);
   fclose(data);
   fclose(master);
   fclose (exceptionfile);
   allnames = fopen("xallnames.txt", "w");
   for (i = 0; i < total_genes; i++) {
	   fprintf(allnames, "%s, %d\n", synonyms[i].name, synonyms[i].mainname);
   }
   fclose(allnames);
   
   master = fopen("datasize.txt", "w");
   fprintf(master, "%u %d %d %d %d\n", origin, number_of_genes, annot_id, excount,
	   total_genes);
   printf("origin = %u, counter = %u, number of genes = %d, number of annotation files = %d\n",
	   origin, counter, number_of_genes, annot_id);
   printf("number of exception words = %d, number of gene names = %d\n", 
	   excount, total_genes);
   fclose(master);
   exit(100);
}

void getLine (char* line, FILE *stream) {
   int i;
   for (i = 0; i < 500; i++) {
      line[i] = getc(stream);
      if (line[i] == '\n' || line[i] == EOF || line[i] == 0) break;
   }
   line[i] = 0;
}
void add(unsigned char x) {
	int y;
	y = x & 0xDF;   //make upper case
	switch (y) {
	case 'A':
		genome <<=2;
		exceptions <<=1;
		break;
	case 'C':
		genome = 1 + (genome <<2);
		exceptions <<=1;
		break;
	case 'G':
		genome = 2 + (genome << 2);
		exceptions <<= 1;
		break;
	case 'T':
		genome = 3 + (genome << 2);
		exceptions <<= 1;
		break;
	default:
		genome <<= 2;	//encode invalid characters like 'A's
		exceptions = 1 + (exceptions << 1);  //insert a bit into exception vector
		break;
	}
	j++;
	if (j ==4) {   //when a genome nibble is full
		fputc(genome, dnadata);	   //write it out
                fileSize++;   //add to size of dnadata file
		j = 0;
		genome = 0;
	}
	k++;
	if (k == 32) {  //a potential word of exceptions
		if (exceptions) { //if there are exception bits
			fprintf(exceptionfile, "%d %8.8llx\n", counter/32, exceptions);
			excount++; //count number of exception words
		}
		k = 0;
		exceptions = 0;
	}
	counter++;
}

void add_name (char *name, int index) {
	unsigned int i;
	char temp[500];
        if (strlen(name) > 499) {
           printf("too long gene name: %s\n", name);
           exit(100);
        }
        if (total_genes >= max_genes) {
           max_genes += 100000;
           synonyms = realloc(synonyms, max_genes*sizeof(synonym));
#ifdef DEBUG
           printf("total number of genes names seen so far is %d\n",
             total_genes);
#endif
        }
	for (i = 0; i < strlen(name); i++) {
		if (name[i] && name[i] != '"') temp[i] = name[i];
		else break;
	}
	temp[i] = 0;

	synonyms[total_genes].name = malloc(strlen(temp) + 1);
	strcpy(synonyms[total_genes].name, temp);
	synonyms[total_genes].mainname = index;
	total_genes++;
}

int end_of_sequence(char *line, FILE *stream, char direction) {
   int k;
   char *ptr, *ptr2;
   char junk[200];
   int start, stop;

   ptr = skipwhitesp(line);
   sscanf(ptr, "%s", junk);
   k = strlen(junk);
   ptr += k;	//point to end of string on mRNA or CDS statement
   ptr = skipwhitesp(ptr);   //skip over whitespace
   if (direction == '-') {  //reading in negative direction
      while( *ptr < '0' || *ptr > '9') *ptr++;  //skip to digits
      readlimits(ptr, &start, &stop);
      return(start);
   }
   else {
	  if (0 == searchstr(ptr, "(")) { //simple limits, no complements or joins
		  readlimits(ptr, &start, &stop);
		  return(stop);
	  }
      while (0 == searchstr(ptr, ")" ) ) {
		  //skip over lines that do not end the join or complement
		  getLine(ptr, stream);
	  }
	  while (ptr2 = searchstr(ptr, ",")) ptr = ptr2;
	  //now we are pointing to after the last comma
	  ptr = skipwhitesp(ptr); //eliminate any white space
      readlimits(ptr, &start, &stop);
      return(stop);
   } 
}
		
