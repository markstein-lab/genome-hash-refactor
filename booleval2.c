//Copyright (c) in silico Labs, LLC, 2001 
#include <stdlib.h>
//#include <malloc.h>
#include "group_defs2.h"
#include <string.h>

int booleval(unsigned char *boolstr, unsigned char * clusseq) {
unsigned char  var[100];
int i, j, k, m, lparen, rparen, mult, blen;
//#include <ctype.h>

   //quickie pseudo boolean

   blen = strlen(boolstr);
   j = 0;
   k=  0;
   lparen = 0;
   rparen = 0;
   mult = 1;
  
   for (i = 0; i<blen; i++) {
     if (isalpha (boolstr[i]) && lparen == 1) {
       var[j] = toupper(boolstr[i]); // seq id 
       j++;
     }
    
     else if (boolstr[i] == '(' ) {
       lparen = 1;
     }
     
     else if (boolstr[i] == ')' && lparen == 1) {
     // check requirements against cluster 
       for (m=0; m<j && mult>=1; m++) {
         for (k=0; k<(int)strlen(clusseq) && mult>=1; k++) {
           if (clusseq[k] == var[m]) {
             mult--;
           }
         }
       }
       if (mult <= 0) {
         mult = 1;
         j=0;
         
       }
       else {   
         
         return (0);
       }
     }    
      
     else if (isdigit (boolstr[i])) {
       // form integer multiplier  num = boolstr[i] - '0' ;  // pick up units
       mult = boolstr[i] - '0';
       
       while (i < blen) {
          if (isdigit (boolstr[i+1])) {
            mult = mult*10 + boolstr[i+1] - '0';
            i++;
          }
          else break;
       }
            
     }
    }
    
    return(1); 
  }  
     
   
 
