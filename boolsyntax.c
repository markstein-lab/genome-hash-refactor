//Copyright (c) in silico Labs, LLC 2005

// Checks the syntax of a "boolean expression"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define DIGITS	 1
#define MULTVAR	 2
#define LPAREN	 3
#define EXPR	 4
#define AND		 5
#define OR		 6
#define NOT		 7
#define STAR	 8
#define SHARP	 9
#define DOLLAR	 10

char* digits(char* message, int i)
{
		while(isdigit(message[i])) i++;
		return message+i-1;
}

char* letters(char* message, int i)
{
		while(isupper(message[i])) i++;
		return message+i-1;
}

int collapse(int stack[], int stkptr)
{
		while (stkptr) {
				if (stkptr == 1) return stkptr;
				if (stack[stkptr-1] == EXPR) {
						if (stack[stkptr-2] == NOT) {
								stack[stkptr - 2] = EXPR;
								stkptr--;
								continue;
						}
						if (stkptr > 2) {
								if (stack[stkptr-3] == EXPR &&
									(stack[stkptr-2] == AND ||
									 stack[stkptr-2] == OR)) {
										stkptr-=2;
										continue;
								}
						}
				}
				return stkptr;
		}
		return stkptr;
}

char* booleansyntax(char message[])
{
		int lparen, stkptr;
		int stack[100];
		int i;
		char* k;

		lparen = 0; // set to true if an open lparen is seen
		stkptr = 0;
		if (strlen(message) == 0) return message;

		for (i = 0; i < strlen(message); i++) {
				if (isspace(message[i])) continue;
				if (isdigit(message[i])) {
						k = digits(message, i);
						stack[stkptr++] = DIGITS;
						i = k - message;
						continue;
				}
				if (isupper(message[i])) {
						k = letters(message, i);
						i = k - message;
						if(message[i+1] == '=') {
								// require all hits to be in same direction
								i++;
						}
						stack[stkptr++] = MULTVAR;
						if (stkptr == 1) stack[stkptr-1] = EXPR;
						if(stkptr > 1 && stack[stkptr-2] == DIGITS) {
								stack[stkptr-2] = EXPR;
								stkptr--;
						}

						if (stkptr > 2 && stack[stkptr-2] == SHARP ) {
								if (stkptr > 3 && stack[stkptr-3] == DIGITS) {
										if (stkptr > 4 && stack[stkptr-4] == STAR &&
											stack[stkptr-5] == DIGITS) {
												stack[stkptr - 5] = EXPR;
												stkptr -= 4;
										}
										else {
												stack[stkptr - 3] = EXPR;
												stkptr-=2;
										}
								}
								else if (stack[stkptr-3] == DIGITS) {
										stack[stkptr - 3] = EXPR;
										stkptr-=2;
								}
								else {
										stack[stkptr - 2] = EXPR;
										stkptr -= 1;
								}
						}
						else if (stkptr > 1 && stack[stkptr - 2] == SHARP) {
								stack[stkptr - 2] = EXPR;
								stkptr--;
						}
						if (stkptr > 1 && stack[stkptr-2] == NOT) {
								stack[stkptr-2] = EXPR;
								stkptr--;
						}
						if (stkptr > 2 && ( // try to collapse expr op expr into expr
									stack[stkptr-2] == AND ||
									stack[stkptr-2] == OR ) &&
							stack[stkptr-3] == EXPR) {
								stkptr-=2;
						}
						continue;
				}
				if (message[i] == '(') {
						k = booleansyntax(message+i+1);
						if (k == 0) return(0); // improper subexpression- total expression bad
						stack[stkptr++] = EXPR; // paren expression valid. put EXPR onto stack
						i = k - message;
						continue;
				}
				if (message[i] == ')') {
						stkptr = collapse (stack, stkptr);
						if (lparen) {
								if (stkptr > 0 && stack[stkptr-1] == LPAREN) {
										stack[stkptr-1] = EXPR;
										lparen = 0;
										continue;
								}
								return 0;
						}
						if (stkptr == 1) return message+i;
						else return 0;
				}
				if (message[i] == '*') {
						stack[stkptr++] = STAR;
						continue;
				}
				if (message[i] == '#') {
						stack[stkptr++] = SHARP;
						continue;
				}
				if (message[i] == '&') {
						stack[stkptr++] = AND;
						continue;
				}
				if (message[i] == '|') {
						stack[stkptr++] = OR;
						continue;
				}
				if (message[i] == '~') {
						stack[stkptr++] = NOT;
						continue;
				}
				if (message[i] == '$') {
						stack[stkptr++] = DOLLAR;
						continue;
				}
				if (message[i] == 'a' && message[i+1] == 'n' && message[i+2] == 'd') {
						stack[stkptr++] = AND;
						i+=2;
						continue;
				}
				if (message[i] == 'n' && message[i+1] == 'o' && message[i+2] == 't') {
						stack[stkptr++] = NOT;
						i+=2;
						continue;
				}
				if (message[i] == 'o' && message[i+1] == 'r') {
						stack[stkptr++] = OR;
						i+=1;
						continue;
				}
				return 0; // improper character in expression
		}
		stkptr = collapse(stack, stkptr);
		if (stkptr == 1 && stack[stkptr-1] == EXPR) return message+i-1;
		else return 0;
}
