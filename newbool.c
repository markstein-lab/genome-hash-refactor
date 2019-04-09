//Copyright (c) in silico Labs, LLC 2004

// Reads input from files named in filelist.txt. Each file listed in filelist
// contains a sequence of nucleotides.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define NUMBER 1
#define NOT	   2
#define HITS   3
#define BOOL   4
#define AND	   5
#define OR	   6
#define STOP   7
#define START  8
#define COLUMNS 9
#define LIMIT  10
#define QUANTUM 11
#define TRUE   1
#define FALSE  0

int newbool(int counts[][3], int arraysize, char boolexp[])
{
	int sum, i, cols, pos;
	int stack[20], oper[20], stkndx, operndx, number;
	unsigned int mask;
	int skip, maxlimit;
	char c;
	int quantum;
	int oneDirection, psum, msum, pcols, mcols;

	i = pos = 0;
	stkndx = 0;
	operndx = 0;
	stack[0] = 1;
	skip = 0;
	quantum = 1;
	while(1) {
		if (operndx > 0 && !skip) { // try to collapse stack
			if (operndx > 1 && oper[operndx-1] == STOP) { // START-STOP -- remove parens
				if (oper[operndx-2] == START) {
					operndx = operndx - 2;
					continue;
				}

				if (oper[operndx-2] == NOT) { // complement top of stack
					stack[stkndx-1] = stack[stkndx-1]? 0 : 1;
					oper[operndx-2] = STOP;
					operndx = operndx -1;
					continue;
				}
				if (oper[operndx-2] == AND) { // AND top of stack items
					stack[stkndx-2] &= stack[stkndx-1];
					oper[operndx-2] = STOP;
					stkndx -= 1;
					operndx -= 1;
					continue;
				}

				if (oper[operndx-2] == OR) { //	OR top of stack items
					stack[stkndx-2] |= stack[stkndx-1];
					oper[operndx-2] = STOP;
					stkndx -= 1;
					operndx -= 1;
					continue;
				}

			}
			if (oper[operndx-1] == NOT) { // NOT - complement top of stack
				stack[stkndx-1] = stack[stkndx-1]? 0 : 1;
				operndx -= 1;
				continue;
			}

			if (oper[operndx-1] == AND) {
				stack[stkndx-2] &= stack[stkndx-1];
				stkndx--;
				operndx--;
				continue;
			}

			if (oper[operndx-1] == OR) {
				if (boolexp[pos] == 0) {
					stack[stkndx-2] |= stack[stkndx-1];
					stkndx--;
					operndx--;
					continue;
				}
			}
		}

		c = boolexp[pos++];
		if (c == 0) {
			skip = 0;
			if (operndx> 0) {
				pos--;
				continue;
			}
			if (operndx < 0) {
				// illegal expression -- return false
				stack[0] = 0;
			}
			return stack[0];
		}
		if (c == ' ') continue;
		skip = 0;
		if (isdigit(c)) { // interpret integer
			number = c - '0';
			while (isdigit(boolexp[pos])) {
				number = 10*number + boolexp[pos++] - '0';
			}
			stack[stkndx++] = number;
			oper[operndx++] = NUMBER;
			continue;
		}
		if (isupper(c)) { // interpret class letter
			mask = 1<<(c-'A');
			oneDirection = FALSE;
			while (isupper(boolexp[pos])) {
				mask |= 1<<(boolexp[pos++] - 'A');
			}
			if (boolexp[pos] == '=') {
				oneDirection = TRUE;
				pos++;
			}
			while (boolexp[pos] == ' ') pos++;
			if (operndx > 0 && oper[operndx-1] == START &&
				boolexp[pos] == ')' ) {
				// parens surround a string of letters -- drop them
				pos++; // drop the input right paren
				operndx--;
			}
			sum = 0; cols = 0; maxlimit = 0;
			psum = 0; pcols = 0;
			msum = 0; mcols = 0;
			for (i = 0; i < arraysize; i++) {
				if (mask & (1<<i)) {
					if (oneDirection == FALSE) {
						sum += counts[i][0];
						if (counts[i][0] >= quantum) cols += 1;
					}
					else {
						psum += counts[i][1];
						msum += counts[i][2];
						if (counts[i][1] >= quantum) pcols += 1;
						if (counts[i][2] >= quantum) mcols += 1;
					}
				}
				else {
					if (oneDirection == FALSE) {
						if (counts[i][0] > maxlimit) maxlimit = counts[i][0];
					}
					else {
						if (counts[i][1] > maxlimit) maxlimit = counts[i][1];
						if (counts[i][2] > maxlimit) maxlimit = counts[i][2];
					}
				}
			}
			if (oneDirection) {
				sum = (psum > msum) ? psum : msum;
				cols = (pcols > mcols) ? pcols : mcols;
			}
			if (operndx > 0 && oper[operndx-1] == COLUMNS) {
			   stack[stkndx++] = cols;
			   operndx--;
			}
			else if (operndx > 0 && oper[operndx-1] == LIMIT) {
			   stack[stkndx++] = maxlimit;
			   operndx--;
			}
			else {
			   stack[stkndx++] = sum;
			}
			oper[operndx++] = HITS;
			if (operndx > 1 && oper[operndx-2] == NUMBER) {
				// specific number of hits wanted.
				if (stack[stkndx-1] >= stack[stkndx-2]) stack[stkndx-2] = 1;
				else stack[stkndx-2] = 0;
				stkndx--;
				operndx-=2;
			}
			else {
				if (stack[stkndx-1] > 0) stack[stkndx-1] = 1;
				operndx-=1;
			}
			quantum = 1;
			continue;
		}

		if (c == ')') {
			oper[operndx++] = STOP;
			continue;
		}


		if (c == '(') { // parenthesized expression
			oper[operndx++] = START;
			skip = 1;
			continue;
		}

		if (c == '|') {
			oper[operndx++] = OR;
			skip = 1;
			continue;
		}

		if (c == '&') {
			oper[operndx++] = AND;
			skip = 1;
			continue;
		}

		if (c == '~') {
			oper[operndx++] = NOT;
			skip = 1;
			continue;
		}

		if (c == '#') {
			if (oper[operndx] == START) return 0; // error condition -return false
			oper[operndx++] = COLUMNS;
			continue;
		}

		if (c == 'n' && boolexp[pos] == 'o' && boolexp[pos+1] == 't') {
			pos +=2;
			oper[operndx++] = NOT;
			skip = 1;
			continue;
		}

		if (c == 'a' && boolexp[pos] == 'n' && boolexp[pos+1] == 'd') {
			pos +=2;
			oper[operndx++] = AND;
			skip = 1;
			continue;
		}

		if (c == 'o' && boolexp[pos] == 'r') {
			pos +=1;
			oper[operndx++] = OR;
			skip = 1;
			continue;
		}

		if (c == '$') {
			if (oper[operndx] == START) return 0; // error condition -return false
			oper[operndx++] = LIMIT;
			continue;
		}

		if (c == '*') {
			if (oper[operndx] == START) return 0; // error condition -return false
			quantum = stack[--stkndx];
			operndx--;
			continue;
		}
	}
}
