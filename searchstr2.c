//Copyright (c) World Internet Productions, LLC, 1999
//Copyright (c) in silico Labs, LLC 2006
#include <stdlib.h>

unsigned char *searchstr(unsigned char *string, unsigned char *pattern)
{
		unsigned char *ptr, *s, *t;
		for (ptr = string; *ptr ; ptr++) {
				if (*ptr != *pattern) continue;
				s = ptr + 1;
				t = pattern + 1;
				while (*t) {
						if (*s != *t) goto nextposition;
						s++; t++;
				}
				return (s);
		nextposition:
				;
		}
		return (unsigned char *)0;
}

unsigned char *skipwhitesp(unsigned char *string)
{
		unsigned char *ptr, c;
		ptr = string;
		while (c = *ptr) {
				if (c > ' ') return ptr;
				ptr++;
		}
		return ptr;
}

unsigned char *findattrib(unsigned char *string, unsigned char *pattern)
{
		unsigned char *c;
		c = searchstr(string, pattern);
		if (c) { // found the attribute
				c = searchstr (c,(unsigned char *) ">"); // skip over closing attribute bracket
				if(c) c = skipwhitesp(c); // skip white space, point to date!
		}
		return c;
}

unsigned char *getstrattrib(unsigned char *string,
							unsigned char *pattern,
							unsigned char *result)
{
		// Finds the attribute name "pattern" within the string "string", and then
		// copies the attribute associated with the name into "result".
		// The function returns a pointer to the first character in the
		// string following the attribute.
		unsigned char *c, *s;
		s = result;
		c = findattrib(string, pattern);
		while(*c > ' ') {
				*s++ = *c++;
		}
		*s = 0; // insert end-of-string character
		return c;

}
