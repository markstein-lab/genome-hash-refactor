// copyright (c) in silico Labs, LLC 1998
void hsort2(unsigned int x[], unsigned int y[], int n)
{	//heap sort
	//essentially a sorting algorithm over 1-based indexing.
	//for that reason, all subscripts are reduced by 1 to conform to
	//     conventions of the C language (zero-based indexing)
	//In this version, the sort is of the items in x[]. The item y[i] is
	//     associated with x[i], so whatever rearrangement is performed on the
	//     x[] array, the same is performed on y[].
	int i, j, k, ii;
	   unsigned int temp, temp2;
       for (ii=2; ii<=n; ii++){
		   i = ii;
		   j = i/2;
		   temp = x[i-1];
		   temp2 = y[i-1];
		   while ((j > 0) && (temp> x[j-1])) {
			   x[i-1]=x[j-1];
			   y[i-1]=y[j-1];
			   i = j;
			   j=i/2;
           };
		   x[i-1] = temp;
		   y[i-1] = temp2;
	   };
	   
       for (i = n; i > 1; i--) {
		   temp = x[i-1];
		   temp2= y[i-1];
		   x[i-1] = x[0];
		   y[i-1] = y[0];
		   j = 1;
		   k = 2;
		   while (k < i) {
			   if ((k+1 < i)) 
				   if (x[k] > x[k-1]) k = k + 1;
				   if (x[k-1] > temp){
					   x[j-1] = x[k-1];
					   y[j-1] = y[k-1];
					   j = k;
					   k=k+k;
				   }
				   else break;
		   }
           x[j-1] = temp;
		   y[j-1] = temp2;
	   }
}

