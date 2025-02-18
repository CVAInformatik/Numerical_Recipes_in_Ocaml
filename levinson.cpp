

#include <stdlib.h>
#include <stdio.h>

typedef unsigned int u32;

void mulAndPrint(char *s,  int N,  double *r, double *v)
{
   for(u32 ir = 0 ; ir< N; ir++)
        {
        	double vs = 0.0;
        	for(u32 ibf = 0; ibf < N; ibf++){
        	    vs+= (r[N-1 + ir-ibf]*v[ibf]);
        	}
        	printf( " %s : (r*v)[%d] : %10.4lf \n",s, ir, vs);
        }
}

void  dump(int N,  double *r, double *y, double *x, double *f, double *b)
{
	printf( "r[] "); for(u32 ix = 0 ; ix < ((2*N)-1); ix++)  printf(" %10.4lf", r[ix]); printf("\n");
	printf( "y[] "); for(u32 ix = 0 ; ix < N; ix++)  printf(" %10.4lf", y[ix]); printf("\n");
	printf( "x[] "); for(u32 ix = 0 ; ix < N; ix++)  printf(" %10.4lf", x[ix]); printf("\n");
	printf( "f[] "); for(u32 ix = 0 ; ix < N; ix++)  printf(" %10.4lf", f[ix]); printf("\n");
	printf( "b[] "); for(u32 ix = 0 ; ix < N; ix++)  printf(" %10.4lf", b[ix]); printf("\n");
}


/*
     r[] is a 2N+1 array containing the  r in the toepliz matrix,
     stored in this order  [r(n-1), r(n-2) ...............r(0).............. r(2-n) r(1-n)]
     r(0) is stored in r[n]
     
     y[]  is a N array containg the right side values
     
     x[]  will, with luck, be the solution to  r * x = y
     
     the function return 0 for succes and various other values, if it runs into numerical problems.
     
*/

int  levinson(u32 N,  double *r, double *y, double *x)
{
	 int res = 0;
    char tf[] = "mulAndPrint f";
    char tb[] = "mulAndPrint b";
    char tx[] = "mulAndPrint x";

	
	 	double *f = new double[N] 	;  //foward 
	 	double *b = new double[N] 	;  //backward
	  for(u32 i = 0; i < N; i++) f[i]=b[i]= 0.0;
	  
		if ( r[N]== 0.0 ) { res = -1 ; goto Exit;};
	 
	  //dump(N, r,y,x,f,b);
	   
		x[0] = y[0]/r[N-1]; //first solution  
		
		if( N > 1)  // at least two 
		{
			double divisor = (r[N-1] *r[N-1])-( r[N]*r[N-2]);
			if (divisor == 0.0) { res = -1 ; goto Exit;};
			x[0] = ((y[0] * r[N-1])-(y[1]*r[N-2]))/divisor;
			x[1] = ((y[1] * r[N-1])-(y[0]*r[N]))/divisor;
			f[0] =  r[N-1]/divisor;
			f[1] =  -r[N]/divisor;
			b[0] = -r[N-2]/divisor;
			b[1] =  r[N-1]/divisor;
	   	printf("divisor %10.4lf\n",divisor);
	  }
	  
	  dump(N, r,y,x,f,b);// to here everything looks OK !
/*
    mulAndPrint(tf, N, r, f);
    mulAndPrint(tb, N, r, b);
    mulAndPrint(tx, N, r, x);
    */
		if( N > 2)  // at least three
		
		for( u32 inx = 2 ;  inx < N; inx ++ ){
			  //dump(N, r,y,x,f,b);
				// calculate error terms
				double enf = 0.0;
			  double enb = 0.0;
				for (u32 t = 0 ; t < inx; t++) enf += r[N-1+inx-t] * f[t] ;
				for (u32 t = 0 ; t < inx; t++) enb += r[N-2-t] * b[t] ;
				//printf("inx %d enf %10.4lf  enb %10.4lf\n",inx,  enf, enb);
					
			  double divisor = 1-(enf*enb);
			  if (divisor == 0.0) { res = -1 ; goto Exit;};
			  for( u32 tf = 0 ; tf < inx; tf++) f[tf] /= divisor ;
			  f[inx] = 0.0;
			  			  
			  for( u32 tb = inx ; tb > 0 ; tb--) b[tb] =  b[tb-1]/ divisor ;
			  b[0] = 0.0;
			  
			  for( u32 tt = 0 ; tt <= inx; tt++) {			  	
			  	double newf = f[tt] - enf*b[tt];
			  	double newb = b[tt] - enb*f[tt];
			  	f[tt] = newf;
			  	b[tt] = newb;
			  }
	/*			dump(N, r,y,x,f,b);

    mulAndPrint(tf, N, r, f);
    mulAndPrint(tb, N, r, b);
			*/	
				for( u32 i = inx; i < N; i++) x[i] = 0;
				
				
		  	double enx = 0.0;
			  
			  /*calculate error term */
			  for( u32 xx = 0 ; xx < inx; xx++){ 
			  		  //printf("delta   %10.4lf ",(r[N-1+inx-xx] * x[xx] )  );
			  	    enx += (r[N-1+inx-xx] * x[xx] ) ;			  
			  		  //printf("r  r[%d]= %10.4lf ",N-1+inx-xx, r[N-1+inx-xx] );
			  			//printf("x  x[%d]= %10.4lf\n",xx, x[xx] );
		   				//printf("inx %d Error term  %10.4lf\n",inx,enx);
				}		
			  /* update solution */ 
			  for( u32 yy = 0; yy < inx ; yy++ ) x[yy] += b[yy] * (y[inx]- enx); 
        x[inx] =  (y[inx]- enx) *b[inx];
        //mulAndPrint(tx, N, r, x);
			  
			 //return 0;//debug
		}
	
	Exit:
	  dump(N, r,y,x,f,b);
    mulAndPrint(tx, N, r, x);
    
		delete  [] f;
		delete  [] b;
		return res;	
}


int main( int argc, char **argv  )
{
	int N = 6 ;
	double r[(2*N)-1] = {5.0, 7.0, 11.0, 23.0, 31.0, 43.0, 53.0, 61.0, 73.0, 83.0, 93.0};
	double y[N]     = {3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
	double x[N]     = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	
	levinson( N, r, y, x);
}