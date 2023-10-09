#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include <iostream>
#include <math.h>

extern "C"{
    double *DRang_IW_c(
            double *xyz1, 
            double *xyz2, 
            double *wht1, 
            double *wht2, 
            long dim1, 
            long dim2, 
            long dim3, 
            double *lsbin, long nsbin, 
            double *labin, long nabin, 
            long isDD)
    {
    // -- xyz1, data # (Ndata, 3)
    // -- xyz2, rand # (Nrand, 3)
    // -- w1, weight of data
    // -- w2, weight of rand
    // -- dim1, Ndata 
    // -- dim2, Nrand  
    // -- dim3, Nw, type number of individual weights
        
    double PI = 3.141592653;
    double sep3d_lo, sep3d_up;
    double sep2d_lo, sep2d_up; 
    long ii, jj, kk; 
        
    // -- initialize 
    const long dim = (nsbin+1)*(nabin+1)*dim3; 
    double *DRang; 
    DRang = new double[dim]; 
    for(ii=0; ii<dim; ii++){DRang[ii]=0.0;}

    sep3d_lo =  pow(10.0, lsbin[0]);
    sep3d_up =  pow(10.0, lsbin[nsbin-1]);
    
    sep2d_lo =  pow(10.0, labin[0]);
    sep2d_up =  pow(10.0, labin[nabin-1]);
        
// -- loop 
#pragma omp parallel
{
    double dx, dy, dz; 
    double x1, y1, z1, r1, w1; 
    double x2, y2, z2, r2, w2; 
    double sep2d, sep3d; 
    long i, j, jstart, k, isbin, iabin;
    
    double *DRang_private;  
    DRang_private = new double[dim];
    for(i=0; i<dim; i++){DRang_private[i]=0;}
    
    
    #pragma omp for
    for (i=0;i<dim1;i++){
        
        x1 = xyz1[3*i+0];
        y1 = xyz1[3*i+1];
        z1 = xyz1[3*i+2];
        r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
        // w1 = wht1[dim3*i+0];
        
        if(isDD!=0)
            {jstart =i+1;}
        else 
            {jstart = 0; }
        
        
        for (j=jstart;j<dim2;j++){
                        
            x2 = xyz2[3*j+0];
            y2 = xyz2[3*j+1];
            z2 = xyz2[3*j+2];
            r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
            
            //  if ((i==1)&&(j==1))
            //    {printf("i==%ld,  %12.8lf  %12.8lf  %12.8lf %12.8lf \n", i, x1, y1, z1, w1);} 
            
            // if ((i==1)&&(j==1))
            //    {printf("j==%ld,  %12.8lf  %12.8lf  %12.8lf %12.8lf \n", j, x2, y2, z2, w2);} 
            
            dx    = x1 - x2;
            dy    = y1 - y2;
            dz    = z1 - z2;
            sep3d = sqrt(dx*dx + dy*dy + dz*dz);


         if ( (sep3d > 0)&&(sep3d < sep3d_lo) )
            {isbin = 0;}
         else if (sep3d >= sep3d_up)
            {isbin = nsbin;}
         else if (sep3d <= 0)
            {continue;} 
         else
            {
            isbin  = ( log10(sep3d) - lsbin[0])/(lsbin[1] - lsbin[0])  + 1;
            }

            dx    = x1/r1 - x2/r2;
            dy    = y1/r1 - y2/r2;
            dz    = z1/r1 - z2/r2;
            sep2d = sqrt(dx*dx + dy*dy + dz*dz);
            sep2d = 2.0*asin(0.5*sep2d)*180/PI; 
            
         if ( sep2d < sep2d_lo)
            {iabin = 0;}
         else if ( sep2d >= sep2d_up )
            {iabin = nabin;}
         else
            { 
            iabin = (log10(sep2d) - labin[0])/(labin[1] - labin[0]) + 1;
            } 
            
            // printf("%ld, %ld, %lf - %ld, %lf - %ld \n",i,j, log10(sep3d), isbin, log10(sep2d), iabin);
            
         for (k=0;k<dim3;k++){
             w1 = wht1[dim3*i+k];
             w2 = wht2[dim3*j+k];
             // (isbin, iabin, iw)
             DRang_private[ isbin*(nabin+1)*dim3 + iabin*dim3 + k ] += w1*w2; 

            } // end for k, looping weights w1, w2, w3,...
        } // end for j, looping rand
    } // end for i, looping data 

    #pragma omp critical
    { 
    for(i=0; i<dim; i++){DRang[i] += DRang_private[i];} 
    }

} // end pragma omp parallel 

if (isDD != 0)
{
for(ii=0; ii<dim; ii++){DRang[ii] = 2*DRang[ii];}     
}
// for (kk=0; kk<dim3; kk++){
// printf("kk == %ld\n", kk); 
//  for(jj=0; jj<nsbin+1; jj++){
//     for(ii=0; ii<nabin+1; ii++){ 
//         printf("%lf ", DRang[ jj*(nabin+1)*dim3 + ii*dim3 + kk]); 
//     }
//     printf("\n"); 
//  }
// printf("\n"); 
// }

return DRang;
}
} 
