#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include <iostream>
#include <math.h>

extern "C"{
    double *DRang_IW_C(
            double *xyz1, 
            double *xyz2, 
            double *w1, 
            double *w2, 
            long dim1, 
            long dim2, 
            long dim3, 
            double *lsbin, long nsbin, 
            double *labin, long nabin)
    {
    // -- xyz1, data # (Ndata, 3)
    // -- xyz2, rand # (Nrand, 3)
    // -- dim1, Ndata 
    // -- dim2, Nrand  
    // -- dim3, Nw, type number of individual weights
    // -- w1, weight of data
    // -- w2, weight of rand
        
    double PI = 3.141592653;
    double sep3d_lo, sep3d_up;
    double sep2d_lo, sep2d_up; 
    long ii, jj; 
        
    // -- initialize 
    const long dim = (nsbin+1)*(nabin+1)*dim2; 
    double *DDang; 
    DDang = new double[dim]; 
    for(ii=0; ii<dim; ii++){DDang[ii]=0.0;}

    sep3d_lo =  pow(10.0, lsbin[0]);
    sep3d_up =  pow(10.0, lsbin[nsbin-1]);
    
    sep2d_lo =  pow(10.0, labin[0]);
    sep2d_up =  pow(10.0, labin[nabin-1]);
        
// -- loop 
#pragma omp parallel
{
    double dx, dy, dz, r1, r2; 
    double sep2d, sep3d; 
    long i, j, k, isbin, iabin;
    
    double *DDang_private;  
    DDang_private = new double[dim];
    for(i=0; i<dim; i++){DDang_private[i]=0;}

    #pragma omp for
    for (i=0;i<dim1;i++){
        for (j=0;j<dim2;j++){
            
            r1 = xyz1[i+dim1*0]*xyz1[i+dim1*0] + xyz1[i+dim1*1]*xyz1[i+dim1*2] + xyz1[i+dim1*0]*xyz1[i+dim1*2];
            r2 = xyz2[j+dim2*0]*xyz2[j+dim2*0] + xyz2[j+dim2*1]*xyz2[j+dim2*2] + xyz2[j+dim2*0]*xyz2[j+dim2*2]; 
            r1 = sqrt(r1);
            r2 = sqrt(r2);
            
            dx    = xyz1[i+dim1*0] - xyz2[j+dim2*0];
            dy    = xyz1[i+dim1*1] - xyz2[j+dim2*1];
            dz    = xyz1[i+dim1*2] - xyz2[j+dim2*2];
            sep3d = sqrt(dx*dx + dy*dy + dz*dz);

            dx    = xyz1[i+dim1*0] - xyz2[j+dim2*0];
            dy    = xyz1[i+dim1*1] - xyz2[j+dim2*1];
            dz    = xyz1[i+dim1*2] - xyz2[j+dim2*2];
            sep2d = sqrt(dx*dx+dy*dy+dz*dz); 
            sep2d = 2.0*asin(0.5*sep2d)*180/PI; 

         if ( sep3d <  sep3d_lo )
            {isbin = 0;}
         else if (sep3d >= sep3d_up)
            {isbin = nsbin;}
         else
            {
            sep3d  = log10(sep3d);
            isbin  = (sep3d - lsbin[0])/(lsbin[1] - lsbin[0])  + 1;
            }

         if ( sep2d <  sep2d_lo)
            {iabin = 0;}
         else if ( sep2d >= sep2d_up )
            {iabin = nabin;}
         else
            { 
            sep2d = log10(sep2d);
            iabin = (sep2d - labin[0])/(labin[1] - labin[0]) + 1;
            } 
            
            for (k=0;k<dim3;k++){

         // if out of range, no count. 
 
         // printf("rank %2d %2d; %2d %2d %2d %2d %lf %lf \n", rank, k, i, j, isbin, iabin, sep3d, sep2d);
        
         DRang_private[ (isbin*(nabin+1) + iabin ] += 1; 
         // WWang_private[]

            } // end for k, looping weights w1, w2, w3,...
        } // end for j, looping rand
    } // end for i, looping data 

    #pragma omp critical
    { 
    for(i=0; i<dim; i++){
        DDang[i] += 2*DDang_private[i];} // end for for 
    }

} // end pragma omp parallel 

// for(jj=0; jj<nsbin+1; jj++){
//    for(ii=0; ii<nabin+1; ii++){ 
//        printf("%lf ", DDang[ jj*(nabin+1) + ii]); 
//    }
//    printf("\n"); 
// }

return DDang;
}
} 
