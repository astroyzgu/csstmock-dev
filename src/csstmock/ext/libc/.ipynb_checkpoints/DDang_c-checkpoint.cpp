#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include <iostream>
#include <math.h>

extern "C"{
    double *DDang_C(double *xyz1, double *xyz2, 
            long dim1, long dim2, 
            double *lsbin, long nsbin, 
            double *labin, long nabin) 
    {
    double PI = 3.141592653;
    double sep3d_lo, sep3d_up;
    double sep2d_lo, sep2d_up; 
    long ii, jj; 
    // -- initialize 
    const long dim = (nsbin+1)*(nabin+1); 
    double *DDang; 
    DDang = new double[dim]; 
    for(ii=0; ii<dim; ii++){DDang[ii]=0.0;}

    sep3d_lo =  pow(10.0, lsbin[0]);
    sep3d_up =  pow(10.0, lsbin[nsbin-1]);
    // printf("sep3d: %lf %lf %lf \n", sep3d_lo, pow(10.0, lsbin[1]), sep3d_up);

    sep2d_lo =  pow(10.0, labin[0]);
    sep2d_up =  pow(10.0, labin[nabin-1]);
    // printf("sep2d: %lf %lf %lf \n", sep2d_lo, pow(10.0, labin[1]), sep2d_up);
    // printf("%ld %ld \n", dim1, dim2); 
    for (ii=0;ii<1;ii++){
        // printf("First line in C: %12.8lf %12.8lf %12.8lf \n", xyz1[ii+dim1*0], xyz1[ii+dim1*1], xyz1[ii+dim1*2] );
    } 
    // printf("\n"); 



        
// -- loop 
#pragma omp parallel
{
    double dx, dy, dz; 
    double sep2d, sep3d; 
    long i, j, k, isbin, iabin;
    // const int size = omp_get_num_threads(); 
    const int rank = omp_get_thread_num(); 
    double *DDang_private;  
    DDang_private = new double[dim];
    for(i=0; i<dim; i++){DDang_private[i]=0;}

    #pragma omp for
    // for (k=0;k<dim1*dim1;k++){
    // i =   k/dim1; 
    //   j = k-i*dim1;
    for (i=0;i<dim1;i++){
        for (j=0;j<dim1;j++){
         if(i<=j){continue;}
         k = i*dim1 + j;
         
         // calculate s of each pair
         dx    = xyz1[i+dim1*0] - xyz1[j+dim1*0];
         dy    = xyz1[i+dim1*1] - xyz1[j+dim1*1];
         dz    = xyz1[i+dim1*2] - xyz1[j+dim1*2];
         sep3d = dx*dx + dy*dy + dz*dz;
            
         // if (sep3d > 0.01){
         //    printf("%10d ; %ld %ld %12.8lf %12.8lf \n", k, i, j, sep3d, sep3d );
         //    printf("%10d ; %12.8lf %12.8lf %12.8lf \n", i, xyz1[i*dim2+0], xyz1[i*dim2+1], xyz1[i*dim2+2] );
         //    printf("%10d ; %12.8lf %12.8lf %12.8lf \n", j, xyz1[j*dim2+0], xyz1[j*dim2+1], xyz1[j*dim2+2] );
         //    printf("%10d ; %12.8lf %12.8lf %12.8lf = %12.8lf  %12.8lf \n", k, dx, dy, dz, sep3d, sep3d );
         // }

         // calculate a of the pair 
         dx    = xyz2[i+dim1*0] - xyz2[j+dim1*0];
         dy    = xyz2[i+dim1*1] - xyz2[j+dim1*1];
         dz    = xyz2[i+dim1*2] - xyz2[j+dim1*2];
         sep2d = sqrt(dx*dx+dy*dy+dz*dz); 
         sep2d = 2.0*asin(0.5*sep2d)*180/PI; 
         sep3d = sqrt(sep3d); 
        
         // if out of range, no count. 
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
         // printf("rank %2d %2d; %2d %2d %2d %2d %lf %lf \n", rank, k, i, j, isbin, iabin, sep3d, sep2d);
        
         DDang_private[isbin*(nabin+1) + iabin] += 1;
         // for(int i=0; i<dim; i++){ printf("%2.0lf ", DDang_private[i]);} // end for for 
         // printf("\n"); 
        }
    } // end for k 

    #pragma omp critical
    { 
    for(i=0; i<dim; i++){
        DDang[i] += 2*DDang_private[i];} // end for for 
    }

} // end pragma omp parallel 

//for(jj=0; jj<nsbin+1; jj++){
//    for(ii=0; ii<nabin+1; ii++){ 
//        printf("%lf ", DDang[ jj*(nabin+1) + ii]); 
//    }
//    printf("\n"); 
//}

return DDang;
}
} 