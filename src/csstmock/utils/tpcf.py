import numpy as np 
from Corrfunc.mocks import DDrppi_mocks, DDsmu_mocks 


class tpcf_corrfunc():  
        """
        tools to calculate the two-point correlation function (Build-in Corrfunc). 
        """
        def smu(data, rand, bins, mumax = 1, nmu = 1, nthreads = 12): 
                """
                This function calculates the correlation function ξ(s,μ) for a given set of data and random points.
                
                Parameters:
                data (array): The data points, ra, dec, comoving-distance, weight.
                rand (array): The random points, ra, dec, comoving-distance, weight.
                bins (array): The bin edges for the correlation function.
                mumax (float): The maximum value of μ. Default is 1.
                nmu (int): The number of μ bins. Default is 1.
                nthreads (int): The number of threads to use. Default is 12.
                
                Returns:
                corr (array): The calculated correlation function ξ(s,μ).
                """  
                DD_counts = DDsmu_mocks(1, 1, nthreads, mumax, nmu, bins,
                                        data[:,0], data[:,1], data[:,2], weights1 = data[:,3], 
                                        is_comoving_dist=True, weight_type = "pair_product")
                DR_counts = DDsmu_mocks(0, 1, nthreads, mumax, nmu, bins,
                                        RA1=data[:,0], DEC1=data[:,1], CZ1=data[:,2], weights1=data[:,3],
                                        RA2=rand[:,0], DEC2=rand[:,1], CZ2=rand[:,2], weights2=rand[:,3],
                                        is_comoving_dist=True, weight_type = "pair_product")
                RR_counts = DDsmu_mocks(1, 1, nthreads,  mumax, nmu, bins,
                                        rand[:,0], rand[:,1], rand[:,2], weights1 = rand[:,3], 
                                        is_comoving_dist=True, weight_type = "pair_product") 

                nsbins  = len(bins) - 1; 
                nmubins = int( len(DD_counts)/nsbins )
                DD = np.array( [ DD[4]*DD[5] for DD in DD_counts] ); total_wnpairs_dd = 1.0*np.sum(DD) 
                DR = np.array( [ DR[4]*DR[5] for DR in DR_counts] )
                RR = np.array( [ RR[4]*RR[5] for RR in RR_counts] ); total_wnpairs_rr = 1.0*np.sum(RR) 
                Ndata   = np.sum(data[:,3]); # np.shape(data)[0]; 
                Nrand   = np.sum(rand[:,3]); # np.shape(rand)[0]; 
                fDD = 1.0/(Ndata*Ndata-np.sum(data[:,3]**2) ); 
                fDR = 1.0/(Ndata*Nrand)
                fRR = 1.0/(Nrand*Nrand-np.sum(rand[:,3]**2) );
                DD = fDD*DD;  
                DR = fDR*DR;
                RR = fRR*RR;       
                corr = 0.0*RR; nonzero = RR != 0.0
                corr[nonzero] = 1.0*(DD[nonzero]-2*DR[nonzero])/RR[nonzero] + 1 
                if nmubins != 1: corr = corr.reshape(nsbins, nmubins)
                return corr
        def wp(data, rand, bins, pimax = 40.0, nthreads = 12):
                """
                This function calculates the projected correlation function wp(rp) for a given set of data and random points.
                
                Parameters:
                data (array): The data points, ra, dec, comoving-distance, weight.
                rand (array): The random points, ra, dec, comoving-distance, weight.
                bins (array): The bin edges for the correlation function.
                pimax (float): The maximum value of pi. Default is 40.0, dpi = 1. 
                nthreads (int): The number of threads to use. Default is 12.
                
                Returns:
                corr (array): The calculated projected correlation function wp(rp).
                """
                nrpbins = len(bins) - 1; 
                npibin  = int( int(pimax)/nrpbins )
                corr = rppi(data, rand, bins, pimax = 40.0, nthreads = 12) 
                corr = corr.reshape(nrpbins, npibin)
                corr = 2*np.sum( corr, axis = 1 )
                return corr
                       
        def rppi(data, rand, bins, pimax = 40.0, nthreads = 12):
                """
                This function calculates the projected correlation function χ(rp, pi) for a given set of data and random points.
                
                Parameters:
                data (array): The data points, ra, dec, comoving-distance, weight.
                rand (array): The random points, ra, dec, comoving-distance, weight.
                bins (array): The bin edges for the correlation function.
                pimax (float): The maximum value of pi. Default is 40.0, dpi = 1. 
                nthreads (int): The number of threads to use. Default is 12.
                
                Returns:
                corr (array): The calculated projected correlation function χ(rp, pi).
                """
                DD_counts = DDrppi_mocks(1, 1, nthreads, pimax, bins, 
                                        data[:,0], data[:,1], data[:,2], weights1 = data[:,3],
                                        is_comoving_dist=True, weight_type = "pair_product") 
                DR_counts = DDrppi_mocks(0, 1, nthreads, pimax, bins,
                                        RA1=data[:,0], DEC1=data[:,1], CZ1=data[:,2], weights1 = data[:,3],
                                        RA2=rand[:,0], DEC2=rand[:,1], CZ2=rand[:,2], weights2 = rand[:,3],
                                        is_comoving_dist=True, weight_type = "pair_product")
                RR_counts = DDrppi_mocks(1, 1, nthreads, pimax, bins, 
                                        rand[:,0], rand[:,1], rand[:,2], weights1 = rand[:,3],
                                        is_comoving_dist=True, weight_type = "pair_product") 
                nrpbins = len(bins) - 1; 
                npibin  = int( len(DD_counts)/nrpbins )
                DD = np.array([DD[4]*DD[5] for DD in DD_counts]); 
                DR = np.array([DR[4]*DR[5] for DR in DR_counts]); 
                RR = np.array([RR[4]*RR[5] for RR in RR_counts]); 
                Ndata   = np.sum(data[:,3]); # np.shape(data)[0]; 
                Nrand   = np.sum(rand[:,3]); # np.shape(rand)[0]; 
                fDD = 1.0/(Ndata*Ndata-np.sum(data[:,3]**2) ); 
                fDR = 1.0/(Ndata*Nrand);
                fRR = 1.0/(Nrand*Nrand-np.sum(rand[:,3]**2) );
                DD = fDD*DD;  
                DR = fDR*DR;
                RR = fRR*RR;
                corr = 0.0*RR; nonzero = RR != 0.0
                corr[nonzero] = 1.0*(DD[nonzero]-2*DR[nonzero])/RR[nonzero] + 1
                #corr = corr.reshape(nrpbins, npibin)
                #corr = 2*np.sum( corr, axis = 1 )
                return corr


class tpcf(tpcf_corrfunc): 
    '''
        tools to calculate the two-point correlation function. 
    '''
    def __init__(self, *arg, **kwargs): 
        super(tpcf, self).__init__(*arg, **kwargs)
