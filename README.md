
Model the fragment size distribution of B73 (Maize) gDNA based on fragment analyzer result and calculate the picomole amount of DNA for that distribution of DNA fragment sizes.
===========================================================================================================================================================================================


FA based gDNA size distribution for B73(Maize) gDNA:
================================================
![alt tag](https://github.com/ffrancis/Model_gDNAsizedistAND3prime_pmol_calc/blob/master/FA_distributionMaize_gDNA.png)


Modelled distribution of gDNA size distribution for B73(Maize) gDNA:
================================================
![alt tag](https://github.com/ffrancis/Model_gDNAsizedistAND3prime_pmol_calc/blob/master/Size_distribution_gDNA121614.png)

Description     :
                        
                        The fragment size of gDNA is a continuous variable, and its probability distribution is a continous probability distribution.
                        (The triangular distribution is often used in ill-defined problems where the underlying distribution is not known, but some knowledge of the limits and mode exists)
                        Also refer to http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.triangular.html for additional information.
                                                                       
                        pmol of 3' ends of X ug of a N bp dsDNA fragment        = 2 * 10^6 * X ug (of dsDNA)) / MW (in Da)

                                                                                = 2 * 10^6 * X ug (of dsDNA)) / N bp * 660 (Da)
                        400 U TdT can be used with approx. 100 pmol 3' ends
                        In NEB TdT, we use 6ul * 20U = 120 U TdT
                        For 120 U we can use: (100/400)* 120 = 30 pmol 3'ends
                        So, the amount of dS gDNA that can be used (with 120U TdT), given that it has Nbp:
                        X ug (of dsDNA) = 30 * N bp * 660 (Da) / 2 * 10^6 

                
Note            :       The probability density function for the triangular distribution can written as a function below.


                        def triang_dist(x, l, m, r):
                                if x >= l and x <= m:
                                    p = (2*(x-l))/((r-l)*(m-l))
                                if x >= m and x <= r:
                                    p = (2*(m-x))/((r-l)*(r-m))
                                if not (x >= l and x <= m) or (x >= m and x <= r):
                                    p = 0
                                #return p
                                print p
                                
                        Where:  l = float(min_value)
                                r = float(max_value)
                                m = float(mode)
                                x = random number for which you need to find the likelihood to be in the triangular distribution.

