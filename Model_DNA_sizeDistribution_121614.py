#!/usr/bin/env python

'''
Description     :       Model the fragment size distribution of B73 (Maize) gDNA based on fragment analyzer result and calculate the picomole amount of DNA for that particular distribution of DNA fragment size.
                        
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
'''

import numpy
import matplotlib.pyplot as plt


min_value       = 1200          # Lower bound of size distribution from fragment analyser
max_value       = 40000         # Upper bound of size distribution from fragment analyser
mode            = 29050         # Peak from fragment analyser
num_data_points = 10000


# Model the distribution of the gDNA fragments based on the available parameters(min_value, mode and max_value)
distribution = numpy.random.triangular(min_value, mode, max_value, size=num_data_points)



###################################################################
# Functions to calculate the amount of DNA (in ug) for a given fragment size, that can be used with 120 U TdT
###################################################################

# Functions to calculate the amount of DNA (in ug) for a given fragment size, that can be used with 120 U TdT (ROCHE PROTOCOL)
'''
pmol of 3' ends of X ug of a N bp dsDNA fragment        = 2 * 10^6 * X ug (of dsDNA)) / MW (in Da)

                                                        = 2 * 10^6 * X ug (of dsDNA)) / N bp * 660 (Da)
400 U TdT can be used with approx. 100 pmol 3' ends
In NEB, we use 6ul * 20U = 120 U TdT
For 120 U = (100/400)* 120 = 30 pmol 3'ends
So, X ug (of dsDNA) = 30 * N bp * 660 (Da) / 2 * 10^6 
'''
#def dna_amount_cal(frag_len):
#        est_dna_amount = (30 * frag_len * 660)/2000000
#        return est_dna_amount



# Functions to calculate the amount of DNA (in ug) for a given fragment size, that can be used with 120 U TdT (ROCHE PROTOCOL Felix modified based on B73 experimental results)
'''
Based on the size distribution and the experiment conducted in lab, the TdT effect (120U) seems to saturate at <= 3ug gDNA (safe amount 2.5ug)
So we have estimated that for 120 U TdT, 0.325 pmol 3' ends of B73 genomic DNA can be used.

So,
X ug (of dsDNA) = 0.325 * N bp * 660 (Da) / 2 * 10^6
'''
def dna_amount_cal(frag_len):
        est_dna_amount = (0.325 * frag_len * 660)/2000000
        return est_dna_amount
###################################################################



###################################################################
# Calculate the average amount of DNA that can be used with 120 U TdT
###################################################################
dna_amount = []
for i in distribution:
        dna_amount.append(dna_amount_cal(i))

average_dna_amount = sum(dna_amount ) / float(num_data_points)
print average_dna_amount
###################################################################


###################################################################
##Plot the distribution
###################################################################
fig = plt.figure()
fig.suptitle('Size distribution of B73 gDNA fragments', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.87)


text1 = "min_value = "+str(min_value)
text2 = "; max_value = "+str(max_value)
text3 = "; mode = "+str(mode)
combined_text = text1+text2+text3


ax.set_title(combined_text) # Title with info regarding the min_value, max_value and mode
ax.set_xlabel('DNA fragment length')
ax.set_ylabel('Count')

L = max_value - min_value
bin = L/80
ax.hist(distribution,color='grey',alpha=0.5,bins=bin)
plt.show()
###################################################################


