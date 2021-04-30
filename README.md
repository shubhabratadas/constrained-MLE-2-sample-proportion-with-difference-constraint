# constrained-MLE-2-sample-proportion-with-difference-constraint
#To find constrained MLE in 2 sample proportion with difference given
# Index populations 1 & 2 such that delta >= 0

# n1 = sample size from the first population -- it should be a positive integer
# n2 = sample size from the second population -- it should be a positive integer

# x1 = number of successes in the sample drawn from the first population -- it should be a nonnegative integer less than or equal to n1
# x2 = number of successes in the sample drawn from the second population -- it should be a nonnegative integer less than or equal to n2

# delta = pi1 - pi2 = difference in two population proportions.
# Index the populations such that delta is nonnegative 
# delta should be between 0 and 1 

# p2_mle_delta function returns the MLE of pi2 based on (x1,n1,x2,n2,delta) when pi1 - pi2 = delta
# Examples
p2_mle_delta(10, 10, 0,	30, 1)
p2_mle_delta(45, 125, 35, 75, 0)
p2_mle_delta(40, 40, 0, 10, 0.2)
p2_mle_delta(40, 40, 0,	10, 0.3)
p2_mle_delta(100, 100,	0, 200, 0.55)
p2_mle_delta(10, 50, 0, 10, 0.2)
p2_mle_delta(10, 10, 19, 20, 0.1)
p2_mle_delta(15, 50, 0,	10, 0.2)
p2_mle_delta(15, 50, 15, 400, 0.2)

# To find the standard error of SE(p1 - p2) where p1, p2 are the constrained MLE of Pi1 and Pi2 given Pi1 - Pi2 = delta 
# use the SE_del function
# Example
SE_del(15, 50, 15, 400, 0.2)
