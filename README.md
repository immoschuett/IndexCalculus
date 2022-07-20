# Index Calculus Algorithm
# by Immo Schütt, Berenike Dieterle

A small programm to compute discrete Logarithms in small Fields where p is a Save-Prime.
I.e. If   F_p=<g>, then log_g(a) = l => g^l = a  
 
1. Precomputation: Sieve out smooth relations to recover logarithms of a factor base.
2. Compute individual logarithms (naive randomized)

# Time for Precomputation:
usualle the longer part:
fast for p as Int64: <1 seconds for p < 2^64

TODO list so far:
 
>> Implement good logger to quick_overview/present performance/correctness
 
>> Debug/improve with @code_warntype 
 
>> Optimize with  @profile
 
>> Implement block wiedemann for faster performanve 
 
>> Implement sieve for faster finding l for individual logs
 
>> Inplace sieving operations as good as possible
 
>> Preprocess >>> get A to be squared and sparse. ? but lose relations then we can save some time in wiedemann.
 
>> And much more.
