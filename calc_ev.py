import itertools
import numpy as np
from scipy.special import factorial


#%%
# [2,1,3,4] means h_1 < h_2 < h_4 < h_4

h_values = [1, 2, 3, 4, 5]  
n = len(h_values) +1 
k = 1 # originallyk is k+1 here

# All Permutations
permutations = list(itertools.permutations(h_values))

def is_valid_permutation(h, k):
    # L_{0,n}
    
    if(k==1):
        if h[k] > max(h[:k]):
            return False    
    # L_{i,n}^k für i von 0 bis n-k-1
        i = 0
        for i in range(n - k - 2):
            if min(h[i], h[i + k+1]) > h[i+1]:
                return False
    
    # L_{n-k,n}
        if h[n - k-2] > h[n-k-1]:
            return False
        return True
    else:
        i = 0
        
        #print(h[n - k-2],(h[n - k-1:]))
        #print(h[k], h[:k])
        if h[k] > max(h[:k]): # überprüft!
            return False
    
    # L_{i,n}^k für i von 0 bis n-k-1
        for i in range(n - k - 2): # überprüft!
            if min(h[i], h[i + k]) > max(h[i+1:i+1+k]):
                return False
    
    # L_{n-k,n}
        if h[n - k-2] > max(h[n - k-1:]): 
            return False 
        return True

#print(is_valid_permutation([4, 1, 2, 3], k))
valid_permutations = [perm for perm in permutations if is_valid_permutation(perm, k)]

print(f"Number of Valid Permutations: {len(valid_permutations)}")
print("Valid Permutations:")
for perm in valid_permutations:
    print(perm)
    
print(len(valid_permutations)/factorial(n-1))
#%%
h_values = [1, 2, 3, 4, 5, 6, 7, 8]  
n = len(h_values) +1 
#k = 1 # originallyk is k+1 here

# All Permutations
permutations = list(itertools.permutations(h_values))

k = [1,2,3,4, 5, 6, 7]
for i in k:
    valid_permutations = [perm for perm in permutations if is_valid_permutation(perm, i)]
    print(len(valid_permutations)) #/factorial(n-1))
print(factorial(n-1))
#%%
#n = 10
161280
204120
258048
241920
259200
272160
282240
362880.0

# n = 11
1612800
2041200
2322432
2419200
2592000
2721600
2822400
2903040
3628800.0
