#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[6]:


def findallprimes(maxnum):
    
    Bool = np.ones(maxnum, dtype = bool)
    Bool[0] = 0
    maxiter = int(np.sqrt(maxnum))
    for i in range(2, maxiter+1):
        if Bool[i-1]:
            # Sieve of Eratosthenes
            Bool[np.arange(2*i-1, maxnum, i)] = False
    
    return np.where(Bool)[0] + 1

def trial_division(n, primeList):
    
    if n < 2:
        return []
        
    prime_factors = []
    for p in primeList:
        if p*p > n: break
        while n % p == 0:
            prime_factors.append(p)
            n //= p
    if n > 1:
        prime_factors.append(n)
    return prime_factors


# In[10]:


def get_witnesses(n):
    assert (n > 2) and (n % 2 == 1)
    if n < 2047:
        # References: [1], [2], [4]
        witnesses = (2,)
    elif n < 1373653:  # ~1.3 million
        # References: [1], [2], [3], [4]
        witnesses = (2, 3)
    elif n < 9080191:  # ~9.0 million
        # References: [3], [4]
        witnesses = (31, 73)
    elif n < 25326001:  # ~25.3 million
        # References: [1], [2], [3], [4]
        witnesses = (2, 3, 5)
    elif n < 3215031751:  # ~3.2 billion
        # References: [1], [2], [3], [4]
        witnesses = (2, 3, 5, 7)
    elif n < 4759123141:  # ~4.7 billion
        # References: [3], [4]
        witnesses = (2, 7, 61)
    elif n < 2152302898747:  # ~2.1 trillion
        # References: [1], [2], [3], [4]
        witnesses = (2, 3, 5, 7, 11)
    elif n < 3474749660383:  # ~3.4 trillion
        # References: [1], [2], [3], [4]
        witnesses = (2, 3, 5, 7, 11, 13)
    elif n < 341550071728321:  # ~341.5 trillion
        # References: [1], [2], [3], [4]
        witnesses = (2, 3, 5, 7, 11, 13, 17)
    elif n < 3825123056546413051:  # ~3.8 million trillion
        # References: [1], [4]
        witnesses = (2, 3, 5, 7, 11, 13, 17, 19, 23)
    elif n <= 2**64:
        witnesses = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    else:
        witnesses = None
        
    return witnesses

def is_miller_rabin_prime(n):
    """
    Returns True for primes, False for non-primes and raises ValueError for
    unsure cases.
    """
    if n <= 1:
        return False
    elif n == 2:
        return True
    elif n % 2 == 0:
        return False
    
    witnesses = get_witnesses(n)
    if witnesses is None:
        msg = 'No definite Miller-Rabin test is available for %d' % n
        raise ValueError(msg)
    
    d, s = factorN(n-1)
    for w in witnesses:
        if is_composite(w, n, d, s):
            return False
    
    return True

def is_composite(witness, n, d, s):
    
    assert d*2**s == n-1
    
    if pow(witness, d, n) == 1:
        return False
    for i in range(s):
        if pow(witness, 2**i*d, n) == n-1:
            return False
    
    return True
  
def factorN(n):
    
    s = 0
    d = n
    while True:
        q, r = divmod(d, 2)
        if r == 1:
            break
        s += 1
        d = q
    
    assert d*2**s == n
    
    return d, s


# In[11]:


def spiralPrimes(ratio = 0.1, maxfactor = 30000, method = 'Miller-Rabin'):
    
    if method == 'trial-division':
        primeList = findallprimes(maxfactor)
    
    n = 3
    count = 0
    while True:
        assert n < maxfactor, 'The primeList is limited, raise maxfactor'
        
        total = 2 * n - 1
        for num in [n**2-n+1, n**2-2*n+2, n**2-3*n+3]:
            if method == 'trial-division':
                primefactor = trial_division(num, primeList)
                if len(primefactor) == 1:
                    count += 1
            elif method == 'Miller-Rabin':
                if is_miller_rabin_prime(num):
                    count += 1
        
        R = count / float(total)
        if R < ratio:
            return n
        
        n += 2

if __name__ == '__main__':
    print(spiralPrimes()) # ~844ms


# In[ ]:




