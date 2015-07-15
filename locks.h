#if defined(_OPENMP)

/**
 * Locks rows
 */
void lock(int i, bool *locks){
    
    // implementation of a spin lock
    for (bool locked = false; locked == false;){
        
        #pragma omp critical (LockRegion)
        {
            locked = !locks[i-1] && !locks[i] && !locks[i+1];
            if (locked){
                locks[i-1] = true; locks[i] = true; locks[i+1] = true;
            }
        }
    
    }
    
}

/**
 * Unlocks rows
 */
void unlock(int i, bool *locks){
    #pragma omp critical (LockRegion)
    {
        locks[i-1] = false; locks[i] = false; locks[i+1] = false;
    }
}

#endif