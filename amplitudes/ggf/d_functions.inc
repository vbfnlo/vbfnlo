*  Include file for the cached Passarino-Veltman 4-point-functions, i.e. D0, Dij
*  Content: Common block with the cache-table
*  Author: Gunnar Klaemke
*  Date: 16.10.07
*  Modified: Michael Kubocz
*  Date: 19.09.07
*  Cache_Size and Cache_values were increased due to new DD-coefficients

c     parameters that define the size of the arrays		
      INTEGER Cache_Size
      PARAMETER (Cache_Size=10000) 
      INTEGER LT_Size
      PARAMETER (LT_Size=10000) ! Size of look-up table - can't be big enough	
c     Cache_index: the current position for a new cash entry	
c     Cache_count: the number of entries in the cache (should be smaller than Cache_Size)
      INTEGER Cache_index,Cache_count
c     the cached function-results	
      COMPLEX*16 Cache_values(1:31,1:Cache_Size)
c     the keys for the cached values
c     use equivalence to strings for faster comparisons
      REAL*8 Cache_keys(7,1:Cache_Size)	
      CHARACTER*(8*7) Cache_keys_Data(1:Cache_Size)
      EQUIVALENCE (Cache_keys,Cache_keys_Data)
c     table and call-counter for fast look-ups
      INTEGER LookupTable(1:LT_Size)
      INTEGER CallCounter 
c     everything in one common block
      COMMON/D_1M_CACHE/Cache_values, Cache_keys, LookupTable,
     &                  Cache_index, Cache_count, CallCounter

