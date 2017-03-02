1. Add hash-counter serialization / deserialization for stream
2. Add OMP tasking to stream
3. Reduce mem usage of reference generation for stream
4. Use hash tables for lookup of hashes, rather than arrays, to make lookup linear time in read length rather than (read len + ref len)
5. Check for page faults and mem leaks.
