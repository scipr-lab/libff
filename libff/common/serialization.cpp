namespace libff {

#ifdef BINARY_OUTPUT
bool binary_output = true;
#else
bool binary_output = false;
#endif
 
#ifdef MONTGOMERY_OUTPUT
bool montgomery_output = true;
#else
bool montgomery_output = false;
#endif

#ifdef NO_PT_COMPRESSION
bool no_pt_compression = true;
#else
bool no_pt_compression = false;
#endif

}