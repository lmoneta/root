const static char *opencl_src = "\n\
typedef ulong uint64_t;\n\
typedef uint uint32_t;\n\
typedef uchar uint8_t;\n\
 inline float u01_closed_closed_32_24(uint32_t i){\n\
    return i*(1.f/4294967296.f);\n\
}\n\
 inline float u01_closed_open_32_24(uint32_t i){\n\
    return (i>>8)*(1.f/16777216.f);\n\
}\n\
 inline float u01_open_closed_32_24(uint32_t i){\n\
    return (1+(i>>8))*(1.f/16777216.f);\n\
}\n\
 inline float u01_open_open_32_24(uint32_t i){\n\
    return (1+(i>>8))*(16777215.f * (1.f/16777216.f) * (1.f/16777216.f));\n\
}\n\
 inline double u01_closed_closed_64_53(uint64_t i){\n\
    return i*(1./(4294967296.*4294967296.));\n\
}\n\
 inline double u01_closed_open_64_53(uint64_t i){\n\
    return (i>>11)*(1./(4294967296.*2097152.));\n\
}\n\
 inline double u01_open_closed_64_53(uint64_t i){\n\
    return (1+(i>>11))*(1./(4294967296.*2097152.));\n\
}\n\
 inline double u01_open_open_64_53(uint64_t i){\n\
    return (1+(i>>11))*(9007199254740991.*(1./(4294967296.*2097152.))*(1./(4294967296.*2097152.)));\n\
}\n\
 inline double u01_closed_closed_32_53(uint32_t i){\n\
    return i*(4294967297.*(1./4294967296.)*(1./4294967296.));\n\
}\n\
 inline double u01_closed_open_32_53(uint32_t i){\n\
    return i*(1./4294967296.);\n\
}\n\
 inline double u01_open_closed_32_53(uint32_t i){\n\
    return (1.+i)*(1./4294967296.);\n\
}\n\
 inline double u01_open_open_32_53(uint32_t i){\n\
    return (0.5+i)*(1./4294967296.);\n\
}\n\
\n\
typedef struct {\n\
    uint64_t input;\n\
    const char *expected[(4 + 8)];\n\
} KatU01Test;\n\
\n\
typedef struct {\n\
    float valuef[4];\n\
    double valued[8];\n\
} KatU01Result;\n\
__kernel void\n\
dev_execute_tests(__global uint64_t *vals, __global KatU01Result *ret) {\n\
    unsigned tid = get_global_id(0);\n\
    uint64_t v64 = vals[tid];\n\
    uint32_t v32 = vals[tid] & 0xffffffff;\n\
    if (v32 == v64) {\n\
 ret[tid].valuef[0] = u01_closed_closed_32_24(v32);\n\
 ret[tid].valuef[1] = u01_closed_open_32_24(v32);\n\
 ret[tid].valuef[2] = u01_open_closed_32_24(v32);\n\
 ret[tid].valuef[3] = u01_open_open_32_24(v32);\n\
\n\
 ret[tid].valued[0] = u01_closed_closed_32_53(v32);\n\
 ret[tid].valued[1] = u01_closed_open_32_53(v32);\n\
 ret[tid].valued[2] = u01_open_closed_32_53(v32);\n\
 ret[tid].valued[3] = u01_open_open_32_53(v32);\n\
    } else {\n\
 int i;\n\
 for (i = 0; i < 4; i++) {\n\
     ret[tid].valuef[i] = 0.f;\n\
     ret[tid].valued[i] = 0.;\n\
 }\n\
    }\n\
\n\
    ret[tid].valued[4] = u01_closed_closed_64_53(v64);\n\
    ret[tid].valued[5] = u01_closed_open_64_53(v64);\n\
    ret[tid].valued[6] = u01_open_closed_64_53(v64);\n\
    ret[tid].valued[7] = u01_open_open_64_53(v64);\n\
}\n\
";
