uint64_t RDTSC_start_clk;
uint64_t RDTSC_end_clk;
uint64_t RDTSC_total_clk;
uint64_t etime;

#define frequency 3.1
#define GHz       1000000000

#define STAMP0 ({\
    uint64_t start; \
    unsigned cycles_low, cycles_high; \
    asm volatile ("CPUID\n\t" \
                  "RDTSCP\n\t" \
                  "mov %%edx, %0\n\t" \
                  "mov %%eax, %1\n\t": "=r" (cycles_high), "=r" (cycles_low) \
                  : \
                  : \
                  "%rax", "%rbx", "%rcx", "%rdx"); \
    start = ( ((uint64_t)cycles_high << 32) | cycles_low ); \
    start;})

#define STAMP1 ({\
    uint64_t end; \
    unsigned cycles_low1, cycles_high1; \
    asm volatile ("RDTSCP\n\t" \
                  "mov %%edx, %0\n\t" \
                  "mov %%eax, %1\n\t" \
                  "CPUID\n\t": "=r" (cycles_high1), "=r" (cycles_low1) \
                  : \
                  : \
                  "%rax", "%rbx", "%rcx", "%rdx"); \
    end = ( ((uint64_t)cycles_high1 << 32) | cycles_low1 ); \
    end;})

#define MEASURE(x) \
	RDTSC_start_clk = STAMP0; \
    {x}; \
	RDTSC_end_clk = STAMP1; \
	RDTSC_total_clk = (RDTSC_end_clk - RDTSC_start_clk); \
    etime = RDTSC_total_clk;
