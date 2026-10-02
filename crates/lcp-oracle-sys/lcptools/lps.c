#include "lps.h"


/**
 * @brief Parses a sequence to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores.
 *
 * This function iterates over a sequence defined by iterators `begin` and `end` and identifies key
 * segments, called "cores," that represent the (LCP) regions. By analyzing
 * character relationships in the sequence (such as equality or relative order), it builds and stores
 * these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the sequence to parse.
 * @param end Iterator pointing to the end of the sequence to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @param offset The distance measure where the indecies of the core will be shifted by.
 * @return Size of the cores identified in the given string.
 */
int parse1(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = begin;
    const char *it2 = end;
    int core_index = 0;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (alphabet[(unsigned char)*it1] == alphabet[(unsigned char)*(it1+1)]) {
            continue;
        }

        // check for RINT core
        if (alphabet[(unsigned char)*(it1+1)] == alphabet[(unsigned char)*(it1+2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 + 2;
            while (temp < end && alphabet[(unsigned char)*(temp-1)] == alphabet[(unsigned char)*temp]) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1) {
                    init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                    core_index++;
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                init_core1(&(cores[core_index]), it1, it2-it1, it1-begin+offset, it2-begin+offset);
                core_index++;

                continue;
            }
        }

        if (alphabet[(unsigned char)*it1] > alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] < alphabet[(unsigned char)*(it1+2)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                core_index++;
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, 3, it1-begin+offset, it2-begin+offset);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            alphabet[(unsigned char)*it1] < alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] > alphabet[(unsigned char)*(it1+2)] &&
            alphabet[(unsigned char)*(it1-1)] <= alphabet[(unsigned char)*(it1)] &&
            alphabet[(unsigned char)*(it1+2)] >= alphabet[(unsigned char)*(it1+3)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                core_index++;
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, 3, it1-begin+offset, it2-begin+offset);
            core_index++;

            continue;
        }
    }

    return core_index;
}

/**
 * @brief Parses a sequence to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores using complement alphabet.
 *
 * This function iterates over a sequence defined by iterators `begin` and `end` and identifies key
 * segments, called "cores," that represent the (LCP) regions. By analyzing
 * character relationships in the sequence (such as equality or relative order based on complement), 
 * it builds and stores these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the sequence to parse.
 * @param end Iterator pointing to the end of the sequence to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @param offset The distance measure where the indecies of the core will be shifted by.
 * @return Size of the cores identified in the given string.
 */
int parse2(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = end - 1;
    const char *it2 = begin - 1;
    int core_index = 0;

    // find lcp cores
    for (; begin <= it1 - 2; it1--) {

        // skip invalid character
        if (rc_alphabet[(unsigned char)*it1] == rc_alphabet[(unsigned char)*(it1-1)]) {
            continue;
        }

        // check for RINT core
        if (rc_alphabet[(unsigned char)*(it1-1)] == rc_alphabet[(unsigned char)*(it1-2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 - 2;
            while (begin <= temp && rc_alphabet[(unsigned char)*(temp+1)] == rc_alphabet[(unsigned char)*temp]) {
                temp--;
                middle_count++;
            }
            if (begin <= temp) {
                // check if there is any SSEQ cores left behind
                if (it1 < it2) {
                    init_core2(&(cores[core_index]), it2+1, it2-it1+2, end-it2-1+offset, end-it1-1+offset);
                    core_index++;
                }

                // create RINT core
                it2 = it1 - 2 - middle_count;
                init_core2(&(cores[core_index]), it1, 2+middle_count, end-it1-1+offset, end-it2-1+offset);
                core_index++;

                continue;
            }
        }

        if (rc_alphabet[(unsigned char)*it1] > rc_alphabet[(unsigned char)*(it1-1)] &&
            rc_alphabet[(unsigned char)*(it1-1)] < rc_alphabet[(unsigned char)*(it1-2)]) {

            // check if there is any SSEQ cores left behind
            if (it1 < it2) {
                init_core2(&(cores[core_index]), it2+1, it2-it1+2, end-it2-1+offset, end-it1-1+offset);
                core_index++;
            }

            // create LMIN core
            it2 = it1 - 3;
            init_core2(&(cores[core_index]), it1, 3, end-it1-1+offset, end-it2-1+offset);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (begin <= it1-3 &&
            rc_alphabet[(unsigned char)*it1] < rc_alphabet[(unsigned char)*(it1-1)] &&
            rc_alphabet[(unsigned char)*(it1-1)] > rc_alphabet[(unsigned char)*(it1-2)] &&
            rc_alphabet[(unsigned char)*(it1+1)] <= rc_alphabet[(unsigned char)*(it1)] &&
            rc_alphabet[(unsigned char)*(it1-2)] >= rc_alphabet[(unsigned char)*(it1-3)]) {

            // check if there is any SSEQ cores left behind
            if (it1 < it2) {
                init_core2(&(cores[core_index]), it2+1, it2-it1+2, end-it2-1+offset, end-it1-1+offset);
                core_index++;
            }

            // create LMAX core
            it2 = it1 - 3;
            init_core2(&(cores[core_index]), it1, 3, end-it1-1+offset, end-it2-1+offset);
            core_index++;

            continue;
        }
    }

    return core_index;
}

/**
 * @brief Parses a array of cores to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores.
 *
 * This function iterates over a array of `core` structures defined by iterators `begin` and `end` and 
 * identifies key segments, called "cores," that represent the (LCP) regions. By analyzing
 * `core` structure relationships in the array (such as equality or relative order), it builds and stores
 * these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the `core` array to parse.
 * @param end Iterator pointing to the end of the `core` array to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @return Size of the cores identified in the given string.
 */
int parse3(struct core *begin, struct core *end, struct core *cores) {

    struct core *it1 = begin;
    struct core *it2 = end;
    int core_index = 0;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (core_eq(it1, it1+1)) {
            continue;
        }

        // check for RINT core
        if (core_eq(it1+1, it1+2)) {

            // count middle characters
            uint32_t middle_count = 1;
            struct core *temp = it1 + 2;
            while (temp < end && core_eq(temp-1, temp)) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1) {
                    init_core3(&(cores[core_index]), it2-1, it1-it2+2);
                    core_index++;
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                init_core3(&(cores[core_index]), it1, it2-it1);
                core_index++;

                continue;
            }
        }

        // check for LMIN
        if (core_gt(it1, it1+1) && core_lt(it1+1, it1+2)) {
            
            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core3(&(cores[core_index]), it2-1, it1-it2+2);
                core_index++;
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1, it2-it1);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            core_lt(it1, it1+1) &&
            core_gt(it1+1, it1+2) &&
            core_leq(it1-1, it1) &&
            core_geq(it1+2, it1+3)) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core3(&(cores[core_index]), it2-1, it1-it2+2);
                core_index++;
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1, it2-it1);
            core_index++;

            continue;
        }
    }
    return core_index;
}

void init_lps(struct lps *lps_ptr, const char *str, int len) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1(str, str+len, lps_ptr->cores, 0);
}

void init_lps_offset(struct lps *lps_ptr, const char *str, int len, uint64_t offset) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1(str, str+len, lps_ptr->cores, offset);
}

void init_lps2(struct lps *lps_ptr, const char *str, int len) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse2(str, str+len, lps_ptr->cores, 0);
}

void init_lps3(struct lps *lps_ptr, FILE *in) {
    // read the level from the binary file
    if (fread(&(lps_ptr->level), sizeof(int), 1, in) != 1) {
        fprintf(stderr, "Error reading level from file\n");
        exit(EXIT_FAILURE);
    }

    // read the size (number of cores)
    if(fread(&(lps_ptr->size), sizeof(int), 1, in) != 1) {
        fprintf(stderr, "Error reading size from file\n");
        exit(EXIT_FAILURE);
    }

    lps_ptr->cores = NULL;

    if (lps_ptr->size) {
        // allocate memory for the cores array
        lps_ptr->cores = (struct core *)malloc(lps_ptr->size * sizeof(struct core));

        // read each core object from the file
        for (int i = 0; i < lps_ptr->size; i++) {
            struct core *cr = &(lps_ptr->cores[i]);

            if (fread(&(cr->bit_size), sizeof(ubit_size), 1, in) != 1) {
                fprintf(stderr, "Error reading bit_size from file at %d\n", i);
                exit(EXIT_FAILURE);
            }
    
            ubit_size block_number = (cr->bit_size + UBLOCK_BIT_SIZE - 1) / UBLOCK_BIT_SIZE;
            cr->bit_rep = (ublock *)malloc(block_number * sizeof(ublock));
            if (fread(cr->bit_rep, block_number * sizeof(ublock), 1, in) != 1) {
                fprintf(stderr, "Error reading bit_rep from file at %d\n", i);
                exit(EXIT_FAILURE);
            }
         
            if (fread(&(cr->label), sizeof(lcp_label), 1, in) != 1) {
                fprintf(stderr, "Error reading label from file at %d\n", i);
                exit(EXIT_FAILURE);
            }
            if (fread(&(cr->start), sizeof(lcp_pos), 1, in) != 1) {
                fprintf(stderr, "Error reading start from file at %d\n", i);
                exit(EXIT_FAILURE);
            }
            if (fread(&(cr->end), sizeof(lcp_pos), 1, in) != 1) {
                fprintf(stderr, "Error reading end from file at %d\n", i);
                exit(EXIT_FAILURE);
            }
        }
    }
}

void init_lps4(struct lps *lps_ptr, const char *str, int len, int lcp_level, int chunk_size) {

    if (lcp_level < 1)
        return;

    lps_ptr->level = 1;
    lps_ptr->size = 0; 
    int estimated_size = (int)(len / pow((double)CONSTANT_FACTOR, lcp_level));
    lps_ptr->cores = (struct core *)malloc(estimated_size*sizeof(struct core));

    int str_index = 0, core_index = 0;

    {
        int str_len = minimum(chunk_size, len);
        struct lps temp_lps;
        init_lps_offset(&temp_lps, str, str_len, 0);
        lps_deepen(&temp_lps, lcp_level);

        if (temp_lps.size) {
            memcpy(lps_ptr->cores, temp_lps.cores, (temp_lps.size)*sizeof(struct core));
            core_index = (temp_lps.size);
            lps_ptr->size = (temp_lps.size);
            if (temp_lps.size>1)
                str_index = lps_ptr->cores[core_index-2].start;
            else 
                str_index = lps_ptr->cores[core_index-1].start;
        }
        free(temp_lps.cores);
    }

    while (str_index < len) {
        int str_len = minimum(chunk_size, len-str_index);
        struct lps temp_lps;
        init_lps_offset(&temp_lps, str+str_index, str_len, str_index);
        lps_deepen(&temp_lps, lcp_level);

        if (1<temp_lps.size) {
            int overlap = 2;
            while (0<overlap) {
                if (lps_ptr->cores[core_index-overlap].start == temp_lps.cores[0].start)
                    break;
                overlap--;
            }
            for(int i=0; i<overlap; i++) {
                free_core(&(temp_lps.cores[i]));
            }
            memcpy(lps_ptr->cores+core_index, temp_lps.cores+overlap, (temp_lps.size-overlap)*sizeof(struct core));
            core_index += (temp_lps.size-overlap);
            lps_ptr->size += (temp_lps.size-overlap);

            if ((lcp_pos)str_index < lps_ptr->cores[core_index-2].start) {
                str_index = lps_ptr->cores[core_index-2].start;
                free(temp_lps.cores);
                continue;
            } 
        }
        
        // find next start point
        for(int i=str_index+str_len-1; str_index <= i; i--) {
            if (alphabet[(unsigned char)*(str+i)] == -1) {
                str_index = i+1;
                break;
            }
        }
        if (alphabet[(unsigned char)*(str+str_index)] != -1) { // all of the characters are valid, so not valid cores found
            str_index += str_len;
        }
        
        free(temp_lps.cores);
    }

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));
}

void free_lps(struct lps *lps_ptr) {
    for(int i=0; i<lps_ptr->size; i++) {
        free(lps_ptr->cores[i].bit_rep);
    }
    free(lps_ptr->cores);
    lps_ptr->size = 0;
}

void write_lps(struct lps *lps_ptr, FILE *out) {
    // write the level field
    fwrite(&(lps_ptr->level), sizeof(int), 1, out);

    // write the size (number of cores)
    fwrite(&(lps_ptr->size), sizeof(int), 1, out);

    // write each core object iteratively
    if (lps_ptr->size) {
        for (int i = 0; i < lps_ptr->size; i++) {
            const struct core *cr = &(lps_ptr->cores[i]);

            fwrite(&(cr->bit_size), sizeof(ubit_size), 1, out);
            
            ubit_size block_number = (cr->bit_size + UBLOCK_BIT_SIZE - 1) / UBLOCK_BIT_SIZE;
            fwrite(cr->bit_rep, sizeof(ublock), block_number, out);
            
            fwrite(&(cr->label), sizeof(lcp_label), 1, out);
            fwrite(&(cr->start), sizeof(lcp_pos), 1, out);
            fwrite(&(cr->end), sizeof(lcp_pos), 1, out);
        }
    }
}

int64_t lps_memsize(const struct lps *lps_ptr) {
    uint64_t total = sizeof(struct lps);
    
    for(int i=0; i<lps_ptr->size; i++) {
        total += core_memsize(&(lps_ptr->cores[i]));
    }

    return total;
}

/**
 * @brief Performs Deterministic Coin Tossing (DCT) compression for a given number of iterations.
 *
 * This function applies DCT compression to the cores stored in the LPS structure. During each
 * iteration, consecutive cores are compared and compressed from right to left. The compression
 * method used depends on the current LPS level: level 1 uses level-1 core compression, while
 * higher levels use upper-level core compression.
 *
 * The number of DCT iterations controls how many leading cores are consumed before the remaining
 * cores are parsed into the next compression level.
 *
 * @param lps_ptr The `lps` object whose cores will be compressed.
 * @param dct_iteration_count The number of DCT iterations to perform.
 * @return 0 if DCT compression was performed, -1 if there are not enough cores.
 */
int lcp_dct_iters(struct lps *lps_ptr, int dct_iteration_count) {
    
    // at least 2 cores are needed for compression
    if (lps_ptr->size < dct_iteration_count + 1) {
        return -1;
    }

    for (int dct_index = 0; dct_index < dct_iteration_count; dct_index++) {
        struct core *it_left = lps_ptr->cores + lps_ptr->size - 2, *it_right = lps_ptr->cores + lps_ptr->size - 1;

        for (; lps_ptr->cores + dct_index <= it_left; it_left--, it_right--) {
            core_compress(it_left, it_right);
        }
    }

    return 0;
}

/**
 * @brief Performs Deterministic Coin Tossing (DCT) compression using the default iteration count.
 *
 * This function applies DCT compression to the cores stored in the LPS structure using
 * `LCP_DCT_ITERATION_COUNT`. It is a convenience wrapper around `lcp_dct_iters`.
 *
 * @param lps_ptr The `lps` object whose cores will be compressed.
 * @return 0 if DCT compression was performed, -1 if there are not enough cores.
 */
int lcp_dct(struct lps *lps_ptr) {
    return lcp_dct_iters(lps_ptr, LCP_DCT_ITERATION_COUNT);
}

int lps_deepen1_dct_iters(struct lps *lps_ptr, int dct_iteration_count) {

    // compress cores
    if (lcp_dct_iters(lps_ptr, dct_iteration_count) < 0) {
        lps_ptr->size = 0;
        lps_ptr->level++;
        return 0;
    }

    // find new cores
    int new_size = parse3(lps_ptr->cores + dct_iteration_count, lps_ptr->cores + lps_ptr->size, lps_ptr->cores);
    int temp = new_size;

    // remove old cores
    while(temp < lps_ptr->size) {
        temp++;
    }
    lps_ptr->size = new_size;

    lps_ptr->level++;

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));

    return 1;
}

int lps_deepen1(struct lps *lps_ptr) {
    return lps_deepen1_dct_iters(lps_ptr, LCP_DCT_ITERATION_COUNT);
}

int lps_deepen_dct_iters(struct lps *lps_ptr, int lcp_level, int dct_iteration_count) {
    if (lcp_level <= lps_ptr->level)
        return 0;

    while (lps_ptr->level < lcp_level && lps_deepen1_dct_iters(lps_ptr, dct_iteration_count))
        ;

    return 1;
}

int lps_deepen(struct lps *lps_ptr, int lcp_level) {
    return lps_deepen_dct_iters(lps_ptr, lcp_level, LCP_DCT_ITERATION_COUNT);
}

void print_lps(const struct lps *lps_ptr) {
    printf("Level: %d \n", lps_ptr->level);
    for(int i=0; i<lps_ptr->size; i++) {
        print_core(&(lps_ptr->cores[i]));
        printf(" ");
    }
}

int lps_eq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 0;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 0;
        }
    }

    return 1;
}

int lps_neq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 1;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 1;
        }
    }

    return 0;
}