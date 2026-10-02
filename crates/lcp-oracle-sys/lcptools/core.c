#include "core.h"


/**
 * @brief Computes the 32-bit MurmurHash3 hash for a given key.
 *
 * This function computes a 32-bit hash of the input data 'key' with the
 * specified length 'len' and an optional seed value. It processes the
 * input in blocks and handles any remaining bytes.
 *
 * @param key Pointer to the data to be hashed.
 * @param len The length of the data in bytes.
 * @param seed An initial seed value for the hash computation.
 * @return The resulting 32-bit hash value.
 */
static inline uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed) {
    const uint8_t *data = (const uint8_t *)key;
    const int nblocks = len / 4;

    uint32_t h1 = seed;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    // Body: Process blocks of 4 bytes at a time
    const uint32_t *blocks = (const uint32_t *)(data + nblocks * 4);

    for (int i = -nblocks; i; i++) {
        uint32_t k1 = blocks[i];

        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;

        h1 ^= k1;
        h1 = (h1 << 15) | (h1 >> (32 - 15));
        h1 = h1 * 5 + 0xe6546b64;
    }

    // Tail: Process remaining bytes
    const uint8_t *tail = (const uint8_t *)(data + nblocks * 4);

    uint32_t k1 = 0;

    switch (len & 3) {
    case 3:
        k1 ^= tail[2] << 16;
        break;
    case 2:
        k1 ^= tail[1] << 8;
        break;
    case 1:
        k1 ^= tail[0];
        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;
        h1 ^= k1;
    }

    // Finalization: Mix the hash to ensure the last few bits are fully mixed
    h1 ^= len;

    /* fmix32 */
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    return h1;
}

void init_core1(struct core *cr, const char *begin, lcp_pos distance, lcp_pos start_index, lcp_pos end_index) {

    cr->start = start_index;
    cr->end = end_index;
    cr->bit_size = alphabet_bit_size * distance;

    /* allocate memory for representation */
    ubit_size block_number = (cr->bit_size + UBLOCK_BIT_SIZE - 1) / UBLOCK_BIT_SIZE;
    cr->bit_rep = (ublock *)malloc(block_number * sizeof(ublock));
    memset(cr->bit_rep, 0, block_number * sizeof(ublock));

    ubit_size shift = 0;
    int block_index = block_number - 1;

    for (const char *it = begin + distance - 1; begin <= it; it--) {

        int o_bit_rep = alphabet[(int)(*it)];

        /* shift and paste */
        cr->bit_rep[block_index] |= (o_bit_rep << shift);

        /* if there is an overflow after shifting, it pastes the */
        /* overfloaw to the left block. */
        if (shift + alphabet_bit_size > UBLOCK_BIT_SIZE) {
            cr->bit_rep[block_index - 1] |= (o_bit_rep >> (UBLOCK_BIT_SIZE - shift));
        }

        if (shift + alphabet_bit_size >= UBLOCK_BIT_SIZE) {
            block_index--;
        }

        shift = (shift + alphabet_bit_size) % UBLOCK_BIT_SIZE;
    }
    
    cr->label = 0;
    cr->label |= ((distance-2) << (3 * alphabet_bit_size));
    cr->label |= (alphabet[(*(begin)) & 0xDF] << (2 * alphabet_bit_size));
    cr->label |= (alphabet[(*(begin+distance-2)) & 0xDF] << alphabet_bit_size);
    cr->label |= (alphabet[(*(begin+distance-1)) & 0xDF]);
}

void init_core2(struct core *cr, const char *begin, lcp_pos distance, lcp_pos start_index, lcp_pos end_index) {

    cr->start = start_index;
    cr->end = end_index;
    cr->bit_size = alphabet_bit_size * distance;

    /* allocate memory for representation */
    ubit_size block_number = (cr->bit_size + UBLOCK_BIT_SIZE - 1) / UBLOCK_BIT_SIZE;
    cr->bit_rep = (ublock *)malloc(block_number * sizeof(ublock));
    memset(cr->bit_rep, 0, block_number * sizeof(ublock));

    ubit_size shift = 0;
    int block_index = block_number - 1;

    for (const char *it = begin - distance + 1; it <= begin; it++) {
        int o_bit_rep = rc_alphabet[(int)(*it)];

        /* shift and paste */
        cr->bit_rep[block_index] |= (o_bit_rep << shift);

        /* if there is an overflow after shifting, it pastes the */
        /* overfloaw to the left block. */
        if (shift + alphabet_bit_size > UBLOCK_BIT_SIZE) {
            cr->bit_rep[block_index - 1] |= (o_bit_rep >> (UBLOCK_BIT_SIZE - shift));
        }

        if (shift + alphabet_bit_size >= UBLOCK_BIT_SIZE) {
            block_index--;
        }

        shift = (shift + alphabet_bit_size) % UBLOCK_BIT_SIZE;
    }

    cr->label = 0;
    cr->label |= ((distance-2) << (3 * alphabet_bit_size));
    cr->label |= (rc_alphabet[(*(begin)) & 0xDF] << (2 * alphabet_bit_size));
    cr->label |= (rc_alphabet[(*(begin+distance-2)) & 0xDF] << alphabet_bit_size);
    cr->label |= (rc_alphabet[(*(begin+distance-1)) & 0xDF]);
}

void init_core3(struct core *cr, struct core *begin, lcp_pos distance) {

    // it is known that other core is placed in cr
    free(cr->bit_rep);

    cr->start = begin->start;
    cr->end = (begin+distance-1)->end;
    cr->bit_size = 0;

    for (struct core *it = begin; it < begin + distance; it++) {
        cr->bit_size += it->bit_size;
    }

    /* allocate memory for representation */
    ubit_size block_number = (cr->bit_size + UBLOCK_BIT_SIZE - 1) / UBLOCK_BIT_SIZE;
    cr->bit_rep = (ublock *)malloc(block_number * sizeof(ublock));
    memset(cr->bit_rep, 0, block_number * sizeof(ublock));

    ubit_size shift = 0;
    int block_index = block_number - 1;

    for (struct core *it = begin + distance - 1; begin <= it; it--) {

        ublock *o_bit_rep = it->bit_rep;

        for (int i = (it->bit_size - 1) / UBLOCK_BIT_SIZE; 0 <= i; i--) {

            ubit_size curr_block_size = (i > 0 ? UBLOCK_BIT_SIZE : it->bit_size % UBLOCK_BIT_SIZE);

            /* shift and paste */
            cr->bit_rep[block_index] |= (o_bit_rep[i] << shift);

            /* if there is an overflow after shifting, it pastes the */
            /* overfloaw to the left block. */
            if (shift + curr_block_size > UBLOCK_BIT_SIZE) {
                cr->bit_rep[block_index - 1] |= (o_bit_rep[i] >> (UBLOCK_BIT_SIZE - shift));
            }

            if (shift + curr_block_size >= UBLOCK_BIT_SIZE) {
                block_index--;
            }

            shift = (shift + curr_block_size) % UBLOCK_BIT_SIZE;
        }
    }

    lcp_label data[4];
    data[0] = (begin)->label;
    data[1] = (begin+distance-2)->label;
    data[2] = (begin+distance-1)->label;
    data[3] = distance-2;
    cr->label = MurmurHash3_32((void*)data, 4 * sizeof(lcp_label), 42);
}

void init_core4(struct core *cr, ubit_size bit_size, ublock *bit_rep, lcp_label label, lcp_pos start, lcp_pos end) {
    cr->bit_size = bit_size;
    cr->bit_rep = bit_rep;
    cr->label = label;
    cr->start = start;
    cr->end = end;
}

void free_core(struct core* cr) {
    free(cr->bit_rep);
    cr->bit_rep = NULL;
}

void core_compress(const struct core *left_core, struct core *right_core) {

    ubit_size index = minimum(left_core->bit_size, right_core->bit_size);
    ubit_size left_block = (left_core->bit_size - 1) / UBLOCK_BIT_SIZE,
              right_block = (right_core->bit_size - 1) / UBLOCK_BIT_SIZE;

    while (index >= UBLOCK_BIT_SIZE && left_core->bit_rep[left_block] == right_core->bit_rep[right_block]) {
        left_block--;
        right_block--;
        index -= UBLOCK_BIT_SIZE;
    }

    left_block = left_core->bit_rep[left_block];
    right_block = right_core->bit_rep[right_block];

    while (index > 0 && left_block % 2 == right_block % 2) {
        left_block /= 2;
        right_block /= 2;
        index--;
    }

    if (right_core->bit_size > UBLOCK_BIT_SIZE) {
        free(right_core->bit_rep);
        right_core->bit_rep = (ublock *)malloc(sizeof(ublock));
    }

    right_core->bit_rep[0] = 0;

    // shift left by 1 bit and set last bit to difference
    right_core->bit_rep[0] = 2 * (minimum(right_core->bit_size, left_core->bit_size) - index) + (right_block % 2);
    right_core->bit_size = 0;

    if (right_core->bit_rep[0] > 0) {
        right_core->bit_size = (32 - __builtin_clz(right_core->bit_rep[0]));
    }

    right_core->bit_size = right_core->bit_size > 1 ? right_core->bit_size : 2;

    // now, the right core is dependent on the left; hence, its coverage spans towards the left
    right_core->start = left_core->start;
}

uint64_t core_memsize(const struct core *cr) {
    return sizeof(struct core) + sizeof(ublock) * ((cr->bit_size + UBLOCK_BIT_SIZE - 1) / UBLOCK_BIT_SIZE);
}

void print_core(const struct core *cr) {
    ubit_size block_number = (cr->bit_size - 1) / UBLOCK_BIT_SIZE + 1;
    for (int index = cr->bit_size - 1; 0 <= index; index--) {
        printf("%d", (cr->bit_rep[block_number - index / UBLOCK_BIT_SIZE - 1] >> (index % UBLOCK_BIT_SIZE)) & 1);
    }
}

// core comparison operator implementation

int core_eq(const struct core *lhs, const struct core *rhs) {

    if (lhs->bit_size != rhs->bit_size) {
        return 0;
    }

    ubit_size index = 0;

    while (index < lhs->bit_size) {
        if (lhs->bit_rep[index / UBLOCK_BIT_SIZE] != rhs->bit_rep[index / UBLOCK_BIT_SIZE])
            return 0;

        index += UBLOCK_BIT_SIZE;
    }

    return 1;
}

int core_neq(const struct core *lhs, const struct core *rhs) {

    if (lhs->bit_size != rhs->bit_size) {
        return 1;
    }

    ubit_size index = 0;

    while (index < lhs->bit_size) {
        if (lhs->bit_rep[index / UBLOCK_BIT_SIZE] != rhs->bit_rep[index / UBLOCK_BIT_SIZE])
            return 1;

        index += UBLOCK_BIT_SIZE;
    }

    return 0;
}

int core_gt(const struct core *lhs, const struct core *rhs) {

    if (lhs->bit_size != rhs->bit_size) {
        return lhs->bit_size > rhs->bit_size;
    }

    ubit_size index = 0;

    while (index < lhs->bit_size) {
        if (lhs->bit_rep[index / UBLOCK_BIT_SIZE] == rhs->bit_rep[index / UBLOCK_BIT_SIZE]) {
            index += UBLOCK_BIT_SIZE;
            continue;
        }

        return lhs->bit_rep[index / UBLOCK_BIT_SIZE] > rhs->bit_rep[index / UBLOCK_BIT_SIZE];
    }

    return 0;
}

int core_lt(const struct core *lhs, const struct core *rhs) {

    if (lhs->bit_size != rhs->bit_size) {
        return lhs->bit_size < rhs->bit_size;
    }

    ubit_size index = 0;

    while (index < lhs->bit_size) {
        if (lhs->bit_rep[index / UBLOCK_BIT_SIZE] == rhs->bit_rep[index / UBLOCK_BIT_SIZE]) {
            index += UBLOCK_BIT_SIZE;
            continue;
        }

        return lhs->bit_rep[index / UBLOCK_BIT_SIZE] < rhs->bit_rep[index / UBLOCK_BIT_SIZE];
    }

    return 0;
}

int core_geq(const struct core *lhs, const struct core *rhs) {

    if (lhs->bit_size != rhs->bit_size) {
        return lhs->bit_size >= rhs->bit_size;
    }

    ubit_size index = 0;

    while (index < lhs->bit_size) {
        if (lhs->bit_rep[index / UBLOCK_BIT_SIZE] != rhs->bit_rep[index / UBLOCK_BIT_SIZE]) {
            return lhs->bit_rep[index / UBLOCK_BIT_SIZE] >= rhs->bit_rep[index / UBLOCK_BIT_SIZE];
        }

        index += UBLOCK_BIT_SIZE;
    }

    return 1;
}

int core_leq(const struct core *lhs, const struct core *rhs) {

    if (lhs->bit_size != rhs->bit_size) {
        return lhs->bit_size <= rhs->bit_size;
    }

    ubit_size index = 0;

    while (index < lhs->bit_size) {
        if (lhs->bit_rep[index / UBLOCK_BIT_SIZE] != rhs->bit_rep[index / UBLOCK_BIT_SIZE]) {
            return lhs->bit_rep[index / UBLOCK_BIT_SIZE] <= rhs->bit_rep[index / UBLOCK_BIT_SIZE];
        }

        index += UBLOCK_BIT_SIZE;
    }

    return 1;
}