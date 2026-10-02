/**
 * @file core.h
 * @brief Header file for the `core` struct, which represents bit-encoded
 * sequences and provides various utilities for encoding, compression, and file
 * I/O operations on bit sequences. (I/O operations will be supported later)
 *
 * This file contains the declaration of the `core` struct along with its member
 * functions, constructors, destructors, and operator overloads. It supports
 * operations for encoding strings into bit sequences, combining multiple
 * sequences, compressing data, and reading/writing from/to files.
 *
 * Key functionalities include:
 * - Encoding sequences of characters (e.g., ACGT) into a compact bit-encoded
 * form.
 * - Supporting reverse complement encoding of DNA sequences.
 * - Saving and loading core from files.
 * - Calculating memory usage of the constructed core structure.
 *
 * Dependencies:
 * - Requires encoding.h and config.h for auxiliary data structures and utilities.
 *
 * @see encoding.h
 * 
 * @struct core
 *
 * @author Akmuhammet Ashyralyyev
 * @version 2.1
 * @date 2026-08-11
 *
 */

#ifndef CORE_H
#define CORE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "encoding.h"
#include "config.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define UBLOCK_BIT_SIZE 32

#define minimum(a, b) ((a) < (b) ? (a) : (b))

typedef unsigned int ublock;
typedef uint32_t ubit_size;

#if LCP_LABEL_BITS == 64
typedef uint64_t lcp_label;
#elif LCP_LABEL_BITS == 32
typedef uint32_t lcp_label;
#else
#error "LCP_LABEL_BITS must be 32 or 64; regenerate config.h"
#endif

#if LCP_POS_BITS == 64
typedef uint64_t lcp_pos;
#elif LCP_POS_BITS == 32
typedef uint32_t lcp_pos;
#else
#error "LCP_POS_BITS must be 32 or 64; regenerate config.h"
#endif

struct core {
    lcp_pos start;      // start offset of lcp core within string
    lcp_pos end;        // end offset of lcp core within string
    ublock *bit_rep;    // open bit representation of lcp core alphabet
    lcp_label label;       // if od the core, computed from either base alphabet or sub-labels
    ubit_size bit_size; // number of bits stored in bit_rep
};

/**
 * @brief Initializes a core structure with the provided string data and index range.
 * 
 * This function processes a given substring starting at `begin` with a specified 
 * distance (length of the substring) and assigns start and end indices for tracking.
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param begin Pointer to the start of the string data.
 * @param distance Length of the substring to process.
 * @param start_index Start index of the substring within the data.
 * @param end_index End index of the substring within the data.
 */
void init_core1(struct core *cr, const char *begin, lcp_pos distance, lcp_pos start_index, lcp_pos end_index);

/**
 * @brief Initializes a core structure with the provided string data and index range.
 * 
 * Similar to `init_core1`, but initalizes the core structure with reverse complement
 * alphabet encoding.
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param begin Pointer to the start of the string data.
 * @param distance Length of the substring to process.
 * @param start_index Start index of the substring within the data.
 * @param end_index End index of the substring within the data.
 */
void init_core2(struct core *cr, const char *begin, lcp_pos distance, lcp_pos start_index, lcp_pos end_index);

/**
 * @brief Initializes a core structure by combining data from other core structures.
 * 
 * This function initializes a new core structure (`cr`) using a sequence of 
 * `core` objects starting from `begin` with the specified `distance` (number 
 * of `core` objects to process).
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param begin Pointer to the start of the sequence of core structures.
 * @param distance Number of core structures to process in the sequence.
 */
void init_core3(struct core *cr, struct core *begin, lcp_pos distance);

/**
 * @brief Directly initializes a core structure with precomputed representation and metadata.
 * 
 * This function allows the initialization of a core structure when the bit 
 * representation, bit size, and label are already computed. Useful for 
 * deserializing or cloning a core structure.
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param bit_size Size of the bit representation in bits.
 * @param bit_rep Pointer to the precomputed bit representation array.
 * @param label Unique label assigned to the core structure.
 * @param start Start index of the substring or sequence represented by the core.
 * @param end End index of the substring or sequence represented by the core.
 */
void init_core4(struct core *cr, ubit_size bit_size, ublock *bit_rep, lcp_label label, lcp_pos start, lcp_pos end);

/**
 * @brief Frees the allocated memory associated with a core structure.
 * 
 * This function ensures that all dynamically allocated resources (e.g., the bit 
 * representation array) associated with the given core structure are properly 
 * released, preventing memory leaks.
 * 
 * @param cr Pointer to the core structure to deallocate.
 */
void free_core(struct core* cr);

/**
 * @brief Compresses the right `core` object by comparing it with
 * left `core` object.
 *
 * This function compresses the core's bit sequence by identifying
 * common patterns between the right `core` object and left `core` object,
 * and reduces the size of the sequence accordingly.
 *
 * @param left_core The `core` object to compare against for compression.
 * @param right_core The `core` object that will be compressed.
 */
void core_compress(const struct core *left_core, struct core *right_core);

/**
 * @brief Calculates the total memory size used by the `core` object.
 *
 * This function computes the memory used by the `core` object,
 * including the bit sequence and metadata.
 *
 * @return The total memory size in bytes.
 */
uint64_t core_memsize(const struct core *cr);

/**
 * @brief Output the bit representation of a `core` object.
 *
 * This function iterates over the bits in the `core` object and outputs
 * them to the given output stream. It prints each bit as either 0 or 1,
 * starting from the most significant bit to the least significant bit.
 *
 * @param cr The `core` object to be output.
 */
void print_core(const struct core *cr);

// core comparison operator implementation

/**
 * @brief Operator overload for equality comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the two objects are equal, 0 otherwise.
 */
int core_eq(const struct core *lhs, const struct core *rhs);

/**
 * @brief Operator overload for greater-than comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is greater, 0 otherwise.
 */
int core_neq(const struct core *lhs, const struct core *rhs);

/**
 * @brief Operator overload for smaller-than comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is smaller, 0 otherwise.
 */
int core_gt(const struct core *lhs, const struct core *rhs);

/**
 * @brief Operator overload for not-equal-to comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the two objects are not equal, 0 otherwise.
 */
int core_lt(const struct core *lhs, const struct core *rhs);

/**
 * @brief Operator overload for greater-than-or-equal-to comparison between
 * two `core` objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is greater than or equal, 0 otherwise.
 */
int core_geq(const struct core *lhs, const struct core *rhs);

/**
 * @brief Operator overload for smaller-than-or-equal-to comparison between
 * two `core` objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is smaller than or equal, 0 otherwise.
 */
int core_leq(const struct core *lhs, const struct core *rhs);

#ifdef __cplusplus
}
#endif

#endif