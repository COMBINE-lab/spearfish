/**
 * @file encoding.h
 * @brief Provides utilities for encoding and reverse complementing DNA
 * sequences using a character-to-bit mapping system.
 *
 * This file defines functions and variables for managing the encoding of DNA
 * sequences into bit-encoded values, as well as their reverse complement
 * encodings. It includes functionality for initializing encoding coefficients,
 * displaying the encoding summary, and loading custom encodings from
 * user-specified files or maps.
 *
 * Example usage:
 * @code
 *   // Initialize standard encoding and reverse complements
 *   LCP_INIT();
 *
 *   // Initialize encoding from a file
 *   LCP_INIT_FILE("encoding_test.txt", 0);
 * @endcode
 *
 * @see core.h
 *
 *
 * @author Akmuhammet Ashyralyyev
 * @version 2.1
 * @date 2026-08-11
 *
 */

#ifndef ENCODING_H
#define ENCODING_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#define maximum(a, b) ((a) > (b) ? (a) : (b))

extern int alphabet[128];
extern int rc_alphabet[128];
extern char characters[128];
extern int alphabet_bit_size;

/**
 * @brief Displays the alphabet encoding summary including coefficients
 * and dictionary bit size.
 */
void LCP_SUMMARY();

/**
 * @brief Initializes the encoding coefficients for standard DNA bases
 * (A, C, G, T) and their reverse complements. Sets default values for
 * coefficients and dictionary bit size.
 */
void LCP_INIT();

/**
 * @brief Initializes the encoding coefficients for standard DNA bases
 * (A, C, G, T) and their reverse complements. Sets default values for
 * coefficients and dictionary bit size.
 * @param verbose If 0, prints the encoding summary after initialization.
 */
void LCP_INIT2(int verbose);

/**
 * @brief Initializes the encoding coefficients by reading them from a
 * file. The file must contain character, encoding, and reverse
 * complement encoding on each line.
 * @param filename Path to the file containing the character encodings.
 * @param verbose If true, prints the encoding summary after
 * initialization.
 * @return Always returns 0 upon successful initialization.
 * @throws std::invalid_argument if any invalid data is found in the
 * file.
 */
int LCP_INIT_FILE(const char *filename, int verbose);

#ifdef __cplusplus
}
#endif

#endif
