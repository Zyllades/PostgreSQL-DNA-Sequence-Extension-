#include "postgres.h"
#include "libpq/pqformat.h"
#include "lib/stringinfo.h"
#include "utils.h"
#include "fmgr.h"
#include "dna.h"


/*****************************************************************************/
/* Helper Functions */

/**
 * Allocates and initializes a DNA object from a nucleotide sequence.
 * Ensures the sequence contains only valid nucleotides: A, T, G, and C.
 */



bool verify_dna(uint8_t *sequence, int packed_length) {
    int i;                       
    uint8_t byte, high_nibble, low_nibble;
    for (i = 0; i < packed_length; i++) {
        byte = sequence[i];
        high_nibble = (byte >> 4) & 0x0F;
        low_nibble = byte & 0x0F;
        if (high_nibble != A && high_nibble != T && high_nibble != G &&
            high_nibble != C && high_nibble != PADDING) {
            return false;
        }
        if (low_nibble != A && low_nibble != T && low_nibble != G &&
            low_nibble != C && low_nibble != PADDING) {
            return false;
        }
    }
    return true;
}

static uint8_t *dna_make(const char *input, int length) {
    int packed_length = (length + 1) / 2;
    uint8_t *sequence = (uint8_t *) malloc(packed_length * sizeof(uint8_t));
    for (int i = 0; i < length; i += 2) {
        uint8_t nucleotide1 = nucleotideToBits(toupper(input[i]));
        uint8_t nucleotide2 = (i + 1 < length) ? nucleotideToBits(toupper(input[i + 1])) : PADDING;
        pack_nucleotides(sequence, i / 2, nucleotide1, nucleotide2);
    }
    if(verify_dna(sequence, packed_length)){
        return sequence;
    }
    else{
        ereport(ERROR,
                    (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                     errmsg("Invalid nucleotide. "
                            "DNA sequence must only contain A, T, G, and C.")));
    }
    return NULL;
}


/*****************************************************************************/
/* DNA Functions */

/**
 * PostgreSQL input function for DNA type.
 */
PG_FUNCTION_INFO_V1(dna_in);
Datum
dna_in(PG_FUNCTION_ARGS) {
    char *str = PG_GETARG_CSTRING(0);
    int length = strlen(str);
    DNA *dna = (DNA *) palloc(sizeof(DNA) + (length+1)/2 * sizeof(uint8_t));
    dna->length = length;
    memcpy(dna->sequence, dna_make(str, length), ((length+1)/2 * sizeof(uint8_t)));
    SET_VARSIZE(dna,sizeof(DNA) + (length+1)/2 * sizeof(uint8_t));
    PG_FREE_IF_COPY(str,0);
    PG_RETURN_POINTER(dna);
}

/**
 * PostgreSQL output function for DNA type.
 */
PG_FUNCTION_INFO_V1(dna_out);
Datum
dna_out(PG_FUNCTION_ARGS) {
    char *result;
    // Retrieve the DNA object
    DNA *dna =(DNA*) PG_GETARG_POINTER(0);
    elog(INFO, "DNA received: %d", dna->length); 
    // Use the generalized `sequence_out` function
    result = sequence_to_string(dna->sequence, dna->length);
    // Free memory if DNA is a copy
    PG_FREE_IF_COPY(dna, 0);

    // Return the string
    PG_RETURN_CSTRING(result);
    
}

/**
 * DNA constructor function.
 */

PG_FUNCTION_INFO_V1(dna_constructor);
Datum
dna_constructor(PG_FUNCTION_ARGS) {
    DNA *dna ;
    // Retrieve the length and sequence arguments
    int32 length = PG_GETARG_SEQUENCE_LENGTH(0);
    // Retrieve the input as a C string
    char *str = PG_GETARG_CSTRING(1);
    if (length != strlen(str))
    {
        ereport(ERROR,
                    (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                     errmsg("Invalid declared length.")));
    }
    // Allocate memory for the DNA object
    dna = (DNA *) palloc(sizeof(DNA) + length * sizeof(uint8_t));
    dna->length = length;
    // Call dna_make directly to handle nucleotide validation and conversion
    memcpy(dna->sequence, dna_make(str, length), length * sizeof(uint8_t));
    SET_VARSIZE(dna,sizeof(DNA) + length * sizeof(uint8_t));
    //elog(INFO, "DNA dna_in length: %d", dna->length);
    //elog(INFO, "DNA dna_in sequence first: %d", dna->sequence[0]);
    PG_FREE_IF_COPY(str,1);
    // Return the DNA object
    PG_RETURN_POINTER(dna);
    }