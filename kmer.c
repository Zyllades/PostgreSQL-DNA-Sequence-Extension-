#include "postgres.h"
#include "libpq/pqformat.h"
#include "lib/stringinfo.h"
#include "utils.h"
#include "fmgr.h"
#include "kmer.h"
#include "qkmer.h"

/*****************************************************************************/
/* Helper Functions */


bool verify_kmer(uint8_t *sequence, int packed_length) {
    for (int i = 0; i < packed_length; i++) {
        uint8_t byte = sequence[i];
        // Extract the high nibble (first 4 bits)
        uint8_t high_nibble = (byte >> 4) & 0x0F;
        // Extract the low nibble (last 4 bits)
        uint8_t low_nibble = byte & 0x0F;
        // Verify both nibbles are valid DNA bases or padding
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

/**
 * Allocates and initializes a k-mer object from a nucleotide sequence.
 */
uint8_t *kmer_make(const char *input, int length) {
    int packed_length = (length + 1) / 2;
    uint8_t *sequence = (uint8_t *) malloc(packed_length * sizeof(uint8_t));
    // Convert and store the nucleotide sequence
    for (int i = 0; i < length; i += 2) {
        uint8_t nucleotide1 = nucleotideToBits(toupper(input[i]));
        uint8_t nucleotide2 = (i + 1 < length) ? nucleotideToBits(toupper(input[i + 1])) : PADDING;
        pack_nucleotides(sequence, i / 2, nucleotide1, nucleotide2);
    }
    if(verify_kmer(sequence, packed_length)){
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

PG_FUNCTION_INFO_V1(kmer_cast_from_text);
Datum
kmer_cast_from_text(PG_FUNCTION_ARGS)
{
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout,PointerGetDatum(txt)));
    int length = strlen(str);
    KMER *kmer = (KMER *) palloc(sizeof(KMER) + length * sizeof(uint8_t));
    kmer->length = length;
    memcpy(kmer->sequence, kmer_make(str, length), length * sizeof(uint8_t));
    SET_VARSIZE(kmer,sizeof(KMER) + length * sizeof(uint8_t));
    PG_FREE_IF_COPY(str,0);
    PG_RETURN_POINTER(kmer);

}

/*****************************************************************************/
/* k-mer Functions */

/**
 * PostgreSQL input function for k-mer type.
 */
PG_FUNCTION_INFO_V1(kmer_in);
Datum
kmer_in(PG_FUNCTION_ARGS) {
    char *str = PG_GETARG_CSTRING(0);
    int length = strlen(str);
    KMER *kmer = (KMER *) palloc(sizeof(KMER) + length * sizeof(uint8_t));
    kmer->length = length;
    memcpy(kmer->sequence, kmer_make(str, length), length * sizeof(uint8_t));
    SET_VARSIZE(kmer,sizeof(KMER) + length * sizeof(uint8_t));
    PG_FREE_IF_COPY(str,0);
    PG_RETURN_POINTER(kmer);
}
/**
 * PostgreSQL output function for DNA type.
 */
PG_FUNCTION_INFO_V1(kmer_out);
Datum
kmer_out(PG_FUNCTION_ARGS) {
    KMER *kmer =(KMER*) PG_GETARG_POINTER(0);
    char *result = sequence_to_string(kmer->sequence, kmer->length);
    PG_FREE_IF_COPY(kmer, 0);
    PG_RETURN_CSTRING(result);
}

/**
 * DNA constructor function.
 */

PG_FUNCTION_INFO_V1(kmer_constructor);
Datum
kmer_constructor(PG_FUNCTION_ARGS) {
    KMER *kmer;
    int32 length = PG_GETARG_SEQUENCE_LENGTH(0);
    char *str = PG_GETARG_CSTRING(1);
    if (length != strlen(str))
    {
        ereport(ERROR,
                    (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                     errmsg("Invalid declared length.")));
    }
    kmer = (KMER *) palloc(sizeof(KMER) + length * sizeof(uint8_t));
    kmer->length = length;
    memcpy(kmer->sequence, kmer_make(str, length), length * sizeof(uint8_t));
    SET_VARSIZE(kmer,sizeof(KMER) + length * sizeof(uint8_t));
    PG_FREE_IF_COPY(str,1);
    PG_RETURN_POINTER(kmer);
    }


/*****************************************************************************/

bool kmer_equals(KMER *kmer1, KMER *kmer2) {
    char * string1 = sequence_to_string(kmer1->sequence, kmer1->length);
    char * string2 = sequence_to_string(kmer2->sequence, kmer2->length);
    if (kmer1->length != kmer2->length) {
//        elog(INFO, "Length differ, returning 0");
        return false;
    }
    for (int i = 0; i < kmer1->length; i++) {
        if (string1[i] != string2[i]) {
  //          elog(INFO, "symbols differ, returning 0: %c vs %c", string1[i], string2[i]);
            return false; // Sequences differ
        }
    }
    elog(INFO, "Sequences are the same, returning 1");
    return true;
}



/*KMER equals*/
PG_FUNCTION_INFO_V1(equals);
Datum 
equals(PG_FUNCTION_ARGS)
{
  KMER *kmer_1 =(KMER *) PG_GETARG_POINTER(0);
  KMER *kmer_2 =(KMER *) PG_GETARG_POINTER(1);
  bool result = kmer_equals(kmer_1,kmer_2);
  PG_FREE_IF_COPY(kmer_1,0);
  PG_FREE_IF_COPY(kmer_2,1);
  PG_RETURN_BOOL(result);
}

/*generate KMER*/
PG_FUNCTION_INFO_V1(generate_kmers);

Datum
generate_kmers(PG_FUNCTION_ARGS)
{
    FuncCallContext *funcctx;
    TupleDesc tuple_desc;
    AttInMetadata *attinmeta;
    MemoryContext oldcontext;
    int k;
    text *dna_text;
    char *dna_sequence;
    int dna_length;
    char *context_data;
    KMER *kmer;
    Datum values[1];
    bool nulls[1] = {false};    
    HeapTuple tuple;
    

    if (SRF_IS_FIRSTCALL())
    {
        funcctx = SRF_FIRSTCALL_INIT();
        oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

        // Extract input arguments
        dna_text = PG_GETARG_TEXT_P(0);
        k = PG_GETARG_INT32(1);

        // Convert DNA sequence to C string
        dna_sequence = text_to_cstring(dna_text);
        dna_length = strlen(dna_sequence);

        // Validate input
        if (k <= 0 || k > dna_length)
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                     errmsg("Invalid k value. It must be > 0 and <= length of the DNA sequence.")));

        // Prepare the user context
        
        context_data = palloc(dna_length + 1 + sizeof(int));
        strcpy(context_data, dna_sequence);
        memcpy(&context_data[dna_length + 1], &k, sizeof(int));
        funcctx->user_fctx = context_data;
        funcctx->max_calls = dna_length - k + 1;

        // Prepare the tuple descriptor
        if (get_call_result_type(fcinfo, NULL, &tuple_desc) != TYPEFUNC_COMPOSITE)
            ereport(ERROR,
                    (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                     errmsg("Function must return a composite type.")));
        attinmeta = TupleDescGetAttInMetadata(tuple_desc);
        funcctx->attinmeta = attinmeta;

        MemoryContextSwitchTo(oldcontext);
    }

    funcctx = SRF_PERCALL_SETUP();

    if (funcctx->call_cntr < funcctx->max_calls)
    {
        // Context retrieval
        int call_idx;
        char *kmer_substr;
        dna_sequence = (char *)funcctx->user_fctx;
        k = *((int *)(dna_sequence + strlen(dna_sequence) + 1));
        call_idx = funcctx->call_cntr;
        // Extract k-mer
        kmer_substr = palloc(k + 1);
        memcpy(kmer_substr, &dna_sequence[call_idx], k);
        kmer_substr[k] = '\0';
        // Create KMER object
        kmer = (KMER *)palloc(sizeof(KMER) + k * sizeof(uint8_t));
        SET_VARSIZE(kmer, sizeof(KMER) + k * sizeof(uint8_t));
        kmer->length = k;
        memcpy(kmer->sequence, kmer_make(kmer_substr, k), k);
        // Build tuple
        values[0] = PointerGetDatum(kmer);

        tuple_desc = funcctx->attinmeta->tupdesc;
        tuple = heap_form_tuple(tuple_desc, values, nulls);

        // Free temporary memory
        pfree(kmer_substr);
        pfree(kmer);

        // Return tuple
        SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(tuple));
    }
    else
    {
        SRF_RETURN_DONE(funcctx);
    }
}


/* Helper function for starts_with logic */
bool kmer_starts_with(KMER *kmer, KMER *prefix)
{
    uint8_t prefix_last;
    uint8_t kmer_last;
    int full_bytes;
    int is_odd;
    full_bytes = prefix->length / 2;
    is_odd = prefix->length % 2;
    if (prefix->length > kmer->length)
    {
        return false;
    }
    if (full_bytes != 0)
    {
        if (memcmp(prefix->sequence, kmer->sequence, full_bytes) != 0)
        {
            return false;
        }
    }
    if (is_odd)
    {
        prefix_last = (prefix->sequence[full_bytes] >> 4) & 0xF;
        kmer_last = (kmer->sequence[full_bytes] >> 4) & 0xF;
        if (prefix_last != kmer_last)
        {
            return false;
        }
    }

    return true; // Prefix matches
}

/* KMER starts_with */
PG_FUNCTION_INFO_V1(starts_with);
Datum
starts_with(PG_FUNCTION_ARGS)
{
    KMER *kmer = (KMER *) PG_GETARG_POINTER(0);
    KMER *prefix = (KMER *) PG_GETARG_POINTER(1);

    // Use a helper function to encapsulate the starts_with logic
    bool result = kmer_starts_with(kmer,prefix);

    PG_FREE_IF_COPY(kmer, 0);
    PG_FREE_IF_COPY(prefix, 1);
    

    PG_RETURN_BOOL(result);
}

bool iupac_matches(uint8_t base, uint8_t pattern)
{
    if (pattern == 0xE)
        return true;
    switch (pattern)
    {
        case 0x0: return base == 0x0; // A
        case 0x1: return base == 0x1; // C
        case 0x2: return base == 0x2; // G
        case 0x3: return base == 0x3; // T
        case 0x4: return base == 0x0 || base == 0x2; // R (A, G)
        case 0x5: return base == 0x1 || base == 0x3; // Y (C, T)
        case 0x6: return base == 0x2 || base == 0x1; // S (G, C)
        case 0x7: return base == 0x0 || base == 0x3; // W (A, T)
        case 0x8: return base == 0x2 || base == 0x3; // K (G, T)
        case 0x9: return base == 0x0 || base == 0x1; // M (A, C)
        case 0xA: return base == 0x1 || base == 0x2 || base == 0x3; // B (C, G, T)
        case 0xB: return base == 0x0 || base == 0x2 || base == 0x3; // D (A, G, T)
        case 0xC: return base == 0x0 || base == 0x1 || base == 0x3; // H (A, C, T)
        case 0xD: return base == 0x0 || base == 0x1 || base == 0x2; // V (A, C, G)
        default: return false;
    }
}


/* KMER contains */
PG_FUNCTION_INFO_V1(contains);

Datum
contains(PG_FUNCTION_ARGS)
{
    QKMER *qkmer = (QKMER *) PG_GETARG_POINTER(0);
    KMER *kmer = (KMER *) PG_GETARG_POINTER(1);

    // Use a helper function to encapsulate the contains logic
    bool result = kmer_contains(qkmer, kmer);

    PG_FREE_IF_COPY(qkmer, 0);
    PG_FREE_IF_COPY(kmer, 1);

    PG_RETURN_BOOL(result);
}

/* Helper function for contains logic */
bool kmer_contains(QKMER *qkmer, KMER *kmer)
{
    uint8_t qkmer_base, kmer_base;
    int start, i;
    bool match;
    if (qkmer->length > kmer->length)
    {
        return false;
    }
    for (start = 0; start <= kmer->length - qkmer->length; start++)
    {
        match = true;
        for (i = 0; i < qkmer->length; i++)
        {
            kmer_base = (kmer->sequence[(start + i) / 2] >> (4 * ((start + i + 1) % 2))) & 0xF;
            qkmer_base = (qkmer->sequence[i / 2] >> (4 * ((i + 1) % 2))) & 0xF;
            if (!iupac_matches(kmer_base, qkmer_base))
            {
                match = false;
                break;
            }
        }
   if (match)
        {
            return true;
        }
    }
    return false;
}

/* KMER contains */
PG_FUNCTION_INFO_V1(kcontains);

Datum
kcontains(PG_FUNCTION_ARGS)
{
    KMER *kmer1 = (KMER *) PG_GETARG_POINTER(0);
    KMER *kmer2 = (KMER *) PG_GETARG_POINTER(1);

    // Use a helper function to encapsulate the contains logic
    bool result = k_ta_mer_contains(kmer1, kmer2);

    PG_FREE_IF_COPY(kmer1, 0);
    PG_FREE_IF_COPY(kmer2, 1);

    PG_RETURN_BOOL(result);
}


bool k_ta_mer_contains(KMER *kmer1, KMER *kmer2)
{
    uint8_t kmer1_base;
    uint8_t kmer2_base;
    int start;
    int i;
    bool match;

    // Ensure the query k-mer is not longer than the k-mer
    if (kmer1->length > kmer2->length)
    {
        return false;
    }

    // Iterate through possible starting positions in the k-mer
    for (start = 0; start <= kmer2->length - kmer1->length; start++)
    {
        match = true;

        // Check if the pattern matches this substring
        for (i = 0; i < kmer1->length; i++)
        {
            kmer2_base = (kmer2->sequence[(start + i) / 2] >> (4 * ((start + i + 1) % 2))) & 0xF;
            kmer1_base = (kmer1->sequence[i / 2] >> (4 * ((i + 1) % 2))) & 0xF;

            if (!iupac_matches(kmer2_base, kmer1_base))
            {
                match = false;
                break;
            }
        }

        // If a match is found, return true
        if (match)
        {
            return true;
        }
    }

    // No match found
    return false;
}