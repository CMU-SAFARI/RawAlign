#include <stdint.h>
#include <string.h> //for memset

#ifndef ROPTIONS_H
#define ROPTIONS_H

#define RI_I_NAIVE		0x1
#define RI_I_MIN		0x2
#define RI_I_BLEND		0x4
#define RI_I_SYNCMER	0x8

#define RI_M_SEQUENCEUNTIL			0x1
#define RI_M_DTW_EVALUATE_CHAINS	0x2
#define RI_M_DTW_OUTPUT_CIGAR		0x4
#define RI_M_DTW_LOG_SCORES			0x8
#define RI_M_DISABLE_CHAININGSCORE_FILTERING 0x10
#define RI_M_OUTPUT_CHAINS			0x20
#define RI_M_LOG_ANCHORS			0x40
#define RI_M_LOG_NUM_ANCHORS		0x80

#define RI_M_DTW_BORDER_CONSTRAINT_GLOBAL	0
#define RI_M_DTW_BORDER_CONSTRAINT_SPARSE	1
#define RI_M_DTW_BORDER_CONSTRAINT_LOCAL	2

#define RI_M_DTW_FILL_METHOD_FULL		0
#define RI_M_DTW_FILL_METHOD_BANDED		1

#ifdef __cplusplus
extern "C" {
#endif

// indexing and mapping options
typedef struct ri_idxopt_s{
	short b, w, e, n, q, lq, k, flag;
	int64_t mini_batch_size;
	uint64_t batch_size;
} ri_idxopt_t;

typedef struct ri_mapopt_s{

	//ONT Device specific parameters
	uint32_t bp_per_sec;
	uint32_t sample_rate;
	uint32_t chunk_size;

	//Chaining parameters
	uint32_t min_events;
	uint32_t max_gap_length;
	uint32_t max_target_gap_length;
	uint32_t chaining_band_length;
	uint32_t max_num_skips;
	uint32_t min_num_anchors;
	uint32_t num_best_chains;
	float min_chaining_score;

	//Mapping parameters
	uint32_t step_size;
	uint32_t max_num_chunk;
	uint32_t min_chain_anchor;
	uint32_t min_chain_anchor_out;
	uint32_t dtw_border_constraint;
	uint32_t dtw_fill_method;
	float dtw_band_radius_frac;
	float dtw_match_bonus;
	float dtw_min_score;

	float min_bestmap_ratio;
	float min_bestmap_ratio_out;

	float min_meanmap_ratio;
	float min_meanmap_ratio_out;

	float t_threshold;
	uint32_t tn_samples;
	uint32_t ttest_freq;
	uint32_t tmin_reads;
	
	int64_t flag;    // see ri_F_* macros
	int64_t mini_batch_size; // size of a batch of query bases to process in parallel

	//Event detector options
	uint32_t window_length1;
	uint32_t window_length2;
	float threshold1;
	float threshold2;
	float peak_height;
} ri_mapopt_t;

/**
 * Initializes the default indexing options
 *
 * @param opt	pointer to the indexing options
 * 
 */
void ri_idxopt_init(ri_idxopt_t *opt);

/**
 * Initializes the default mapping options
 *
 * @param opt	pointer to the mapping options
 * 
 */
void ri_mapopt_init(ri_mapopt_t *opt);

#ifdef __cplusplus
}
#endif
#endif //ROPTIONS_H