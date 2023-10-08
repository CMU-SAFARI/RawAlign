#include "roptions.h"

#include <stdlib.h>

void ri_mapopt_init(ri_mapopt_t *opt)
{
	memset(opt, 0, sizeof(ri_mapopt_t));

	opt->bp_per_sec = 450;
	opt->sample_rate = 4000;
	opt->chunk_size = 4000;

	//chaining
	opt->max_gap_length = 2000;
	opt->max_target_gap_length = 5000;
	opt->chaining_band_length = 5000;
	opt->max_num_skips = 25;
	opt->min_num_anchors = 2;
	opt->num_best_chains = 3;
	opt->min_chaining_score = 10.0f;

	opt->step_size = 1; //read_seeding_step_size
	opt->min_events = 50;
	opt->max_num_chunk = 30;//max_num_chunks
	opt->min_chain_anchor = 2; //stop_mapping_min_num_anchors
	opt->min_chain_anchor_out = 2;//output_mapping_min_num_anchors

	opt->min_bestmap_ratio = 1.2f; //stop_mapping_ratio
	opt->min_bestmap_ratio_out = 1.2f; //output_mapping_ratio

	opt->min_meanmap_ratio = 5; //stop_mapping_mean_ratio
	opt->min_meanmap_ratio_out = 5; //output_mapping_mean_ratio

	opt->mini_batch_size = 500000000;

	//Default options for event detection. TODO: Make it flexible so that we can change them to RNA values as well
	opt->window_length1 = 3;
    opt->window_length2 = 6;
    opt->threshold1 = 4.30265f;  // 4.60409f,//1.4f,
    opt->threshold2 = 2.57058f;  // 3.16927f,//9.0f,
    opt->peak_height = 1.0f;      // 0.2f

	opt->t_threshold = 1.5f;
	opt->tn_samples = 5;
	opt->ttest_freq = 500;
	opt->tmin_reads = 500;

	//dtw
	opt->dtw_border_constraint = RI_M_DTW_BORDER_CONSTRAINT_SPARSE;
	opt->dtw_fill_method = RI_M_DTW_FILL_METHOD_BANDED;
	opt->dtw_band_radius_frac = 0.10f;
	opt->dtw_match_bonus = 0.4f;
	opt->dtw_min_score = 20.0f;

	//TODO: RNA values:
	// opt->window_length1 = 7,
	// opt->window_length2 = 14,
	// opt->threshold1 = 2.5f,
	// opt->threshold2 = 9.0f,
	// opt->peak_height = 1.0f;
}
