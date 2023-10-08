#include "rmap.h"
#include "kthread.h"
#include "kvec.h"
#include "rutils.h"
#include "rsketch.h"
#include "revent.h"
#include "sequence_until.h"
#include "dtw.hpp"
#include <algorithm> //for sort
#include <math.h>
#include <sstream>

static ri_tbuf_t *ri_tbuf_init(void)
{
	ri_tbuf_t *b;
	b = (ri_tbuf_t*)calloc(1, sizeof(ri_tbuf_t));
	b->km = ri_km_init();
	return b;
}

static void ri_tbuf_destroy(ri_tbuf_t *b)
{
	if (b == 0) return;
	ri_km_destroy(b->km);
	free(b);
}

// void ri_seed_mz_flt(void *km, mm128_v *mv, int32_t q_occ_max, float q_occ_frac)
// {
// 	mm128_t *a;
// 	size_t i, j, st;
// 	if (mv->n <= q_occ_max || q_occ_frac <= 0.0f || q_occ_max <= 0) return;
// 	KMALLOC(km, a, mv->n);
// 	for (i = 0; i < mv->n; ++i)
// 		a[i].x = mv->a[i].x, a[i].y = i;
// 	radix_sort_128x(a, a + mv->n);
// 	for (st = 0, i = 1; i <= mv->n; ++i) {
// 		if (i == mv->n || a[i].x != a[st].x) {
// 			int32_t cnt = i - st;
// 			if (cnt > q_occ_max && cnt > mv->n * q_occ_frac)
// 				for (j = st; j < i; ++j)
// 					mv->a[a[j].y].x = 0;
// 			st = i;
// 		}
// 	}
// 	ri_kfree(km, a);
// 	for (i = j = 0; i < mv->n; ++i)
// 		if (mv->a[i].x != 0)
// 			mv->a[j++] = mv->a[i];
// 	mv->n = j;
// }

std::string anchors_to_string(const ri_anchor_t *anchors, uint32_t n_anchors){
	std::stringstream ss;
	for(size_t i = 0; i < n_anchors; i++){
		ss << "(";
		ss << anchors[i].query_position;
		ss << ",";
		ss << anchors[i].target_position;
		ss << ")";
	}
	return ss.str();
}

void comp_mapq(std::vector<ri_chain_t> &chains, const ri_mapopt_t *opt) {
  if (chains.size() == 1) {
    chains[0].mapq = 60;
    return;
  } else {
    int mapq;
	if(opt->flag & RI_M_DTW_EVALUATE_CHAINS){
		//TODO: Possibly need to adapt this mapq calculation, since the alignment score might behave differently than the chaining score
		//currently, this is the same as the chaining score calculation from rawhash
		mapq = 40 * (1 - chains[1].alignment_score / chains[0].alignment_score);
	}
	else{
		mapq = 40 * (1 - chains[1].chaining_score / chains[0].chaining_score);
	}

	if (mapq > 60) {
      mapq = 60;
    }
    if (mapq < 0) {
      mapq = 0;
    }
    chains[0].mapq = (uint8_t)mapq;
  }
}

void gen_primary_chains(void *km, std::vector<ri_chain_t> &chains, const ri_mapopt_t *opt) {
	std::sort(chains.begin(), chains.end(), std::greater<ri_chain_t>());
	std::vector<ri_chain_t> primary_chains;
	primary_chains.reserve(chains.size());
	primary_chains.emplace_back(chains[0]);
	for (uint32_t ci = 1; ci < chains.size(); ++ci) {
		bool is_primary = true;
		//if the score is less than 1/3 of the best score, it's not primary
		//TODO: For some reason this break is not a memory leak, even though it skeps the ri_kfree below.
		//Why?
		if(opt->flag & RI_M_DTW_EVALUATE_CHAINS){
			if (chains[ci].alignment_score < primary_chains.back().alignment_score / 3) {
				break;
			}
		}
		else{
			if (chains[ci].chaining_score < primary_chains.back().chaining_score / 3) {
				break;
			}
		}

		//check for any overlap with previous primary chains
		//if it's an overlap, it's not primary
		for (uint32_t pi = 0; pi < primary_chains.size(); ++pi) {
			if (chains[ci].reference_sequence_index == primary_chains[pi].reference_sequence_index) {
				if (std::max(chains[ci].start_position, primary_chains[pi].start_position) <= std::min(chains[ci].end_position, primary_chains[pi].end_position)) {
					is_primary = false;
					break;
				}
			}
		}
		if (is_primary) {
			primary_chains.emplace_back(chains[ci]);
		}else{
			if(chains[ci].anchors)ri_kfree(km, chains[ci].anchors);
		}
	}
	chains.swap(primary_chains);
}

void traceback_chains(void *km,
    int min_num_anchors, int strand, size_t chain_end_anchor_index,
    uint32_t chain_target_signal_index,
    const std::vector<float> &chaining_scores,
    const std::vector<size_t> &chaining_predecessors,
    const std::vector<std::vector<ri_anchor_t> > &anchors_fr,
    std::vector<bool> &anchor_is_used, std::vector<ri_chain_t> &chains) {

  if (!anchor_is_used[chain_end_anchor_index]) {
    std::vector<ri_anchor_t> anchors;
    anchors.reserve(100);
    bool stop_at_an_used_anchor = false;
    size_t chain_start_anchor_index = chain_end_anchor_index;
    // Add the end anchor
    anchors.push_back(anchors_fr[chain_target_signal_index][chain_start_anchor_index]);
    // The end anchor has an used predecessor
    if (chaining_predecessors[chain_start_anchor_index] != chain_start_anchor_index && anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
      stop_at_an_used_anchor = true;
    }
    anchor_is_used[chain_start_anchor_index] = true;
    uint32_t chain_num_anchors = 1;
    while (chaining_predecessors[chain_start_anchor_index] != chain_start_anchor_index && !anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
      chain_start_anchor_index = chaining_predecessors[chain_start_anchor_index];
      anchors.push_back(anchors_fr[chain_target_signal_index][chain_start_anchor_index]);
      if (chaining_predecessors[chain_start_anchor_index] != chain_start_anchor_index && anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
        stop_at_an_used_anchor = true;
      }
      anchor_is_used[chain_start_anchor_index] = true;
      ++chain_num_anchors;
    }
    if (chain_num_anchors >= (uint32_t)min_num_anchors) {
      float adjusted_chaining_score = chaining_scores[chain_end_anchor_index];
      if (stop_at_an_used_anchor) {
        adjusted_chaining_score -= chaining_scores[chaining_predecessors[chain_start_anchor_index]];
      }
	  ri_anchor_t* a_anchors = (ri_anchor_t*)ri_kmalloc(km,anchors.size()*sizeof(ri_anchor_t));
	  std::copy(anchors.begin(), anchors.end(), a_anchors);
      chains.emplace_back(ri_chain_t{adjusted_chaining_score, 0.0f, chain_target_signal_index, 
		  				  anchors_fr[chain_target_signal_index][chain_start_anchor_index].target_position,
                          anchors_fr[chain_target_signal_index][chain_end_anchor_index].target_position,
                          chain_num_anchors, 0, strand, a_anchors});
    }
  }
}

bool compare(const std::pair<float, size_t> &left, const std::pair<float, size_t> &right) {
	if (left.first > right.first) return true;
	else if (left.first == right.first) return (left.second > right.second);
	else return false;
}

void align_chain(ri_chain_t &chain, const ri_idx_t *ri, const float* read_events, const uint32_t n_read_events, const uint32_t chunk_start, const ri_mapopt_t *opt, bool cigar=false, float min_score=-1e10){
	float *ref_events;
	if(chain.strand == 1){
		ref_events = ri->forward_signals[chain.reference_sequence_index];
	}
	else{
		ref_events = ri->reverse_signals[chain.reference_sequence_index];
	}

	float dtw_cost = 0.0f;
	uint32_t num_aligned_read_events = 0;
	if(opt->dtw_border_constraint == RI_M_DTW_BORDER_CONSTRAINT_GLOBAL){
		//note that these indices are not flipped
		//end_anchor is anchors[0] due to chaining's traceback from right to left
		const ri_anchor_t &start_anchor = chain.anchors[chain.n_anchors-1];
		const ri_anchor_t &end_anchor = chain.anchors[0];

		float *ref_region = ref_events + start_anchor.target_position;
		uint32_t ref_region_size = end_anchor.target_position - start_anchor.target_position + 1;

		const float *read_region = read_events + start_anchor.query_position;
		uint32_t read_region_size = end_anchor.query_position - start_anchor.query_position + 1;

		if(!cigar){
			float max_attainable_score = read_region_size*opt->dtw_match_bonus;
			if(max_attainable_score < min_score){
				chain.alignment_score = -1e10;
				return;
			}
			if(opt->dtw_fill_method == RI_M_DTW_FILL_METHOD_FULL){
				dtw_cost = DTW_global(read_region, read_region_size, ref_region, ref_region_size);
			}
			else{
				int band_radius = std::max(1, (int)(read_region_size*opt->dtw_band_radius_frac));
				dtw_cost = DTW_global_slantedbanded_antidiagonalwise(read_region, read_region_size, ref_region, ref_region_size, band_radius);
			}
		}
		else{
			dtw_result res;
			if(opt->dtw_fill_method == RI_M_DTW_FILL_METHOD_FULL){
				res = DTW_global_tb(read_region, read_region_size, ref_region, ref_region_size);
			}
			else{
				assert(false); //not implemented
			}
			dtw_cost = res.cost;
			new(&chain.dtw_result) dtw_result; //ensure that c++ objects (currently, a vector) inside chain are properly initialized
												//this might cause a memory leak if we overwrite the dtw_result multiple times like this
												//not sure if there's a nice fix
			for(size_t i=0; i<res.alignment.size(); i++){
				res.alignment.back().position.i += start_anchor.query_position;
				res.alignment.back().position.j += start_anchor.target_position;
			}
			chain.dtw_result = res;
		}
		num_aligned_read_events = read_region_size;
	}
	else if(opt->dtw_border_constraint == RI_M_DTW_BORDER_CONSTRAINT_SPARSE){
		uint32_t alignment_parts = chain.n_anchors-1;
		std::vector<alignment_element> alignment;

		//these are only needed for the early termination condition
		const ri_anchor_t &chain_start_anchor = chain.anchors[chain.n_anchors-1];
		const ri_anchor_t &chain_end_anchor = chain.anchors[0];
		uint32_t chain_read_region_size = chain_end_anchor.query_position - chain_start_anchor.query_position + 1;
		float current_max_attainable_score = chain_read_region_size*opt->dtw_match_bonus; //max score if the remaining part of the read were matching perfectly

		for(size_t alignment_part=0; alignment_part<alignment_parts; alignment_part++){
			//TODO: some of these could be merged when they are very short, for performance

			//note that these indices are not flipped
			//due to chaining's traceback from right to left
			const ri_anchor_t &start_anchor = chain.anchors[alignment_parts-alignment_part];
			const ri_anchor_t &end_anchor = chain.anchors[alignment_parts-alignment_part-1];

			float *ref_region = ref_events + start_anchor.target_position;
			uint32_t ref_region_size = end_anchor.target_position - start_anchor.target_position + 1;

			const float *read_region = read_events + start_anchor.query_position;
			uint32_t read_region_size = end_anchor.query_position - start_anchor.query_position + 1;
			
			if(!cigar){
				//exclude last element for all but the last alignment part
				//to avoid double-counting the last element
				if(current_max_attainable_score < min_score){
					chain.alignment_score = -1e10;
					return;
				}

				bool exclude_last_element = (alignment_part != alignment_parts-1);
				float sub_dtw_cost;
				if(opt->dtw_fill_method == RI_M_DTW_FILL_METHOD_FULL){
					sub_dtw_cost = DTW_global(read_region, read_region_size, ref_region, ref_region_size, exclude_last_element);
				}
				else{
					int band_radius = std::max(1, (int)(read_region_size*opt->dtw_band_radius_frac));
					sub_dtw_cost = DTW_global_slantedbanded_antidiagonalwise(read_region, read_region_size, ref_region, ref_region_size, band_radius, exclude_last_element);
				}
				dtw_cost += sub_dtw_cost;
				current_max_attainable_score -= sub_dtw_cost;
			}
			else{
				bool exclude_last_element = (alignment_part != alignment_parts-1);
				dtw_result sub_res = DTW_global_tb(read_region, read_region_size, ref_region, ref_region_size);
				for(size_t i=0; i<sub_res.alignment.size(); i++){
					alignment.push_back(sub_res.alignment[i]);
					alignment.back().position.i += start_anchor.query_position;
					alignment.back().position.j += start_anchor.target_position;
				}
				dtw_cost+=sub_res.cost;
			}
			num_aligned_read_events += read_region_size;
		}
		if(cigar){
			new(&chain.dtw_result) dtw_result; //ensure that c++ objects (currently, a vector) inside chain are properly initialized
												//this might cause a memory leak if we overwrite the dtw_result multiple times like this
												//not sure if there's a nice fix
			chain.dtw_result = dtw_result{dtw_cost, alignment};
		}
	}
	else{
		fprintf(stderr, "ERROR: invalid border constraint\n");
		exit(EXIT_FAILURE);
	}

	chain.alignment_score = num_aligned_read_events*opt->dtw_match_bonus - dtw_cost;

	if(opt->flag & RI_M_DTW_LOG_SCORES){
		char str[256];
		sprintf(str, "chaining_score=%f alignment_score=%f\n", chain.chaining_score, chain.alignment_score);
		fprintf(stderr, str);
	}
}

void gen_chains(void *km, const pipeline_mt *p, const ri_idx_t *ri, const float* chunk_events, const uint32_t l_chunk_events, const uint32_t chunk_start, const size_t n_seq, ri_reg1_t* reg, const ri_mapopt_t *opt){

	// Chaining parameters
	int max_gap_length = opt->max_gap_length;
	int max_target_gap_length = opt->max_target_gap_length;
	int chaining_band_length = opt->chaining_band_length;
	int max_num_skips = opt->max_num_skips;
	int min_num_anchors = opt->min_num_anchors;
	int num_best_chains = opt->num_best_chains;
	float min_chaining_score = opt->min_chaining_score;

	uint32_t mask = (1ULL<<31)-1;

	// kvec_t(ri_anchor_t)* anchors_fr[2] = {0};
	// kvec_t(ri_anchor_t)* anchors_f = anchors_fr;
	// kvec_t(ri_anchor_t)* anchors_r = anchors_fr+1;

	// anchors_f = (kvec_t(ri_anchor_t)*)calloc(n_seq, sizeof(kvec_t(ri_anchor_t)));
	// anchors_r = (kvec_t(ri_anchor_t)*)calloc(n_seq, sizeof(kvec_t(ri_anchor_t)));
	// memset(anchors_f, 0, n_seq*sizeof(kvec_t(ri_anchor_t)));
	// memset(anchors_r, 0, n_seq*sizeof(kvec_t(ri_anchor_t)));

	std::vector<std::vector<std::vector<ri_anchor_t> > > anchors_fr(2);
	anchors_fr[0] = std::vector<std::vector<ri_anchor_t> >(n_seq);
	anchors_fr[1] = std::vector<std::vector<ri_anchor_t> >(n_seq);
	std::vector<std::vector<ri_anchor_t> > &anchors_f = anchors_fr[0];
	std::vector<std::vector<ri_anchor_t> > &anchors_r = anchors_fr[1];

	// Get anchors in previous chains
	ri_chain_t* previous_chains = reg->chains;
	if (reg->n_chains > 0) {
		for (uint32_t c_ind = 0; c_ind < reg->n_chains; ++c_ind) {
			int strand = previous_chains[c_ind].strand;
			uint32_t reference_sequence_index = previous_chains[c_ind].reference_sequence_index;
			// anchors_fr[strand][reference_sequence_index].reserve(previous_chains[c_ind].n_anchors);
			anchors_fr[strand][reference_sequence_index].reserve(previous_chains[c_ind].n_anchors);
			// kv_resize(ri_anchor_t, 0, anchors_fr[strand][reference_sequence_index], anchors_fr[strand][reference_sequence_index].n + previous_chains[c_ind].n_anchors)
			for (uint32_t a_ind = 0; a_ind < previous_chains[c_ind].n_anchors; ++a_ind) {
				anchors_fr[strand][reference_sequence_index].emplace_back(previous_chains[c_ind].anchors[a_ind]);
				// kv_push(ri_anchor_t, 0, anchors_fr[strand][reference_sequence_index], previous_chains[c_ind].anchors[a_ind]);
			}
		}
	}

	if(reg->chains){
		for(uint32_t i = 0; i < reg->n_chains; ++i) if(reg->chains[i].anchors) ri_kfree(km, reg->chains[i].anchors);
		ri_kfree(km, reg->chains); reg->n_chains = 0; reg->chains=NULL;
	}

	uint32_t i, pi;
	mm128_v riv = {0,0,0};
	uint64_t hashVal, keyval;
	ri_sketch(km, chunk_events, 0, 0, l_chunk_events, ri->w, ri->e, ri->n, ri->q, ri->lq, ri->k, &riv);
	//   ri_seed_mz_flt(0, &riv, 1000, 0.01f);

	uint64_t chunk_seed_hits = 0;
	for (i = 0; i < riv.n; ++i) {
		hashVal = riv.a[i].x>>RI_HASH_SHIFT;
		pi = (uint32_t)riv.a[i].y>>RI_POS_SHIFT;
		if(pi >= l_chunk_events){
			fprintf(stderr, "i=%d pi=%d l_sig=%d\n", i, pi, l_chunk_events);
		}

		const uint64_t *cr;
		int t;
		cr = ri_idx_get(ri, hashVal, &t);
		chunk_seed_hits += t;

		if(t == 0) continue;

		for(int s = 0; s < t; ++s){
			keyval = cr[s];
			uint32_t t_ind = (uint32_t)(keyval>>RI_ID_SHIFT), target_signal_position = (uint32_t)(keyval>>RI_POS_SHIFT)&mask;
			anchors_fr[(keyval&1)][t_ind].emplace_back(ri_anchor_t{target_signal_position,  pi + chunk_start});
			// kv_push(ri_anchor_t, 0, anchors_fr[(keyval&1)][t_ind], ri_anchor_t{target_signal_position, pi + q_offset})
		}
	}

	ri_kfree(km, riv.a);

	// Sort the anchors based on their occurrence on target signal
	for (size_t t_ind = 0; t_ind < n_seq; ++t_ind) {
		// std::sort(anchors_f[t_ind].a, anchors_f[t_ind].a + anchors_f[t_ind].n);
		// std::sort(anchors_r[t_ind].a, anchors_r[t_ind].a + anchors_r[t_ind].n);
		std::sort(anchors_f[t_ind].begin(), anchors_f[t_ind].end());
		std::sort(anchors_r[t_ind].begin(), anchors_r[t_ind].end());
	} 

	if(opt->flag & RI_M_LOG_ANCHORS){
		for(size_t t_ind = 0; t_ind < n_seq; ++t_ind){
			for(int strand = 0; strand <= 1; strand++){
				std::string anchors = anchors_to_string(&anchors_fr[strand][t_ind][0], anchors_fr[strand][t_ind].size());
				std::stringstream anchor_log_stream;
				anchor_log_stream << "readname=" << reg->read_name;
				anchor_log_stream << " refname=" << ri->seq[t_ind].name;
				anchor_log_stream << " strand=" << strand;
				anchor_log_stream << " anchors=" << anchors;
				anchor_log_stream << "\n";
				std::cerr << anchor_log_stream.str();
			}
		}
	}

	if(opt->flag & RI_M_LOG_NUM_ANCHORS){
		std::stringstream num_anchor_log_stream;
		num_anchor_log_stream << "readname=" << reg->read_name;
		num_anchor_log_stream << " pos=[" << chunk_start << "," << chunk_start+l_chunk_events-1 << "]"; 
		num_anchor_log_stream << " num_anchors=" << chunk_seed_hits;
		num_anchor_log_stream << "\n";
		std::cerr << num_anchor_log_stream.str();
	}

	// Chaining DP done on each individual  target signal
	float max_chaining_score = 0;  // std::numeric_limits<float>::min();
	std::vector<ri_chain_t> chains;
	for (size_t target_signal_index = 0; target_signal_index < n_seq; ++target_signal_index) {
		for (int strand_i = 0; strand_i < 2; ++strand_i) {
			std::vector<float> chaining_scores;
			chaining_scores.reserve(anchors_fr[strand_i][target_signal_index].size());
			std::vector<size_t> chaining_predecessors;
			chaining_predecessors.reserve(anchors_fr[strand_i][target_signal_index].size());
			std::vector<bool> anchor_is_used;
			anchor_is_used.reserve(anchors_fr[strand_i][target_signal_index].size());
			std::vector<std::pair<float, size_t> > end_anchor_index_chaining_scores;
			end_anchor_index_chaining_scores.reserve(10);
			for (size_t anchor_index = 0; anchor_index < anchors_fr[strand_i][target_signal_index].size(); ++anchor_index) {
				float distance_coefficient = 1; /*- 0.2 * anchors_fr[strand_i][target_signal_index][anchor_index].distance / search_radius;*/
				chaining_scores.emplace_back(distance_coefficient * ri->e);
				chaining_predecessors.emplace_back(anchor_index);
				anchor_is_used.push_back(false);
				int32_t current_anchor_target_position = anchors_fr[strand_i][target_signal_index][anchor_index].target_position;
				int32_t current_anchor_query_position = anchors_fr[strand_i][target_signal_index][anchor_index].query_position;
				int32_t start_anchor_index = 0;
				if (anchor_index > (size_t)chaining_band_length) start_anchor_index = anchor_index - chaining_band_length;

				int32_t previous_anchor_index = anchor_index - 1;
				int32_t num_skips = 0;
				for (; previous_anchor_index >= start_anchor_index; --previous_anchor_index) {
					int32_t previous_anchor_target_position = anchors_fr[strand_i][target_signal_index][previous_anchor_index].target_position;
					int32_t previous_anchor_query_position = anchors_fr[strand_i][target_signal_index][previous_anchor_index].query_position;

					if (previous_anchor_query_position == current_anchor_query_position) continue;
					if (previous_anchor_target_position == current_anchor_target_position) continue;
					if (previous_anchor_target_position + max_target_gap_length < current_anchor_target_position) break;

					int32_t target_position_diff = current_anchor_target_position - previous_anchor_target_position;
					assert(target_position_diff > 0);
					int32_t query_position_diff = current_anchor_query_position - previous_anchor_query_position;
					float current_chaining_score = 0;

					if (query_position_diff < 0) continue;
					
					float matching_dimensions = std::min(std::min(target_position_diff, query_position_diff), ri->e) * distance_coefficient;
					int gap_length = std::abs(target_position_diff - query_position_diff);
					float gap_scale = target_position_diff > 0? (float)query_position_diff / target_position_diff: 1;
					// float gap_cost = 0;
					// if (gap_length != 0) {
					if (gap_length < max_gap_length && gap_scale < 5 && gap_scale > 0.75) {
						current_chaining_score = chaining_scores[previous_anchor_index] + matching_dimensions;  // - gap_cost;
					}

					if (current_chaining_score > chaining_scores[anchor_index]) {
						chaining_scores[anchor_index] = current_chaining_score;
						chaining_predecessors[anchor_index] = previous_anchor_index;
						--num_skips;
					} else {
						++num_skips;
						if (num_skips > max_num_skips) break;
					}
				}
				// Update chain with max score
				if (chaining_scores[anchor_index] > max_chaining_score) {
					max_chaining_score = chaining_scores[anchor_index];
				}
				if (opt->flag & RI_M_DISABLE_CHAININGSCORE_FILTERING ||
				    chaining_scores.back() >= min_chaining_score &&
					chaining_scores.back() > max_chaining_score / 2) {
					end_anchor_index_chaining_scores.emplace_back(chaining_scores.back(), anchor_index);
				}
			}
			// Sort a vector of <end anchor index, chaining score>
			std::sort(end_anchor_index_chaining_scores.begin(), end_anchor_index_chaining_scores.end(), compare);
			// Traceback all reg->chains from higest score to lowest
			for (size_t anchor_index = 0; anchor_index < end_anchor_index_chaining_scores.size() && anchor_index < (size_t)num_best_chains; ++anchor_index) {
				traceback_chains(km, min_num_anchors, strand_i, end_anchor_index_chaining_scores[anchor_index].second, target_signal_index, chaining_scores, 
								chaining_predecessors, anchors_fr[strand_i], anchor_is_used, chains);
	
				if(!(opt->flag & RI_M_DISABLE_CHAININGSCORE_FILTERING)){
					if (chaining_scores[end_anchor_index_chaining_scores[anchor_index].second] < max_chaining_score / 2) break;
				}
			}
		}
	}

	if(opt->flag & RI_M_DTW_EVALUATE_CHAINS || opt->flag & RI_M_DTW_LOG_SCORES){
		//the sole purpose of this sort is to ideally evaluate the best chain first, so worse chains can be eliminated quicker
		//maybe a different sorting method yields better performance (e.g., according to length, number of anchors, or a combination thereof)
		std::sort(chains.begin(), chains.end(), [](ri_chain_t &a, ri_chain_t &b){return a.chaining_score > b.chaining_score;});
		
		std::vector<ri_chain_t> post_alignment_chains;
		float best_found_alignment = 0.0f; //this could be slightly more agressive and be set to opt->dtw_min_score immediately, but starting with 0 if clearer
		for(ri_chain_t &chain: chains){
			align_chain(chain, ri, p->events[reg->read_id].values, chunk_start + l_chunk_events, chunk_start, opt, false, best_found_alignment);
			if(chain.alignment_score >= opt->dtw_min_score){
				if(chain.alignment_score > best_found_alignment){
					best_found_alignment = chain.alignment_score;
				}
				post_alignment_chains.push_back(chain);
			}
		}
		if(opt->flag & RI_M_DTW_EVALUATE_CHAINS){
			chains = post_alignment_chains;
			//sorting again is likely unnecessary due to gen_primary_chains sorting immediately after this
			//std::sort(chains.begin(), chains.end(), [](ri_chain_t &a, ri_chain_t &b){return a.alignment_score > b.alignment_score;});
		}
	}

	if (chains.size() > 0) {
		// Generate primary reg->chains
		gen_primary_chains(km, chains, opt);
		// Compute MAPQ
		comp_mapq(chains, opt);

		reg->n_chains = chains.size();
		reg->chains = (ri_chain_t*)ri_kmalloc(km, chains.size()*sizeof(ri_chain_t));
		std::copy(chains.begin(), chains.end(), reg->chains);
	}
}

//returns n_regs // currently we report one mapping
void ri_map_frag(pipeline_mt *p, const ri_idx_t *ri, const uint32_t s_len, const float *sig, ri_reg1_t* reg, ri_tbuf_t *b, const ri_mapopt_t *opt, const char *qname){
	
	uint32_t n_chunk_events = 0;
	// ri_reg1_t* reg = 0;
	// uint32_t event_off = 0;

	float* chunk_events = detect_events(b->km, s_len, sig, opt, &n_chunk_events);

	//update the globally accumulated events of the read by concatenating the new chunk to earleir events
	ri_events_t &global_events = p->events[reg->read_id];
	if(!global_events.values){
		global_events.values = (float *)malloc(n_chunk_events*sizeof(float));
		global_events.length = n_chunk_events;
		std::copy(chunk_events, chunk_events+n_chunk_events, global_events.values);
	}
	else{
		float* new_events = (float*)malloc((global_events.length+n_chunk_events)*sizeof(float));
		std::copy(global_events.values, global_events.values+global_events.length, new_events);
		std::copy(chunk_events, chunk_events+n_chunk_events, new_events+global_events.length);
		free(global_events.values);
		global_events.values = new_events;
		global_events.length += n_chunk_events;
	}

	if(n_chunk_events < opt->min_events) {
		if(chunk_events)ri_kfree(b->km, chunk_events);
		return;
	}

	gen_chains(b->km, p, ri, chunk_events, n_chunk_events, reg->offset, ri->n_seq, reg, opt);
	reg->offset += n_chunk_events;

	if(chunk_events)ri_kfree(b->km, chunk_events);
}

std::string dtwresult_to_string(const dtw_result& res){
	std::stringstream ss;
	for(size_t i = 0; i < res.alignment.size(); i++){
		ss << "(";
		ss << res.alignment[i].position.i;
		ss << ",";
		ss << res.alignment[i].position.j;
		ss << ",";
		ss << res.alignment[i].difference;
		ss << ")";
	}
	return ss.str();
}

bool is_mapped_with_high_confidence(const ri_reg1_t* reg0, const ri_mapopt_t *opt, bool log_reason=false){
	//if there are no chains, the read is unmapped
	//unsure if this is explicitly necessary, but copied from where rawhash did this check previously
	uint32_t n_anchors0 = (reg0->n_chains)?reg0->chains[0].n_anchors:0;
	if(n_anchors0==0) return false;
	// Since alignment scores and chaining scores are not on the same scale,
	// the following conditions might need to be specialized for each case.
	// Currently we use the conditions from the original rawhash for both cases.
	if(opt->flag & RI_M_DTW_EVALUATE_CHAINS){
		if (reg0->n_chains >= 2) {
			if (reg0->chains[0].alignment_score / reg0->chains[1].alignment_score >= opt->min_bestmap_ratio){
				if(log_reason)
					std::cerr << "mapping_reason=bestmap_ratio\n";
				return true;
			}

			float mean_alignment_score = 0;
			for (uint32_t c_ind = 0; c_ind < reg0->n_chains; ++c_ind)
				mean_alignment_score += reg0->chains[c_ind].alignment_score;

			mean_alignment_score /= reg0->n_chains;
			if (reg0->chains[0].alignment_score >= opt->min_meanmap_ratio * mean_alignment_score){
				if(log_reason)
					std::cerr << "mapping_reason=meanmap_ratio\n";
				return true;
			}
		} else if (reg0->n_chains == 1 && reg0->chains[0].n_anchors >= opt->min_chain_anchor){
			if(log_reason)
				std::cerr << "mapping_reason=unique\n";
			return true;
		}
		else if(reg0->n_chains > 1){
			if(log_reason)
				std::cerr << "non_mapping_reason=ambiguous\n";
			return false;
		}
		else if(reg0->n_chains == 1 && reg0->chains[0].n_anchors < opt->min_chain_anchor){
			if(log_reason){
				std::stringstream ss;
				ss << "non_mapping_reason=short_chain (" << reg0->chains[0].n_anchors << "<" << opt->min_chain_anchor << ")\n";
				std::cerr << ss.str();
			}
			return false;
		}
		else{
			if(log_reason)
				std::cerr << "non_mapping_reason=no_chains\n";
			return false;
		}
	}
	else{
		//original rawhash never outputs multiple chains
		//if case there are multiple chains remaining
		//- terminate if the best one (index 0) is significantly (i.e., at least opt->min_bestmap_ratio) better than the second best one (index 1)
		//- and terminate if the best one (index 0) is significantly (i.e., at least opt->min_meanmap_ratio) better than the mean score of all chains
		//if there is only a single chain remaining
		//- terminate if the chain has at least opt->min_chain_anchor anchors
		if (reg0->n_chains >= 2) {
			if (reg0->chains[0].chaining_score / reg0->chains[1].chaining_score >= opt->min_bestmap_ratio) return true;

			float mean_chain_score = 0;
			for (uint32_t c_ind = 0; c_ind < reg0->n_chains; ++c_ind)
				mean_chain_score += reg0->chains[c_ind].chaining_score;

			mean_chain_score /= reg0->n_chains;
			if (reg0->chains[0].chaining_score >= opt->min_meanmap_ratio * mean_chain_score) return true;
		} else if (reg0->n_chains == 1 && reg0->chains[0].n_anchors >= opt->min_chain_anchor){
			return true;
		}
	}
	return false;
}

static void map_worker_for(void *_data, long i, int tid) // kt_for() callback 
{
    step_mt *s = (step_mt*)_data; //s->sig and s->n_sig (signals read in this step and num of them)
	const ri_mapopt_t *opt = s->p->opt;
	double t = 0.0;
	ri_tbuf_t* b = s->buf[tid];
	ri_reg1_t* reg0 = s->reg[i];
	ri_sig_t* sig = s->sig[i];

	uint32_t l_chunk = opt->chunk_size;
	uint32_t max_chunk =  opt->max_num_chunk;
	uint32_t qlen = sig->l_sig;
	uint32_t chunk_start; //inclusve
	uint32_t chunk_end; //exclusive

	uint32_t current_chunk = 0;

	t = ri_realtime();
	for (chunk_start = current_chunk = 0; chunk_start < qlen && current_chunk < max_chunk; chunk_start += l_chunk, ++current_chunk) {
		chunk_end = chunk_start + l_chunk;
		if(chunk_end > qlen) chunk_end = qlen;

		ri_map_frag(s->p, s->p->ri, (const uint32_t)chunk_end-chunk_start, (const float*)&(sig->sig[chunk_start]), reg0, b, opt, sig->name);

		//early termination conditions for the mapping of a read
		if(is_mapped_with_high_confidence(reg0, opt)) break;
	}
	double mapping_time = ri_realtime() - t;

	if (current_chunk > 0 && (chunk_start >= qlen || current_chunk == max_chunk)) --current_chunk;

	float read_position_scale = ((float)(current_chunk+1)*l_chunk/reg0->offset) / ((float)opt->sample_rate/opt->bp_per_sec);
	// Save results in vector and output PAF
	
	ri_chain_s* chains = reg0->chains;

	if(!chains) {reg0->n_chains = 0;}

	uint32_t n_anchors0 = (reg0->n_chains)?chains[0].n_anchors:0;
	float mean_chain_score = 0;
	for (uint32_t c_ind = 0; c_ind < reg0->n_chains; ++c_ind)
		mean_chain_score += chains[c_ind].chaining_score;

	mean_chain_score /= reg0->n_chains;
	if (is_mapped_with_high_confidence(reg0, opt)) {

		//maybe we should align based on reg0->read_start_position and reg0->read_end_position
		//calculated below
		if(opt->flag & RI_M_DTW_OUTPUT_CIGAR){
			align_chain(chains[0], s->p->ri, s->p->events[reg0->read_id].values, qlen, chunk_start, opt, true);
		}

		float anchor_ref_gap_avg_length = 0;
		float anchor_read_gap_avg_length = 0;
		for (size_t ai = 0; ai < n_anchors0; ++ai) {
			if (ai < n_anchors0 - 1) {
				anchor_ref_gap_avg_length += chains[0].anchors[ai].target_position - chains[0].anchors[ai + 1].target_position;
				anchor_read_gap_avg_length += chains[0].anchors[ai].query_position - chains[0].anchors[ai + 1].query_position;
			}
		}

		anchor_ref_gap_avg_length /= n_anchors0;
		anchor_read_gap_avg_length /= n_anchors0;
		std::string tags;
		tags.append("mt:f:" + std::to_string(mapping_time * 1000));
		tags.append("\tci:i:" + std::to_string(current_chunk + 1));
		tags.append("\tsl:i:" + std::to_string(qlen));
		tags.append("\tcm:i:" + std::to_string(n_anchors0));
		tags.append("\tnc:i:" + std::to_string(reg0->n_chains));
		tags.append("\ts1:f:" + std::to_string(chains[0].chaining_score));
		tags.append("\ts2:f:" + std::to_string(reg0->n_chains > 1 ? chains[1].chaining_score : 0));
		tags.append("\tsm:f:" + std::to_string(mean_chain_score));
		tags.append("\tat:f:" + std::to_string((float)anchor_ref_gap_avg_length));
		tags.append("\taq:f:" + std::to_string((float)anchor_read_gap_avg_length));
		if(opt->flag & RI_M_DTW_OUTPUT_CIGAR){ //output only the alignment for chain 0, RawHash only outputs one mapping
			tags.append("\talns:f:" + std::to_string(chains[0].alignment_score));
			tags.append("\taln:s:" + dtwresult_to_string(chains[0].dtw_result));
		}
		if(opt->flag & RI_M_OUTPUT_CHAINS){
			tags.append("\tanchors:s:" + anchors_to_string(chains[0].anchors, chains[0].n_anchors));
		}

		reg0->ref_id = chains[0].reference_sequence_index;
		reg0->read_length = (uint32_t)(read_position_scale * chains[0].anchors[0].query_position);
		reg0->read_start_position = (uint32_t)(read_position_scale*chains[0].anchors[n_anchors0-1].query_position);
		reg0->read_end_position = (uint32_t)(read_position_scale * chains[0].anchors[0].query_position);
		reg0->fragment_start_position = chains[0].strand?(uint32_t)(s->p->ri->seq[chains[0].reference_sequence_index].len+1-chains[0].end_position):chains[0].start_position;
		reg0->fragment_length = (uint32_t)(chains[0].end_position - chains[0].start_position + 1);
		reg0->mapq = chains[0].mapq;
		reg0->rev = chains[0].strand;
		reg0->mapped = 1; 
		reg0->tags = strdup(tags.c_str());
	} else {
		std::string tags;
		tags.append("mt:f:" + std::to_string(mapping_time * 1000));
		tags.append("\tci:i:" + std::to_string(current_chunk + 1));
		tags.append("\tsl:i:" + std::to_string(qlen));
		if (reg0->n_chains >= 1) {
			float anchor_ref_gap_avg_length = 0;
			float anchor_read_gap_avg_length = 0;
			for (size_t ai = 0; ai < n_anchors0; ++ai) {
				if (ai < n_anchors0 - 1) {
					anchor_ref_gap_avg_length += chains[0].anchors[ai].target_position - chains[0].anchors[ai + 1].target_position;
					anchor_read_gap_avg_length += chains[0].anchors[ai].query_position - chains[0].anchors[ai + 1].query_position;
				}
			}
			if(n_anchors0)anchor_ref_gap_avg_length /= n_anchors0;
			if(n_anchors0)anchor_read_gap_avg_length /= n_anchors0;
			tags.append("\tcm:i:" + std::to_string(n_anchors0));
			tags.append("\tnc:i:" + std::to_string(reg0->n_chains));
			tags.append("\ts1:f:" + std::to_string(chains[0].chaining_score));
			tags.append("\ts2:f:" + std::to_string(reg0->n_chains > 1 ? chains[1].chaining_score: 0));
			tags.append("\tsm:f:" + std::to_string(mean_chain_score));
			tags.append("\tat:f:" + std::to_string((float)anchor_ref_gap_avg_length));
			tags.append("\taq:f:" + std::to_string((float)anchor_read_gap_avg_length));
		}else {
			tags.append("\tcm:i:0");
			tags.append("\tnc:i:0");
			tags.append("\ts1:f:0");
			tags.append("\ts2:f:0");
			tags.append("\tsm:f:0");
			tags.append("\tat:f:0");
			tags.append("\taq:f:0");
		}

		reg0->ref_id = 0;
		reg0->read_length = (uint32_t)(read_position_scale * reg0->offset);
		reg0->read_start_position = 0;
		reg0->read_end_position = 0;
		reg0->fragment_start_position = 0;
		reg0->fragment_length = 0;
		reg0->mapq = 0;
		reg0->rev = 0;
		reg0->mapped = 0;
		reg0->tags = strdup(tags.c_str());
	}

	for (uint32_t c_ind = 0; c_ind < reg0->n_chains; ++c_ind){
			if(reg0->chains[c_ind].anchors) ri_kfree(b->km, reg0->chains[c_ind].anchors);
	}
	//call destructor on dtw_result of the first chain of a mapped read
	//(otherwise the dtw_result was not initialized)
	//Yes, I know it's ugly. Cost of miximg C and C++.
	if(reg0->mapped && opt->flag & RI_M_DTW_OUTPUT_CIGAR){
		reg0->chains[0].dtw_result.~dtw_result();
	}
	ri_kfree(b->km, reg0->chains);

	if (b->km) {
		ri_km_stat_t kmst;
		ri_km_stat(b->km, &kmst);
		// assert(kmst.n_blocks == kmst.n_cores);
		ri_km_destroy(b->km);
		b->km = ri_km_init();
	}
}

ri_sig_t** ri_sig_read_frag(pipeline_mt *pl, int64_t chunk_size, int *n_){

	int64_t size = 0;
	// int n_eof = 0; //at least one file with not end of file
	// kvec_t(ri_sig_t) a = {0,0,0};
	std::vector<ri_sig_t*> sigvec;
	*n_ = 0;
	if (pl->n_fp < 1) return 0;
	// if (a.m == 0) kv_resize(ri_sig_t, 0, a, 256);
	while (pl->fp) {
		//Reading data in bulk if the buffer is emptied
		while(pl->fp && pl->fp->cur_read == pl->fp->num_read){
			ri_sig_close(pl->fp);
			if(pl->cur_f < pl->n_f){
				if((pl->fp = open_sig(pl->f[pl->cur_f++])) == 0) break;
			}else if(pl->cur_fp < pl->n_fp){
				if(pl->f){
					for(int i = 0; i < pl->n_f; ++i) if(pl->f[i])free(pl->f[i]);
					free(pl->f);
				}
				pl->n_f = 0; pl->cur_f = 0;

				ri_char_v fnames = {0,0,0};
				kv_resize(char*, 0, fnames, 256);
				find_fast5(pl->fn[pl->cur_fp++], &fnames);
				pl->f =  fnames.a;
				if(!fnames.n || ((pl->fp = open_sig(pl->f[pl->cur_f++])) == 0)) break;
				pl->n_f = fnames.n;
				// ++n_read;
			}else {pl->fp = 0; break;}
		}

		if(!pl->fp || pl->fp->cur_read == pl->fp->num_read) break;
		
		ri_sig_t *s = (ri_sig_t*)calloc(1, sizeof(ri_sig_t));
		sigvec.push_back(s);
		ri_read_sig(pl->fp, s);
		size += s->l_sig;

		if(size >= chunk_size) break;
	}

	ri_sig_t** a = 0;

	if(sigvec.size()){
		a = (ri_sig_t**)calloc(sigvec.size(), sizeof(ri_sig_t**));
		std::copy(sigvec.begin(), sigvec.end(), a);
	}
	
	*n_ = sigvec.size();
	return a;
}

static void *map_worker_pipeline(void *shared, int step, void *in)
{
	int i, k;
    pipeline_mt *p = (pipeline_mt*)shared;
    if (step == 0) { // step 0: read sequences
        step_mt *s;
        s = (step_mt*)calloc(1, sizeof(step_mt));
		s->sig = ri_sig_read_frag(p, p->mini_batch_size, &s->n_sig);
		if (s->n_sig && !p->su_stop) {
			s->p = p;
			for (i = 0; i < s->n_sig; ++i){
				ri_events_t events;
				
				s->sig[i]->rid = p->n_processed++;
				p->signals.push_back(*s->sig[i]);

				events.rid = s->sig[i]->rid;
				events.name = s->sig[i]->name;
				events.values = NULL;
				events.length = 0;

				p->events.push_back(events);
			}
			s->buf = (ri_tbuf_t**)calloc(p->n_threads, sizeof(ri_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = ri_tbuf_init();
			s->reg = (ri_reg1_t**)calloc(s->n_sig, sizeof(ri_reg1_t*));
			for(i = 0; i < s->n_sig; ++i){
				s->reg[i] = (ri_reg1_t*)calloc(1, sizeof(ri_reg1_t));
				s->reg[i]->read_id = s->sig[i]->rid;
				s->reg[i]->read_name = s->sig[i]->name;
			}
			return s;
		} else if(s){
			if(s->sig) free(s->sig);
			free(s);
		}
    } else if (step == 1) { // step 1: detect events
		step_mt *s = (step_mt*)in;
		if(!p->su_stop) kt_for(p->n_threads, map_worker_for, in, s->n_sig);

		if(p->opt->flag & RI_M_SEQUENCEUNTIL && !p->su_stop){
			const ri_idx_t *ri = p->ri;
			for (k = 0; k < s->n_sig; ++k) {
				if(s->reg[k] && s->reg[k]->ref_id < ri->n_seq && s->reg[k]->read_name && s->reg[k]->mapped){
					ri_reg1_t* reg0 = s->reg[k];
					p->su_c_estimations[reg0->ref_id] += reg0->fragment_length;
					p->ab_count += reg0->fragment_length;
					p->su_nreads++;
					if(p->su_nreads > p->opt->tmin_reads && !(p->su_nreads%p->opt->ttest_freq)){
						//calculate abundance
						for(uint32_t ce = 0; ce < ri->n_seq; ++ce){
							p->su_estimations[p->su_cur][ce] =  (float)p->su_c_estimations[ce]/p->ab_count;
						}
						
						if(++p->su_cur >= p->opt->tn_samples) p->su_cur = 0;
						if(p->su_nestimations++ >= p->opt->tn_samples){
							if(find_outlier((const float**)p->su_estimations, ri->n_seq, p->opt->tn_samples) <= p->opt->t_threshold){
								//sending the stop signal.
								p->su_stop = k+1;
								fprintf(stderr, "[M::%s] Sequence Until is activated, stopping sequencing after processing %d mapped reads\n", __func__, p->su_nreads);
								break;
							}
						}
					}
				}
			}
		}
		return in;
    } else if (step == 2) { // step 2: align
		// void *km = 0;
        step_mt *s = (step_mt*)in;
		return in;
	} else if (step == 3) { // step 3: output
		// void *km = 0;
        step_mt *s = (step_mt*)in;
		const ri_idx_t *ri = p->ri;
		for (k = 0; k < s->n_sig; ++k) {
			if(s->reg[k]){
				ri_reg1_t* reg0 = s->reg[k];
				char strand = reg0->rev?'-':'+';
				
				if(reg0->ref_id < ri->n_seq && reg0->read_name){
					if(reg0->mapped && (!p->su_stop || k < p->su_stop) ){
						fprintf(stdout, "%s\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%u\t%u\t%d\t%s\n", 
										reg0->read_name, reg0->read_length, reg0->read_start_position, reg0->read_end_position, 
										strand, ri->seq[reg0->ref_id].name, ri->seq[reg0->ref_id].len, reg0->fragment_start_position, reg0->fragment_start_position + reg0->fragment_length, reg0->read_end_position-reg0->read_start_position-1, reg0->fragment_length, reg0->mapq, reg0->tags);
					}else{
						fprintf(stdout, "%s\t%u\t*\t*\t*\t*\t*\t*\t*\t*\t*\t%d\t%s\n", reg0->read_name, reg0->read_length, reg0->mapq, reg0->tags);
					}
				}

				if(reg0->tags) free(reg0->tags);
				free(reg0);
			}
		}

		for (i = 0; i < p->n_threads; ++i) ri_tbuf_destroy(s->buf[i]);
		free(s->buf);

		if(s->reg)free(s->reg);

		for(int i = 0; i < s->n_sig; ++i){
			ri_sig_t *curS = s->sig[i];
			free(p->events[s->sig[i]->rid].values);
			free(curS->sig);
			free(curS->name);
			free(curS);
		} if(s->sig)free(s->sig); free(s);
	}
    return 0;
}

int ri_map_file(const ri_idx_t *idx, const char *fn, const ri_mapopt_t *opt, int n_threads){
	return ri_map_file_frag(idx, 1, &fn, opt, n_threads);
}

int ri_map_file_frag(const ri_idx_t *idx, int n_segs, const char **fn, const ri_mapopt_t *opt, int n_threads){

	int pl_threads;
	pipeline_mt pl;
	if (n_segs < 1) return -1;
	memset((void *)&pl, 0, sizeof(pipeline_mt));
	pl.n_fp = n_segs;
	pl.n_f = 0; pl.cur_f = 0;
	ri_char_v fnames = {0,0,0};
	kv_resize(char*, 0, fnames, 256);
	find_fast5(fn[0], &fnames);
	pl.f =  fnames.a;
	if(!fnames.n || ((pl.fp = open_sig(pl.f[0])) == 0)) return -1;
	if (pl.fp == 0) return -1;
	pl.fn = fn;
	pl.n_f = fnames.n;
	pl.cur_fp = 1;
	pl.cur_f = 1;
	pl.opt = opt, pl.ri = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	pl_threads = pl.n_threads == 1?1:2;
	pl.su_stop = 0;

	if(opt->flag & RI_M_SEQUENCEUNTIL){
		pl.su_nreads = 0;
		pl.su_nestimations = 0;
		pl.ab_count = 0;
		pl.su_cur = 0;
		pl.su_estimations = (float**)calloc(opt->tn_samples, sizeof(float*));
		for(uint32_t i = 0; i < opt->tn_samples; ++i) pl.su_estimations[i] = (float*)calloc(idx->n_seq, sizeof(float));

		pl.su_c_estimations = (uint32_t*)calloc(idx->n_seq, sizeof(uint32_t));
	}
	
	pl.events.reserve(10000000); //reserve space for 10M reads
	 							 //hack-fix, due to crash likely from internal re-allocation during push_back
								 //this wastes a lot of memory and may still break if more than 10M reads are processed
								 //TODO: fix properly with locks or custom data structure
	kt_pipeline(pl_threads, map_worker_pipeline, &pl, 4);

	if(opt->flag & RI_M_SEQUENCEUNTIL){
		// pl.su_nreads = 0;
		// pl.su_nestimations = 0;
		// pl.su_stop = 0;
		if(pl.su_estimations){
			for(uint32_t i = 0; i < opt->tn_samples; ++i) if(pl.su_estimations[i]) free(pl.su_estimations[i]);
			free(pl.su_estimations);
		}

		if(pl.su_c_estimations)free(pl.su_c_estimations);
	}

	ri_sig_close(pl.fp);
	for(size_t i = 0; i < fnames.n; ++i) free(kv_A(fnames, i));
	kv_destroy(fnames);

	return 0;
}
