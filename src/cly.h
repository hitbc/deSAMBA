/*
 * classify.hpp
 *
 *  Created on: 2018-3-19
 *
 *  Author: fenghe
 */

#ifndef SRC_CLY_H_
#define SRC_CLY_H_

#include <stdint.h>
#include "idx.h"
#include "lib/utils.h"

typedef struct {
	int L_min_matching;
	int thread_num;
	int max_sec_N;
	int out_format;
	int show_anchor;
	FILE * outfile;
	int min_score;
}MAP_opt;

typedef struct
{
	uint32_t offset;
	uint32_t len;
	uint8_t  top;
}CLY_seed;

typedef struct
{
	uint16_t mtch_len;
	int16_t score;
	uint8_t left_len;
	uint8_t left_ED;
	uint8_t rigt_len;
	uint8_t rigt_ED;
}Anchor_map;

typedef struct
{
	//map info
	Anchor_map a_m;
	uint8_t direction;
	//position
	uint64_t global_offset;//global_offset of anchor in all reference string
	uint32_t ref_ID;
	uint32_t ref_offset;//offset in ref
	uint32_t index_in_read;
	void *chain_anchor_pre;
	uint16_t seed_ID;

	//other part
	uint16_t chain_id;
	uint8_t anchor_useless;
	uint8_t duplicate;// use as normal anchor but score 1
}Anchor;

kvec_T(Anchor, Anchor_V)

#define PRIMARY 1
#define SECONDARY 2
#define SUPPLYMENTARY 3

typedef struct
{
	uint32_t ref_ID;
	int32_t  q_t_dis;//equal: q_st - t_st, distance between query and target string, near distance will be regarded as within same chain
	uint32_t sum_score;
	uint32_t anchor_number;
	uint8_t  direction;
	uint8_t  with_top_anchor;
	uint8_t  primary;
	uint8_t  pri_index;

	uint32_t t_st;//target start
	uint32_t t_ed;//target end
	uint32_t q_st;//target start
	uint32_t q_ed;//target end

	uint32_t indel;//total length of insert and delete

	uint32_t chain_id;//
	Anchor *chain_anchor_cur;//target start
}chain_item;

kvec_T(chain_item, chain_set)

typedef struct
{
	kseq_t 		*read;
	chain_set 	hit;
	//anchor list
	uint32_t	fast_classify;
	Anchor_V 	anchor_v;
} cly_r;

typedef struct
{
	uint32_t kmer;
	uint32_t next;
	uint32_t pos;
}sparse_align_HASH;



typedef struct
{
	sparse_align_HASH *sp_hash;
	int q_st;
	int q_ed;
	int t_st;
	int t_ed;
}sp_HASH_table;

typedef struct
{
	uint16_t next;
	uint16_t seed_ID:15,
	s_or_e:1;//the start or end of a seed
}seed_con_hash;//seed conflict hash

typedef struct
{
	uint32_t t_pos;
	uint32_t q_pos;
	uint32_t len;
	uint32_t score;
}spd_match;

kvec_T(spd_match, spd_match_set)

typedef struct
{
	uint8_t * 	bin_read;
	uint32_t	m_bin_read;
	uint64_t* 	kmer_buff;
	uint32_t 	m_kmer_buff;
	CLY_seed*	seed_v;
	uint32_t	m_seed_v;
	//sparse hash buff
	sp_HASH_table *sp_table;
	uint32_t m_sp_table;
	sparse_align_HASH * sp_hash;
	sparse_align_HASH * sa_hash[2];//used for M2
	uint32_t m_sp_hash;
	//seed conflict hash
	seed_con_hash *sc_hash;
	uint32_t m_sc_hash;
	//
	spd_match_set sms;

	int max_read_l;//when max read length < 500, treat as 2-generation reads, otherwise 3 generation reads
}Classify_buff_pool;

#endif /* SRC_CLY_H_ */
