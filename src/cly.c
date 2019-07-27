/*
 * classify.cpp
 *
 *  Created on: 2018-4-2
 *      Author: fenghe
 */
#include "cly.h"

#include "lib/utils.h"
#include <pthread.h>
#include <string.h>

#include "bwt.h"
#include "lib/ksw/ksw2.h"

//same as Bit, but all UNKNOWN characters are used as 'C'
uint8_t CLY_Bit[] =
{
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,+
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 0,/*A*/ 1, 1,  1, 1, 1, 2,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  3, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 0,/*A*/ 1, 1,  1, 1, 1, 2,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  3, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1
};

//return -1 when a > b; return 1 when a < b
static int inline chain_cmp_by_score(const void * a_, const void * b_)
{
	chain_item *a = (chain_item *)a_, *b = (chain_item *)b_;
	if(a->with_top_anchor != b->with_top_anchor)
		return (a->with_top_anchor)?(-1):(1);
	int score_a = a->sum_score + ((a->q_ed - a->q_st)<<1);
	score_a -= (a->indel << 2);
	int score_b = b->sum_score + ((b->q_ed - b->q_st)<<1);
	score_b -= (b->indel << 2);
	if(score_a < score_b)
		return 1;
	if(score_a > score_b)
		return -1;
	return 0;
}

static int inline chain_cmp_by_MEM_score(const void * a_, const void * b_)
{
	chain_item *a = (chain_item *)a_, *b = (chain_item *)b_;
	int score_a = (a->sum_score << 5);// + (a->q_ed - a->q_st);// - a->indel;
	int score_b = (b->sum_score << 5);// + (b->q_ed - b->q_st);// - a->indel;
	if(score_a < score_b)
		return 1;
	if(score_a > score_b)
		return -1;
	return 0;
}

typedef struct
{
	uint16_t chain_ID;//data
	uint16_t next;
}CI_hash;//chain insert hash

static inline void chain_insert_meta(Anchor * anchor, chain_item *c, bool new_chain, int dis_minus)
{
	uint32_t 	ref_l = anchor->ref_offset;
	uint32_t	ref_r = ref_l + anchor->a_m.mtch_len;
	uint32_t 	read_l = anchor->index_in_read;
	uint32_t	read_r = read_l + anchor->a_m.mtch_len;
	if(new_chain)
	{
		anchor->chain_id = c->chain_id;
		anchor->chain_anchor_pre = NULL;
		c->ref_ID = anchor->ref_ID;
		c->direction = anchor->direction;
		c->q_t_dis = anchor->ref_offset - anchor->index_in_read;
		c->t_st = ref_l;
		c->t_ed = ref_r;
		c->q_st = read_l;
		c->q_ed = read_r;
		c->with_top_anchor = !anchor->anchor_useless;
		c->anchor_number = 1;
		c->sum_score = (anchor->duplicate)?1:anchor->a_m.score;
		c->indel = 0;
		c->chain_anchor_cur = anchor;
	}
	else
	{
		anchor->chain_id = c->chain_id;
		c->with_top_anchor |= (!anchor->anchor_useless);
		//delete same anchor
		if(c->q_ed >= read_r)
			return;
		c->t_ed = MAX(ref_r, c->t_ed);
		c->q_ed = read_r;
		anchor->chain_anchor_pre = c->chain_anchor_cur;
		c->chain_anchor_cur = anchor;
		//re-set dis
		c->q_t_dis = anchor->ref_offset - anchor->index_in_read;
		c->indel += dis_minus;
		c->anchor_number++;
		c->sum_score += (anchor->duplicate)?1:anchor->a_m.score;
	}
}

#define MAX_dis_MINUS 30
#define MAX_waiting_len 400
int chain_insert_M1_idx[3] = {0, 1, -1};//0->1->-1

void chain_insert_M1_new_hash_node(CI_hash *ci_hash, uint16_t c_key, int chain_ID, int *kmer_con_index)
{
	if(*kmer_con_index >= 2047)
		return;
	while(0 != ci_hash[c_key].next)//skip nodes with data
	{
		if(ci_hash[c_key].chain_ID == chain_ID)//skip when conflict
			return;
		c_key = ci_hash[c_key].next;
	}

	//store data and next
	ci_hash[c_key].chain_ID = chain_ID;
	ci_hash[c_key].next = *kmer_con_index;
	//new next node
	ci_hash[*kmer_con_index].next = 0;
	ci_hash[*kmer_con_index].chain_ID = 0;
	(*kmer_con_index)++;
}

static inline void chain_insert_M1(
		Anchor * anchor,
		chain_set *c,
		CI_hash *ci_hash,
		int *kmer_con_index)
{
	//get hash key
	uint32_t 	ref_ID = anchor->ref_ID;
	int32_t 	dis = anchor->ref_offset - anchor->index_in_read;
	uint16_t hash_key = (anchor->direction << 8) + ((ref_ID & 0x7) << 5);
	uint16_t hash_dis = dis >> 5;

	bool chain_FOUND = false;
	int chain_ID;
	int dis_minus;
	int hash_offset = 0;
	//search hash table
	for(; hash_offset < 3; hash_offset++)
	{
		uint16_t c_key = hash_key + ((hash_dis + chain_insert_M1_idx[hash_offset]) & 0x1f);
		while(ci_hash[c_key].next != 0)
		{
			chain_ID = ci_hash[c_key].chain_ID;
			if(
					ref_ID == c->a[chain_ID - 1].ref_ID &&//chain with same ref-ID
					(dis_minus = ABS(dis - c->a[chain_ID - 1].q_t_dis)) < MAX_dis_MINUS && //
					ABS_U(c->a[chain_ID - 1].t_ed, anchor->ref_offset) < MAX_waiting_len//waiting length at most 500 bp
			)
			{//store, new_chain = false; chain_ID keep same
				chain_FOUND = true;
				break;
			}
			c_key = ci_hash[c_key].next;
		}
		if(chain_FOUND == true)
			break;
	}
	//store
	if(chain_FOUND == false)
	{
		if(anchor->anchor_useless == true && c->n > 200)
		{
			anchor->chain_id = 65535;
			return;
		}
		uint16_t c_key = hash_key + (hash_dis & 0x1f);
		chain_insert_M1_new_hash_node(ci_hash, c_key, c->n + 1, kmer_con_index);
		chain_item *new_c;
		kv_pushp_2(chain_item, c, new_c);
		new_c->chain_id = c->n - 1;
		chain_insert_meta(anchor, new_c, true, dis_minus);
	}
	else
	{
		if(hash_offset != 0)
		{
			uint16_t c_key = hash_key + (hash_dis & 0x1f);
			chain_insert_M1_new_hash_node(ci_hash, c_key, chain_ID, kmer_con_index);
		}
		chain_insert_meta(anchor, c->a + chain_ID - 1, false, dis_minus);
	}
}

static inline void chain_insert_M2(Anchor * anchor, chain_set *c)
{
	chain_item *c_s = c->a, *c_e = c_s + c->n;
	uint8_t 	direction = anchor->direction;
	uint32_t 	ref_ID = anchor->ref_ID;
	int32_t 	dis = anchor->ref_offset - anchor->index_in_read;
	int dis_minus;
	for(;c_s < c_e; c_s++)
	{
		if(
				c_s->direction == direction &&
				c_s->ref_ID == ref_ID &&
				(dis_minus = ABS(dis - c_s->q_t_dis)) < MAX_dis_MINUS &&
				ABS_U(c_s->t_ed, anchor->ref_offset) < MAX_waiting_len)//same chain
		{
			chain_insert_meta(anchor, c_s, false, dis_minus);
			return;
		}
	}
	chain_item *new_c;
	kv_pushp_2(chain_item, c, new_c);
	new_c->chain_id = c->n - 1;
	chain_insert_meta(anchor, new_c, true, dis_minus);
}

#define MAX_hit_len 500
void resolve_tree(cly_r *results)
{
	int kmer_con_index;
	results->hit.n = 0;
	Anchor * anchor_st = results->anchor_v.a, *anchor_ed = anchor_st + results->anchor_v.n;
	if(results->anchor_v.n < 300)
		for(Anchor * anchor = anchor_st; anchor < anchor_ed; anchor++)
			chain_insert_M2(anchor, &(results->hit));
	else
	{
		kmer_con_index = 512;
		CI_hash ci_hash[2048];
		memset(ci_hash, 0, 512*sizeof(CI_hash) + 1);
		for(Anchor * anchor = anchor_st; anchor < anchor_ed; anchor++)
			chain_insert_M1(anchor, &(results->hit), ci_hash, &kmer_con_index);
	}
	if(results->hit.n > 1)
		qsort(results->hit.a, results->hit.n, sizeof(chain_item), chain_cmp_by_score);
	int rst_num = MIN(5, results->hit.n);//delete result without top anchor
	while(rst_num < results->hit.n && results->hit.a[rst_num].with_top_anchor == true)
		rst_num++;
	results->hit.n = rst_num;
}

void transIntoBytes( char *str, int len, uint8_t *byteFormat)
{
	for(int i=0; i < len; ++i)
		byteFormat[i] = CLY_Bit[(int8_t)str[i]];
}

extern uint64_t kmerMask[33];

//store kmers into a buff, set kmers to be 0 when it failed the filter. kmer 0 (all A)is no use normally
void store_kmers(
		uint8_t *bin_read,
		uint32_t kmer_len,
		uint8_t l_e_kmer,
		int single_base_max,
		uint64_t* 	kmer_buff)
{
	//filter
	int base_ACGT[4] = {0};
	for(uint32_t i = 0; i < l_e_kmer ; i++)
		base_ACGT[bin_read[i]]++;
	uint8_t *bin_read_p = bin_read;
	uint64_t MASK = kmerMask[l_e_kmer];
	uint64_t kmer =  bit2_nextKmer_init(bin_read_p, l_e_kmer);
	//i== 0
	int filter_failed = false;
	if(		base_ACGT[0] >= single_base_max ||
			base_ACGT[1] >= single_base_max ||
			base_ACGT[2] >= single_base_max ||
			base_ACGT[3] >= single_base_max)
		filter_failed = true;
	kmer = bit2_nextKmerMASK(bin_read_p++, kmer, l_e_kmer);
	kmer_buff[0] = (filter_failed)?0:kmer;
	for(uint32_t i = 1; i < kmer_len; i++)
	{
		//FILTER
		filter_failed = false;
		base_ACGT[bin_read_p[-1]]--;
		base_ACGT[bin_read_p[l_e_kmer - 1]]++;
		if(		base_ACGT[0] >= single_base_max ||
				base_ACGT[1] >= single_base_max ||
				base_ACGT[2] >= single_base_max ||
				base_ACGT[3] >= single_base_max)
			filter_failed = true;
		//Normal
		kmer = bit2_nextKmerMASK(bin_read_p++, kmer, l_e_kmer);
		kmer_buff[i] = (filter_failed)?0:kmer;
	}
}

//DEBUG function:
void debug_print_all_e_kmer(
		uint32_t read_index,
		uint32_t l_seed_v,
		uint64_t* kmer_buff,
		uint32_t l_kmer_buff,
		uint8_t* e_kmer)
{
	fprintf(stderr,
			"read number: %d\t"
			" island number : %d\n",
			read_index,
			l_seed_v);
	for(uint32_t i = 0; i < l_kmer_buff; i++)
	{
		//search e_k_table
		uint64_t kmer = kmer_buff[i];//200ns
		uint8_t rst = (e_kmer[kmer >> 3] >> (7 - (kmer & 0x7))) & 0x1;
		fprintf(stderr,"%d",rst);
	}
	fprintf(stderr,"\n");
}

void debug_print_ACGTstring(
		uint8_t* string, uint32_t length)
{
	//return;
	fprintf(stderr,"length: %d\t\t", length);
	for(uint32_t i = 0; i < length; i++)
		fprintf(stderr, "%c", "ACGT"[string[i]]);
	fprintf(stderr,"\n");
}

//for forward, load: uni_offset to uni_offset + length - 1
//for reverse, load: uni_offset - length + 1 to uni_offset
static void inline get_ref(
		uint8_t* unitig_str,
		uint8_t* ref_str,
		uint64_t uni_offset,
		uint32_t length,
		bool isForward)
{
	uint64_t offset = uni_offset >> 2;
	uint8_t  odd    = uni_offset & 0x3;
	if(isForward)
		for(uint32_t k = 0; k < length; k++)
			switch(odd)
			{
				case 0: ref_str[k] = (unitig_str[offset] >> 6) & 0x3 ; odd = 1; break;
				case 1: ref_str[k] = (unitig_str[offset] >> 4) & 0x3 ; odd = 2; break;
				case 2: ref_str[k] = (unitig_str[offset] >> 2) & 0x3 ; odd = 3; break;
				case 3: ref_str[k] = (unitig_str[offset]     ) & 0x3 ; odd = 0; offset++; break;
			}
	else
		for(uint32_t k = 0; k < length; k++)
			switch(odd)
			{
				case 0: ref_str[k] = (unitig_str[offset] >> 6) & 0x3 ; odd = 3; offset--; break;
				case 1: ref_str[k] = (unitig_str[offset] >> 4) & 0x3 ; odd = 0; break;
				case 2: ref_str[k] = (unitig_str[offset] >> 2) & 0x3 ; odd = 1; break;
				case 3: ref_str[k] = (unitig_str[offset]     ) & 0x3 ; odd = 2; break;
			}
}

//sp + search_l --> global_offset + sp_offset
//global_offset: means the global position of the TRUE sp
//sp: the sp contains SA INFO
static inline UNITIG * get_uni(
		DA_IDX * idx,
		uint64_t bwt_pos,
		int search_l,
		uint64_t *global_offset,
		uint32_t *uni_offset_)
{
	UNITIG *uni = idx->unitig_v.a + idx->bt.sa_taxon[bwt_pos >> SA_OFF].unitig_ID;//get uni
	uint32_t uni_offset = idx->bt.sa_taxon[bwt_pos >> SA_OFF].offset + search_l + 1;
	if(search_l > 0) //minus offset first, add UNI second
		for(;uni_offset >= uni->length;)
		{
			 uni_offset -= (uni->length + 1);
			 uni++;
		}
	else 			 //minus UNI first, add offset second
		for(;uni_offset < 0;)
		{
			  uni--;
			  uni_offset += (uni->length + 1);
		}
	REF_POS * r_p = idx->r_p_v.a + uni->ref_list;	//get ref
	*global_offset = r_p->global_offset + uni_offset;
	*uni_offset_ = uni_offset;//offset in the uni
	return uni;
}

#define LV_ERROR 4
#define LV_BASE LV_ERROR
#define RETURN_PRO(i) 							\
{												\
	ref[ref_length] = old_ref_end;				\
	query[query_length] = old_query_end;		\
	return i;									\
}
//version:4.1:landauVishkin_semi_extend_ED
//LV direct match, global alignment, when match are found return edit distance, else return query_length(LV_FALSE)
//only edit distance, no CIGAR
//max err num are defined by 'LV_ERROR'
static int32_t inline lv_extd(
		uint8_t* ref,
		int32_t ref_length,
		uint8_t* query,
		int32_t query_length)
{
	if(ref_length < query_length)//swap
	{
		int32_t ref_length_ = ref_length;
		ref_length = query_length;
		query_length = ref_length_;

		uint8_t* ref_ = ref;
		ref = query;
		query = ref_;
	}

	int32_t match_num_data[99];
	int32_t edit_dist_data[99];//99 is enough
	int32_t *mn = match_num_data + LV_BASE + 1, *ed = edit_dist_data + LV_BASE + 1;
	int32_t prev_mn, cur_mn, next_mn;
	int32_t prev_ed, cur_ed, next_ed;

	uint8_t old_ref_end = ref[ref_length], old_query_end = query[query_length];
	ref[ref_length] = '#';//the end of text and pattern can`t be same, besides, they can`t appear in other part of text and pattern
	query[query_length] = '$';
	int32_t best_score = query_length;

    //initialization of match number and edit distance

	for(int i = - LV_BASE - 1; i <= LV_BASE + 1; i++)
	{//check
		mn[i] = -1;
		ed[i] = (i > 0)?(i):(-i);
	}

	for(int i=0; i <= LV_ERROR; i++)
	{
		//step1: init prev, cur and next
        prev_mn = -1;
        cur_mn  = (i - 1);
        next_mn = mn[-i + 1];
		prev_ed = i + 1; cur_ed = i;
		next_ed = ed[-i + 1];

    	for(int j = - i; j <= LV_ERROR; j++)
    	{
            //step2: dynamic programming
    		if(cur_mn + j < ref_length-1)
    		{//CONSIDER BOTH MN AND ED
    			int MAX_mn_ed = cur_mn  + 1 - cur_ed;
    			mn[j] = cur_mn + 1;
    			ed[j] = cur_ed  + 1;
    			if(MAX_mn_ed  < next_mn + 1 - next_ed)
    			{
    				mn[j] = next_mn + 1;
    				ed[j] = next_ed + 1;
    				MAX_mn_ed = next_mn - next_ed;
    			}
    			if(MAX_mn_ed  < prev_mn 	- prev_ed)
    			{
    				mn[j] = prev_mn + 1;
    				ed[j] = prev_ed + 1;
    			}
    		}
    		else{
    			int MAX_mn_ed = cur_mn  - cur_ed;
    			mn[j] = cur_mn;
    			ed[j] = cur_ed  + 1;
    			if(MAX_mn_ed  < prev_mn - prev_ed)
    			{
    				mn[j] = prev_mn;
    				ed[j] = prev_ed + 1;
    				MAX_mn_ed = prev_mn - prev_ed;
    			}
    			if(MAX_mn_ed  < next_mn + 1 - next_ed)
    			{
    				mn[j] = next_mn + 1;
    				ed[j] = next_ed + 1;
    			}
            }
    		int mn_j = MIN(mn[j], query_length);
    		mn_j = MIN(mn_j, ref_length - j);
            //step4: in-line match
   			for(; ref[mn_j + j] == query[mn_j]; mn_j++);
   			mn[j] = mn_j;
            //step5: return when one match are found
            if(query[mn_j] == '$' || ref[mn_j + j] == '#')
            {
            	best_score = MIN(ed[j] - 1, best_score);
            	if(j <= i + 1)
            		RETURN_PRO(best_score);//otherwise research
            }
            //step6: update prev, cur and next
            prev_mn = cur_mn;	cur_mn = next_mn;	next_mn = mn[j + 2];
    		prev_ed = cur_ed;	cur_ed = next_ed; 	next_ed = ed[j + 2];
        }
    }
    RETURN_PRO(best_score);
}

typedef struct
{
	uint8_t *bin_read;
	uint32_t read_L;
	uint16_t seed_ID;
	bool direction;
}SEED_INFO;

typedef struct
{
	int match_len;
	uint64_t sp;
	uint64_t sa_sp;
	int sa_sp_l;
	int kmer_index;
	int read_offset;
}MEM_rst;

static inline void get_new_ed(
		uint32_t *e_d,//edit distance
		uint32_t *len_,//compare length
		uint32_t  *l_mem_ext,//additional MEM length extended
		int32_t  q_off,
		uint64_t t_off,
		uint32_t l_read,
		uint8_t  *q_b,
		uint8_t  *t_b,
		bool is_FWD)
{
	uint8_t q_buff[13], *q = q_buff;
	uint8_t t_buff[13], *t = t_buff;//query string
	uint32_t len, max_len;//query and ref length

	//load read
	if(is_FWD)
	{
		if(q_off < 0)
			q_off = 0;
		max_len = q_off;
		len = MIN(12, max_len);
		for(uint8_t k = 0; k < len; k++)
			q[k] = q_b[q_off - k];
	}
	else
	{
		max_len = l_read - q_off;
		len = MIN(12, max_len);
		q = q_b + q_off;
	}
	//load ref
	get_ref(t_b, t, t_off, len, !is_FWD);
	if(len > 0 && t[0] == q[0])
	{//MEM
		int mtc;
		do
		{
			for(mtc = 0; mtc < len; mtc++)
				if(t[mtc] != q[mtc])
					break;
			if(mtc > 0)
			{
				*l_mem_ext += mtc;
				max_len -= mtc;
				len = MIN(12, max_len);
				if(is_FWD)
				{
					q_off -= mtc;
					t_off -= mtc;
					for(uint8_t k = 0; k < len; k++)
						q[k] = q_b[q_off - k];
				}
				else
				{
					t_off += mtc;
					q += mtc;
				}
				get_ref(t_b, t, t_off, len, !is_FWD);
			}
		}while(mtc > 0);
	}
	//do lv
	*e_d  = lv_extd(t, len, q, len);
	*len_ = len;
}

#define MIN_S_1 12//
#define MIN_S_2 20//
void getUniStr_suf(uint8_t *unitig_str, uint8_t *ref, uint64_t offset, uint32_t len)
{
	for(uint8_t k = 0; k < len; k++)
		ref[k] = unitig_str[offset + k];
}

#define LV_L 12//max query length used in LV

int32_t map_seed(DA_IDX * idx, MEM_rst* m_r, SEED_INFO *s_i, Anchor_V *anchor_v)
{
	MAPQ 	 *m_Q 	= &(idx->mapQ);
	bwt		 *bt 	= &(idx->bt);
	uint64_t b_p 	= m_r->sp;//current postion in the BWT string
	int32_t  q_off 	= m_r->read_offset;
	uint32_t l_m 	= m_r->match_len;//match length of MEM
	uint8_t  *q_b 	= s_i->bin_read;//bin read
	uint8_t  *t_b 	= idx->ref_bin.a;

	UNITIG 	 *uni = NULL;
	uint32_t u_off;//offset in a unitig
	uint64_t t_off;//target offset in all reference string
	uint32_t l_pre, l_suf;
	uint32_t d_pre, d_suf;//edit distance of pre and suf lv
	int32_t  s = 0,	max_s = 0;;//score
	///part 1: mapping seed to UNITIG:CORE
	do{
		///step1: get the prefix
		uint8_t q_pre[LV_L + 1];
		uint8_t t_pre[LV_L + 1];
		uint8_t *q_suf;
		uint8_t t_suf[LV_L + 1];
		//at most read 12 bp
		if(q_off < -1)
			printf(" ");
		l_pre = MIN(q_off + 1, LV_L);
		for(uint8_t k = 0; k < l_pre; k++)
			q_pre[k] = q_b[q_off - k];
		int s_l = 0;//search length
		//locate, stop when get the location of the seed or get enough data to run LV
		if(m_r->sa_sp != MAX_uint64_t)//sa_sp exist
			uni = get_uni(idx, m_r->sa_sp, m_r->sa_sp_l, &t_off, &u_off);
		else
		{
			uint8_t c;
			uint64_t new_sp;
			while(1)
			{
				if((b_p & SA_MASK) == 0)//when SA can be get
					break;
				//get the new sp and c
				c = 0xff;
				new_sp = occ(bt, b_p, &c) + bt->rank[c];
				if(c == 4)//reach the "begin" of a UNITIG
					break;
				t_pre[s_l++] = c;
				b_p = new_sp;
				//when get enough data to compare
				if(s_l >= l_pre)
					break;
			}
			//load ref_prefix when SA are found
			if((b_p & SA_MASK) == 0)
				uni = get_uni(idx, b_p, s_l, &t_off, &u_off);
			else
				l_pre = s_l;
		}
		//get ref_prefix when UNI are known
		if(uni != NULL)
		{
			if(uni->length < MIN_UNI_L)// ignore too short UNITIG
				break;
			l_pre = MIN(l_pre, u_off);
			get_ref(t_b, t_pre, t_off - 1, l_pre, false);
		}
		d_pre = lv_extd(t_pre, l_pre, q_pre, l_pre);
		s = m_Q->Q_MEM[l_m] + m_Q->Q_LV[d_pre][l_pre];
		//fail
		if(s < MIN_S_1 && l_pre == LV_L && uni == NULL)//full-lv but score not enough
		{
			s = 0;
			break;
		}
		///step2: get suffix
		if(uni == NULL)
		{
			//get the UNI
			while(b_p & SA_MASK)
			{
				uint8_t c = 0xff;
				b_p = occ(bt, b_p, &c) +  bt->rank[c];
				s_l++;
			}
			uni = get_uni(idx, b_p, s_l, &t_off, &u_off);
			if(uni->length < MIN_UNI_L)// ignore too short UNITIG
			{
				s = 0;
				break;
			}
		}
		int32_t q_off_r = q_off + l_m + 1;
		uint32_t l_max_suf = MIN(uni->length - u_off - l_m, s_i->read_L - q_off_r);//suf_l = MIN(LV_L, s_i->read_L - q_off, uni->length - u_off - l_m)

		if(l_max_suf != 0)
		{
			l_suf = MIN(l_max_suf, LV_L);
			q_suf = q_b + q_off_r;
			get_ref(t_b, t_suf, t_off + l_m, l_suf, true);
			//MEM
			if(t_suf[0] == q_suf[0])
			{
				int mtc;
				do
				{
					for(mtc = 0; mtc < l_suf; mtc++)
						if(t_suf[mtc] != q_suf[mtc])
							break;
					if(mtc > 0)
					{
						l_m += mtc;
						s = m_Q->Q_MEM[l_m] + m_Q->Q_LV[d_pre][l_pre];
						l_max_suf -= mtc;
						l_suf = MIN(l_max_suf, LV_L);
						q_suf += mtc;
						get_ref(t_b, t_suf, t_off + l_m, l_suf, true);
					}
				}while(mtc > 0);
			}
			d_suf = lv_extd(t_suf, l_suf, q_suf, l_suf);
			s += m_Q->Q_LV[d_suf][l_suf];
		}
		else
			l_suf = d_suf = 0;
		//when FAIL
		if(s <= MIN_S_2 && l_suf == LV_L)
		{
			s = 0;
			break;
		}
	}while(0);

	///part 2: mapping seed to reference
	//store mapping detail informations
	if(s > 0)
	{
		Anchor_map a_m = {l_m, s, l_pre, d_pre, l_suf, d_suf};
		REF_POS * r_p_s = idx->r_p_v.a + uni->ref_list, *r_p_e = idx->r_p_v.a + (uni + 1)->ref_list;//get ref
		bool ref_search_l = (l_pre < LV_L || d_pre == 0)?true:false;//when edit distance is 0, means that extension is not over
		bool ref_search_r = (l_suf < LV_L || d_suf == 0)?true:false;
		uint8_t duplicate = false;
		if(r_p_e - r_p_s > 35)
		{
			duplicate = true;
			r_p_s += ((anchor_v->n && 0xf))*((r_p_e - r_p_s)>>4);//get 10 random results (by current anchor_v->n)
			r_p_e = MIN(r_p_s + 10, r_p_e);
		}
		//if(r_p_e - r_p_s > 40 && anchor_v->n > 30)//ignore results when too many of then
		//	r_p_e = r_p_s;
		for(REF_POS * c_r_p = r_p_s; c_r_p < r_p_e; c_r_p++)
		{
			//for:anchor reach the bound of unitig:
			uint32_t ed_l, ed_r;//edit distance
			uint32_t len_l, len_r;//edit distance
			uint32_t l_m_ext_l = 0, l_m_ext_r;
			if(ref_search_l || ref_search_r)
			{
				if(ref_search_l)
				{
					get_new_ed(&ed_l, &len_l, &l_m_ext_l,
							q_off, c_r_p->global_offset + u_off - 1,
							s_i->read_L, q_b, t_b,true);
					a_m.left_len = len_l;
					a_m.left_ED = ed_l;
				}
				a_m.mtch_len = l_m + l_m_ext_l;
				if(ref_search_r)
				{
					l_m_ext_r = 0;
					get_new_ed(&ed_r, &len_r, &l_m_ext_r,
							q_off + l_m + 1, c_r_p->global_offset + u_off + l_m,
							s_i->read_L, q_b, t_b,false);
					a_m.rigt_len = len_r;
					a_m.rigt_ED = ed_r;
					a_m.mtch_len += l_m_ext_r;
				}
				a_m.score =
						m_Q->Q_MEM[a_m.mtch_len] +
						m_Q->Q_LV[a_m.left_ED][a_m.left_len] +
						m_Q->Q_LV[a_m.rigt_ED][a_m.rigt_len];
				if(a_m.score < MIN_S_2)
					continue;
			}
			max_s = MAX(max_s, a_m.score);
			Anchor * a; kv_pushp_2(Anchor, anchor_v, a);
			//store anchors
			a->direction = s_i->direction;
			a->index_in_read = q_off + 1 - l_m_ext_l;//+2
			a->global_offset = c_r_p->global_offset + u_off - l_m_ext_l;
			a->ref_ID = c_r_p->ref_ID;
			a->ref_offset = a->global_offset - idx->r_i_v.a[a->ref_ID].seq_offset;
			a->a_m = a_m;
			a->seed_ID = s_i->seed_ID;
			a->duplicate = duplicate;
		}
	}
	return max_s;
}

int CLY_seed_cmp(const void*a_,const void*b_)
{
	return ((CLY_seed *)b_)->len - ((CLY_seed *)a_)->len;
}

typedef struct{
	//seed v
	CLY_seed* seed_v_f;
	uint32_t l_seed_v_f;
	uint8_t * bin_read;
	uint64_t* kmer;
	uint32_t direction;
	uint32_t total_score;
}SEARCH_DIR;//search direction

static inline int get_exist_kmer(
		uint8_t* e_kmer_1,
		uint8_t* e_kmer_2,
		uint64_t kmer,
		uint64_t hash_mask)
{
	//filter out kmers when less complex
	if(kmer == 0)
		return 0;
	uint64_t hash_1 = hash64_1(kmer) & hash_mask;//36bit mask
	uint8_t rst_1 = (e_kmer_1[hash_1 >> 3] >> (7 - (hash_1 & 0x7))) & 0x1;
	if(rst_1 == 0)
		return 0;
	uint64_t hash_2 = hash64_2(kmer) & hash_mask;//36bit mask
	uint8_t rst_2 = (e_kmer_2[hash_2 >> 3] >> (7 - (hash_2 & 0x7))) & 0x1;
	return rst_2;
}

#define STEP_EK 3 //step used for exist kmer
uint32_t search_exist_kmer(
		uint64_t* kmer_v,
		uint32_t l_kmer_v,
		CLY_seed *seed_v,
		uint8_t* e_1,
		uint8_t* e_2,
		uint64_t hash_mask)
{
	uint32_t l_seed_v = 0;
	for(uint32_t i = STEP_EK - 1; i < l_kmer_v; i += STEP_EK)
	{
		//search e_k_table	//200ns
		if(get_exist_kmer(e_1, e_2, kmer_v[i], hash_mask) == 1)
		{//kmer exist
			uint32_t offset = i, len = 1;
			for(int j = 1; j < STEP_EK; ++j)
			{// 1 2 3 4	//search backward
				if(get_exist_kmer(e_1, e_2, kmer_v[i - j], hash_mask) == 1)
				{
					offset--;
					len++;
				}
				else
					break;
			}
			for(int j = 1; i + j < l_kmer_v; ++j)
			{// 1 2 3 4 //200ns //search forward, as long as possible
				if(get_exist_kmer(e_1, e_2, kmer_v[i + j], hash_mask) == 1)
				{
					len++;
					if(len > 60)
					{
						i += 50;
						break;
					}
				}
				else
					break;
			}
			seed_v[l_seed_v].offset = offset;
			seed_v[l_seed_v].len = len;//store rst
			l_seed_v ++;
			i = offset + len; //reset i, ignore covered kmers
		}
	}
	return l_seed_v;
}

//get-seed-vector
#define SEED_RANGE 100
void get_seed_vector(E_KMER* ek, uint8_t * bin_read, uint64_t* kmer_buff, uint32_t l_kmer_buff,
		CLY_seed* seed_v, uint32_t direction, SEARCH_DIR *search_dir_)
{
	store_kmers(bin_read, l_kmer_buff, ek->len_e_kmer, ek->single_base_max, kmer_buff);
	uint32_t l_seed_v = search_exist_kmer(kmer_buff, l_kmer_buff, seed_v, ek->e_kmer0, ek->e_kmer1, ek->e_kmer_hash_mask);
	uint32_t total_score = 0;
	//top 1 in 100 bp will be labeled
	int max_index = 0;	uint32_t max_length = 0;	uint32_t index_end = SEED_RANGE;
	for(uint32_t m = 0; m < l_seed_v; m++)
	{
		if(seed_v[m].offset < index_end)
		{
			if(max_length < seed_v[m].len)
			{
				max_length = seed_v[m].len;
				max_index = m;
			}
			seed_v[max_index].top = false;
		}
		else
		{
			seed_v[max_index].top = true;
			index_end += SEED_RANGE;
			total_score += max_length;
			max_length = 0;
			max_index = m;
			max_length = seed_v[m].len;
		}
	}
	seed_v[max_index].top = true;
	total_score += max_length;//the last one
	SEARCH_DIR search_dir = {seed_v, l_seed_v, bin_read, kmer_buff, direction, total_score};
	search_dir_[0] = search_dir;
}

void getIsland(kseq_t *read, Classify_buff_pool *buff, E_KMER* ek, SEARCH_DIR *search_dir)
{
	uint8_t len_e_kmer = ek->len_e_kmer;//set-to-be 20bp for whole meta
	uint32_t read_len = read->seq.l;
	char *useRead = read->seq.s;
	BUFF_REALLOC(buff->bin_read, buff->m_bin_read, read_len + read_len);
	BUFF_REALLOC(buff->kmer_buff, buff->m_kmer_buff, read_len + read_len);
	BUFF_REALLOC(buff->seed_v, buff->m_seed_v, (read_len >> 1));//4 base pair per seed = read_len*2/4
	//get the number of total kmers///forward and reverse, store kmer in kmer buff and search from exist kmer table
	uint32_t l_kmer_buff 	= read_len - len_e_kmer + 1;
///forward
	uint8_t * bin_read_F 	= buff->bin_read;
	uint64_t* kmer_buff_F 	=  buff->kmer_buff;
	CLY_seed* seed_v_F 		= buff->seed_v;
	for(uint32_t k = 0; k < read_len; ++k)
		bin_read_F[k] 		= CLY_Bit[(int8_t)useRead[k]];//get bin read string
	get_seed_vector(ek, bin_read_F, kmer_buff_F, l_kmer_buff, seed_v_F, FORWARD, search_dir);
///reverse
#ifndef CONSIDER_BOTH_ORIENTATION
	uint8_t * bin_read_R 	= buff->bin_read + read_len;
	uint64_t* kmer_buff_R 	= buff->kmer_buff + read_len;
	CLY_seed* seed_v_R 		= buff->seed_v + (read_len >> 2);
	for(uint32_t k = 0; k < read_len; ++k)
		bin_read_R[read_len - k - 1] = 3 - bin_read_F[k];//get the bin reverse read string
	get_seed_vector(ek, bin_read_R, kmer_buff_R, l_kmer_buff, seed_v_R, REVERSE, search_dir + 1);
	if(search_dir[0].total_score < search_dir[1].total_score)
	{
		SEARCH_DIR swap = search_dir[0];
		search_dir[0] = search_dir[1];
		search_dir[1] = swap;
	}
#endif
}

int anchor_cmp_index(const void*a_,const void*b_)
{
	return ((Anchor *)a_)->index_in_read - ((Anchor *)b_)->index_in_read;
}//small to big

int anchor_cmp_score(const void*a_,const void*b_)
{
	return ((Anchor *)b_)->a_m.score - ((Anchor *)a_)->a_m.score;
}//big to small

//insert fail return 0
typedef struct{
	uint64_t *set;
	int l, m;
}SP_SET;

static inline int sp_set_insert(uint64_t node, SP_SET *sp_set)
{
	if(sp_set->l == sp_set->m)
		sp_set->l = 0;//reset
	uint64_t *set = sp_set->set;
	int i = 0;
	for(; i < sp_set->l; i++)
		if(set[i] == node)
			return 0;
	set[i] = node;
	sp_set->l++;
	return 1;
}

#define TAXNOMY_SET_MAX 50
typedef struct
{
	uint32_t tax[TAXNOMY_SET_MAX];
	int l;
}TAXNOMY_SET;

int tax_set_insert(uint64_t node, TAXNOMY_SET *tax_set)
{
	if(tax_set->l == TAXNOMY_SET_MAX)
		tax_set->l = 0;//reset
	uint32_t *set = tax_set->tax;
	int i = 0;
	for(; i < tax_set->l; i++)
		if(set[i] == node)
			return 0;
	set[i] = node;
	tax_set->l++;
	return 1;
}

#define NEXT_SP_UNKNOW 0
#define NEXT_SP_NULL -1
typedef struct
{
	uint64_t sp;
	uint8_t next_char;
	uint8_t search_char;
}SP_LIST;

static int inline MEM_rst_cmp_by_match_len(const void*a_,const void*b_)
{
	return ((MEM_rst *)b_)->match_len - ((MEM_rst *)a_)->match_len;
}

void print_sp_list(SP_LIST *sp,int match_len)
{
	fprintf(stderr, "->%ld:%c:%c:%d ",
			sp->sp,
			"ACTG#"[sp->next_char],
			"ACTG#"[sp->search_char],
			match_len + 20);
}

static inline void bwt_single_search(
		bwt* bt,
		uint64_t sp,
		uint8_t* string,//TODO: change "sa_sp: to accelerate speed.
		int max_match_len,
		SP_SET *sp_set, MEM_rst* mem_rst)
{
	uint64_t *rank = bt->rank;
	uint64_t new_sp, sa_sp = MAX_uint64_t;
	int match_len = 0, sa_sp_l = 0;
	while(1)
	{
		if(match_len >= max_match_len)
			break; //when reach the end of read
		if((sp & SA_MASK) == 0)
		{
			sa_sp = sp;
			sa_sp_l = 0;
		}
		else
			sa_sp_l--;
		uint8_t c = 0xff;
		new_sp = occ(bt, sp, &c) + rank[c];
		if(c != *string) //end
			break;
		match_len++;
		string--;
		//store new sps
		if(sp_set_insert(new_sp, sp_set) == 0)
		{
			mem_rst->match_len = -1000;
			return;
		}
		sp = new_sp;
	}
	mem_rst->sp = sp;
	mem_rst->match_len = match_len;
	mem_rst->sa_sp = sa_sp;
	mem_rst->sa_sp_l = sa_sp_l;
}

//, int MIN_MATCH_LEN, int MAX_MATCH_LEN
//given a string, get all match at least MIN_MATCH_LEN
//#define mem_rst_MAX_search 5
static inline int bwt_MEM_search(
		bwt* bt,
		uint8_t* string,
		uint64_t pre_v,//prefix value
		int max_rst,//max number of result
		int l_min_mth,
		int l_max_mth,
		SP_SET *sp_set,
		MEM_rst* mem_rst)
{
	int n_rst = 0;
	uint64_t *rank = bt->rank;
	uint64_t sp = bt->hash_index[pre_v], ep = bt->hash_index[pre_v + 1], new_sp, new_ep;
	string -= L_PRE_IDX;
	int match_len = L_PRE_IDX;
	uint8_t c;
	while(1)
	{
		c = *string; //get a new char//line 1
		string--;//LINE 2
		new_sp = rank[c] + occ(bt, sp, &c);//re-get new sp/ep
		new_ep = rank[c] + occ(bt, ep, &c);
		if(match_len >= l_min_mth - 1)
		{
			if(new_sp + max_rst >= new_ep)
				break;
			if(match_len >= l_max_mth)
				return 0;
		}
		if(new_sp + 1 >= new_ep)
			break;//LINE 3
		match_len++;
		sp = new_sp;
		ep = new_ep;//
	}
	if(new_sp >= new_ep)//no match long enough
		return 0;
	if(new_sp + 1 == new_ep){//only one
		if(sp_set_insert(new_sp, sp_set) == 0)
			return 0;
		bwt_single_search(bt, new_sp, string, MAX(0, l_max_mth - match_len), sp_set, mem_rst + n_rst);
		mem_rst[n_rst].match_len += match_len + 1;//MEM search part + single search part, should be "+="
		if(mem_rst[n_rst].match_len >= l_min_mth)
			n_rst++;
	}
	else
	{//over one
		for(uint64_t c_sp = new_sp; c_sp < new_ep; c_sp++)
		{
			if(sp_set_insert(c_sp, sp_set) == 0)//try to store 20mer
				continue;
			//single search
			bwt_single_search(bt, c_sp, string, MAX(0, l_max_mth - match_len), sp_set, mem_rst + n_rst);
			mem_rst[n_rst].match_len += match_len + 1;
			if(mem_rst[n_rst].match_len >= l_min_mth)
				n_rst++;
		}
	}
	return n_rst;
}

//debug function
void show_anchor(Anchor *anchor)
{
	Anchor* anchor_b = anchor; //Anchor* anchor_e = anchor_b + p_rst->anchor_l;
	fprintf(stderr,
			"++"
			"%d\t"
			"%d\t"
			"%d\t"
			"%c\t",
			anchor_b->index_in_read,
			anchor_b->ref_ID,
			anchor_b->a_m.score,
			(anchor_b->direction)?'F':'R');
	fprintf(stderr,
			"%d\t"
			"%d\t"
			"%d\t"
			"%d\t"
			"%d\n",
			anchor_b->a_m.mtch_len,
			anchor_b->a_m.left_len,
			anchor_b->a_m.left_ED,
			anchor_b->a_m.rigt_len,
			anchor_b->a_m.rigt_ED);
}

#define MEM_search_FAST 2
#define MIN_MEM_LEN_FAST 21
void fast_classify(
		DA_IDX* idx,
		SEARCH_DIR *s_d,//search direction
		uint32_t read_len, cly_r *results)
{
	uint8_t  	l_ek = idx->ek.len_e_kmer;//length of e kmer
	int 		min_index = MIN_MEM_LEN_FAST - l_ek;
	uint64_t	*kmer = s_d->kmer;
	uint8_t 	*bin_read = s_d->bin_read;
	uint64_t 	sp_set_BUFF[500];
	SP_SET 		sp_set = {sp_set_BUFF, 0, 500};
	MEM_rst 	m_r[MEM_search_FAST];
	CLY_seed	*sv_b = s_d->seed_v_f, *sv_e = sv_b + s_d->l_seed_v_f;//begin and end of seed vector
	SEED_INFO 	s_i = {bin_read, read_len, 0, s_d->direction};

	//extend each kmer island
	for(CLY_seed* c_sv = sv_b; c_sv < sv_e; c_sv++){
		if(c_sv->top == false)
			continue;
		sp_set.l = 0;
		s_i.seed_ID = c_sv - sv_b;
		uint32_t a_b_idx = results->anchor_v.n;
		for(int j = c_sv->len - 1; j >= min_index;)
		{
			//STEP1: get seeds
			int kmer_index = c_sv->offset + j;
			uint64_t prefixValue = kmer[kmer_index] & PRE_IDX_MASK;
			int string_index = kmer_index + l_ek - 1;
			int c_mem_rst_num = bwt_MEM_search(
					&(idx->bt), bin_read + string_index,
					prefixValue, MEM_search_FAST,
					MIN_MEM_LEN_FAST - 1, string_index,
					&sp_set, m_r);
			if(c_mem_rst_num == 0)
			{
				j-= 2;
				continue;
			}
			j -= 3;
			//STEP2: map seeds
			int max_score = 0;
			for(MEM_rst* c_mr = m_r; c_mr < m_r + c_mem_rst_num; ++c_mr)
			{
				c_mr->read_offset = string_index - c_mr->match_len;
				int c_score = map_seed(idx, c_mr, &s_i, &(results->anchor_v));
				max_score = MAX(c_score, max_score);//TODO:::
			}
			if(max_score > 35)
				j -= 7;
			if(max_score > 256)
			{
				if(max_score > 512)//skip next and break, to handle 2rd-G reads
					c_sv++;
				break;//at least 50 match
			}
		}
		//STEP3: delete results
		Anchor * a_b = results->anchor_v.a + a_b_idx;
		Anchor * a_e = results->anchor_v.a + results->anchor_v.n;
		int top_score = 35;
		for(Anchor * anc_c = a_b; anc_c < a_e; anc_c++)
			top_score = MAX(top_score, anc_c->a_m.score);
		for(Anchor * anc_c = a_b; anc_c < a_e; anc_c++)
			anc_c->anchor_useless = (anc_c->a_m.score < top_score)?true:false;
	}
}

#define MEM_search_SLOW 8
#define MIN_MEM_LEN_SLOW 20
void slow_classify(DA_IDX* idx, SEARCH_DIR *search_dir, uint32_t read_len, cly_r *results)
{
	int l_ek = idx->ek.len_e_kmer;//length of e kmer
	uint8_t * bin_read = search_dir->bin_read;
	uint64_t* kmer = search_dir->kmer;
	CLY_seed* sv_f = search_dir->seed_v_f;
	uint64_t sp_set_BUFF[500];
	SP_SET sp_set = {sp_set_BUFF, 0, 500};
	MEM_rst mem_rst[MEM_search_SLOW*800 + 1];
	int mem_rst_num;
	SEED_INFO seed_info = {bin_read, read_len, 0, search_dir->direction};
	for(uint32_t i = 0; i < search_dir->l_seed_v_f; i++)
	{
		if((int)(sv_f[i].len) < 3 && sv_f->top == false)
			continue;
		int min_match_len = MIN(MIN_MEM_LEN_SLOW - 1, l_ek + 1);//at most 19
		sp_set.l = 0;
		mem_rst_num = 0;
		//STEP1: get seed
		for(int j = sv_f[i].len - 1; j >= 1; j-=2)
		{
			int k_idx = sv_f[i].offset + j;//kmer index
			uint64_t pre_v = kmer[k_idx] & PRE_IDX_MASK;
			int s_idx = k_idx + l_ek - 1;//string index
			int c_mem_rst_num = bwt_MEM_search(
					&(idx->bt),
					bin_read + s_idx,
					pre_v,
					MEM_search_SLOW,
					min_match_len,
					s_idx,
					&sp_set,
					mem_rst + mem_rst_num);
			for(int i = mem_rst_num; i < mem_rst_num + c_mem_rst_num; i++)
			{
				mem_rst[i].read_offset = k_idx + l_ek - 1 -
				mem_rst[i].match_len;
			}
			mem_rst_num += c_mem_rst_num;
		}
		//STEP2: map seed
		if(mem_rst_num == 0)
			continue;
		if(mem_rst_num > 1)
			qsort(mem_rst, mem_rst_num, sizeof(MEM_rst), MEM_rst_cmp_by_match_len);//sort
		seed_info.seed_ID = i;
		uint32_t a_b_idx = results->anchor_v.n;
		int max_search = MIN(mem_rst_num, MEM_search_SLOW);
		for(MEM_rst* c_mem_rst = mem_rst;c_mem_rst < mem_rst + max_search;++c_mem_rst)
			map_seed(idx, c_mem_rst, &seed_info, &(results->anchor_v));
		//STEP3: delete results
		Anchor * a_b = results->anchor_v.a + a_b_idx;
		Anchor * a_e = results->anchor_v.a + results->anchor_v.n;
		int top_score = 35;
		for(Anchor * anc_c = a_b; anc_c < a_e; anc_c++)
			top_score = MAX(top_score, anc_c->a_m.score);
		for(Anchor * anc_c = a_b; anc_c < a_e; anc_c++)
			anc_c->anchor_useless = (anc_c->a_m.score < top_score)?true:false;
	}
	results->fast_classify = false;
}

void print_bin_string(uint8_t * string, int len)
{
	for(int i = 0; i < len; i++)
		fprintf(stderr, "%c", "ACGT#$"[string[i]]);
	fprintf(stderr, "\n");
}

extern uint64_t kmerMask[33];
void sp_hash_idx(
		chain_set *hit,
		sp_HASH_table *sp_table,
		sparse_align_HASH * sp_hash_buff,
		SEARCH_DIR *search_dir,
		uint32_t read_l,
		int kmer_shift)
{
	int block_bg_q, block_ed_q;
	int both_dir = 0;
	for(int i = 0; i < hit->n; i++)//check dir
	{
		both_dir |= (hit->a[i].direction == FORWARD)?0x2:0x1;
		if(both_dir == 3)
			break;
	}
	for(int i = 2; i >= 1; i--)
	{
		if((i & both_dir) == 0)
			continue;
		uint32_t direction = (i == 1)?REVERSE:FORWARD;
		SEARCH_DIR *c_search_dir = ((search_dir->direction == direction)?0:1) + search_dir;
		int block_ed = read_l >> 8;
		//set info for each block
		int buff_id = direction, buff_id_ed = (block_ed << 1) + direction;
		int q_st = -35, t_st = 0;
		for(; buff_id <= buff_id_ed; buff_id += 2)
		{
			sp_table[buff_id].q_st = q_st;
			sp_table[buff_id].q_ed = q_st + 255 + 70;
			sp_table[buff_id].sp_hash = sp_hash_buff + (buff_id << 9);
			sp_table[buff_id].t_st = t_st;
			sp_table[buff_id].t_ed = t_st + 255;
			q_st += 256;
			t_st += 256;
		}
		//first one and last one
		sp_table[direction].q_st = 0;
		sp_table[buff_id_ed].t_ed = read_l - 10;
		sp_table[buff_id_ed].q_ed = read_l - 10;
		for(int block_ID = 0; block_ID <= block_ed; block_ID++)
		{
			int buff_id = (block_ID << 1) + direction;
			sparse_align_HASH *sa_hash = sp_table[buff_id].sp_hash;
			memset(sa_hash, 0, 256 * sizeof(sparse_align_HASH));
			int kmer_con_index = 256;
			uint16_t key;
			block_bg_q = sp_table[buff_id].q_st;
			block_ed_q = sp_table[buff_id].q_ed;
			int block_len_q = block_ed_q - block_bg_q + 1;
			uint64_t *read_kmer = c_search_dir->kmer + block_bg_q;
			///step1: build hash index
			for(int i = 0; i < block_len_q; i += 4)//travel reference
			{
				uint32_t kmer = read_kmer[i] >> kmer_shift;
				key = kmer & 0xff;//mask 1024, 10 bit
				while(sa_hash[key].next != 0)//skip block with data
					key = sa_hash[key].next;
				//store data:
				sa_hash[key].kmer = kmer;
				sa_hash[key].pos = i;
				//new next
				sa_hash[key].next = kmer_con_index;
				sa_hash[kmer_con_index].next = 0;
				kmer_con_index++;
			}
		}
	}
}

void sc_hash_idx(seed_con_hash *sc_hash, uint32_t n_sc_hash, chain_set *hit)
{
	memset(sc_hash, 0, 256 * sizeof(seed_con_hash));
	chain_item *st_hit = hit->a, *ed_hit = st_hit + hit->n;
	int sc_con_index = 256;
	for(chain_item * c_h = st_hit; c_h < ed_hit; c_h++)//travel reference
		for(int i = 1; i >= 0; i--)
		{	//1:left, 0:right
			uint16_t c_key = ((i == 1)?(c_h->t_st - c_h->q_st):(c_h->t_ed - c_h->q_ed)) & 0xff;
			//skip block with data
			while(sc_hash[c_key].next != 0)
				c_key = sc_hash[c_key].next;
			//store data:
			sc_hash[c_key].seed_ID = c_h - st_hit + 1;
			sc_hash[c_key].s_or_e = i;
			sc_hash[c_key].next = sc_con_index;
			//get next
			sc_hash[sc_con_index++].next = 0;
		}
}

#define S_A_KEMR_L 9
#define MIN_SCORE_MEM 12//S_A_KEMR_L + 3
int sparse_align(
		sparse_align_HASH * sp_hash,
		uint8_t *ref_p,
		int ref_len,
		uint32_t *q_v,
		uint32_t *t_v
		)
{
	uint8_t *ref = ref_p;
	int count = 0;
	uint16_t key;
	uint64_t MASK = kmerMask[S_A_KEMR_L];
	uint64_t kmer = bit2_nextKmer_init(ref, S_A_KEMR_L);
	//sparse alignment
	for(int i = 0; i < ref_len; i++)
	{
		kmer = bit2_nextKmerMASK(ref, kmer, S_A_KEMR_L);
		ref ++;
		key = kmer & 0xff;
		while(sp_hash[key].next != 0)
		{
			if(kmer == sp_hash[key].kmer)
			{
				q_v[count] = sp_hash[key].pos;
				t_v[count++] = i;
				if(count > 99)
					return count;
			}
			key = sp_hash[key].next;
		}
	}
	return count;
}

bool search_anchor(uint8_t *q, uint8_t *t, int q_st, int t_st, MAPQ *mapQ)
{
	//MAPQ 	 *m_Q 	= &(idx->mapQ);
	//left
	int l_l = MIN(12, q_st);
	l_l = MIN(l_l, t_st);
	int s_l = lv_extd(q + q_st - l_l, l_l, t + t_st - l_l, l_l);
	int s_r = lv_extd(q + q_st + S_A_KEMR_L, 12, t + t_st + S_A_KEMR_L, 12);
	int score = mapQ->Q_LV[s_l][l_l] + mapQ->Q_LV[s_r][12];
	if(score >= 30)
		return true;
	return false;
}

//combine two chains into one
bool combine_chain(
		chain_item * c_st,
		int chain_ID,
		seed_con_hash *sc_hash,
		int dis,//for match, distance between q_ed and t_ed
		bool isleft,
		int c_q_pos,//q_st for right, q_ed for left
		chain_item ** combined_chain)
{
	uint16_t key = (dis)&0xff;
	chain_item *c, *c_h = c_st + chain_ID;
	//search sc_hash table:
	while(sc_hash[key].next != 0)
	{
		uint16_t seed_ID = sc_hash[key].seed_ID;
		c = c_st + seed_ID - 1;
		int dis_con = (isleft)?(c->t_ed - c->q_ed):(c->t_st - c->q_st);
		int q_pos_con = (!isleft)?(c->q_st):(c->q_ed - S_A_KEMR_L);//when left to get end; right to get start
		if(
				dis == dis_con &&//same dis
				c_h != c &&//not same chain
				isleft != sc_hash[key].s_or_e && //left + right or right + left
				ABS_U(c_q_pos, q_pos_con) < 8 &&
				c_h->ref_ID == c->ref_ID &&//same ref
				c_h->direction == c->direction &&//same direction
				c->sum_score != 0 &&//not blank
				seed_ID - 1 > chain_ID//rst afterward
		  )
		{
			//combine
			c_h->sum_score += c->sum_score;
			c_h->anchor_number += c->anchor_number;
			c_h->indel += c->indel;
			c_h->q_st = MIN(c_h->q_st, c->q_st);
			c_h->t_st = MIN(c_h->t_st, c->t_st);
			c_h->q_ed = MAX(c_h->q_ed, c->q_ed);
			c_h->t_ed = MAX(c_h->t_ed, c->t_ed);
			c->sum_score = 0;//delete this chain
			c->t_st = c->t_ed = c->q_st = c->q_ed = 0;
			*combined_chain = c;
			return true;
		}
		key = sc_hash[key].next;
	}
	return false;
}

inline static int MEM_search(uint8_t *q, uint8_t *t, bool forward, int max)
{
	int len = 0;
	if(forward)
		for(;len < max && *q++ == *t++; len++);
	else
		for(;len < max && *q-- == *t-- ; len++);
	return len;
}

///------------------------------------END Get SCORE M1--------------------------------------------------


int get_middle_score(
		sparse_align_HASH * sp_hash,
		uint8_t *ref_p,

		uint8_t *read_p,
		int ref_len,
		int read_base,
		int *dis_minus,
		int *overflow
)
{
	int score = 0;
	*overflow = 0;
	uint8_t *ref = ref_p;//begin from 3
	uint16_t key;
	uint64_t MASK = kmerMask[S_A_KEMR_L];
	uint64_t kmer = bit2_nextKmer_init(ref, S_A_KEMR_L);
	//sparse alignment
	int i = 0;
	for(; i < ref_len; i++, ref++)
	{
		kmer = bit2_nextKmerMASK(ref, kmer, S_A_KEMR_L);
		key = kmer & 0xff;
		while(sp_hash[key].next != 0)
		{
			if(kmer == sp_hash[key].kmer)
			{
				int pos_in_read = sp_hash[key].pos - read_base;
				if(pos_in_read > 0 && ABS_U(i, pos_in_read + *dis_minus) < 25)
				{
					//MEM search, left
					*dis_minus = i - pos_in_read;
					int left_len  = MEM_search(read_p + pos_in_read - 1, ref_p + i - 1, false, 3);
					int right_len = MEM_search(read_p + pos_in_read + S_A_KEMR_L, ref_p + i + S_A_KEMR_L, true, 1000);
					int c_score = left_len + S_A_KEMR_L + right_len;
					if(c_score >= MIN_SCORE_MEM)
					{
						//fprintf(stderr, "%d %d\n", pos_in_read, c_score);
						score += c_score - S_A_KEMR_L;
						i += right_len + S_A_KEMR_L + 1;
						ref += right_len + S_A_KEMR_L + 1;
						if(i > ref_len)
							*overflow = i - ref_len;
						else
							kmer = bit2_nextKmer_init(ref + 1, S_A_KEMR_L);
						break;
					}
				}
			}
			key = sp_hash[key].next;
		}
	}
	return score;
}

int sparse_align_middle_extern(Anchor * c_a, sp_HASH_table * sp_table, DA_IDX * idx, uint8_t *read_p)
{
	int score = 0, over_flow;
	uint64_t t_offset = idx->r_i_v.a[c_a->ref_ID].seq_offset;
	while(c_a != NULL)
	{
		//add score
		score += c_a->a_m.mtch_len - S_A_KEMR_L + 1;
		//fprintf(stderr,"\n%d %d\n", c_a->index_in_read, c_a->a_m.mtch_len);//DEBUG+++
		//add middle score
		Anchor * pre_a = c_a->chain_anchor_pre;
		if(pre_a != NULL)
		{
			int pre_mch = pre_a->a_m.mtch_len;
			int total_ref_len = c_a->ref_offset - (pre_a->ref_offset + pre_mch);
			if(total_ref_len > MIN_SCORE_MEM)
			{
				uint32_t extern_bg = pre_a->index_in_read + pre_mch;
				uint32_t extern_ed = c_a->index_in_read - S_A_KEMR_L;
				int block_bg = extern_bg >> 8;
				int block_ed = extern_ed >> 8;
				int dis_minus = 0;
				uint8_t ref[800];//get ref
				uint64_t ref_offset = pre_a->ref_offset + t_offset + pre_mch;
				int ref_st = 0;
				get_ref(idx->ref_bin.a, ref, ref_offset, total_ref_len + 50, true);
				ref[total_ref_len + 50] = 254;//end mark
				for(int block_ID = block_bg; block_ID <= block_ed; block_ID ++)
				{
					sp_HASH_table * c_sp_table = sp_table + (block_ID << 1) + c_a->direction;
					int ref_len = (block_ID == block_ed)?(total_ref_len - ref_st):(c_sp_table->t_ed - extern_bg);
					if(ref_len <= 0)
						break;
					score += get_middle_score(
							c_sp_table->sp_hash, ref + ref_st,
							read_p + extern_bg, ref_len, extern_bg - c_sp_table->q_st,
							&dis_minus, &over_flow);
					//reset begin for query"extern_bg" and ref"ref_st"
					extern_bg = c_sp_table->t_ed + 1 + over_flow;
					ref_st += ref_len + over_flow;
				}
			}
		}
		c_a = pre_a;
	}
	return score;
}

int sparse_align_edge_extern_left(
		sp_HASH_table *sp_table,
		chain_item * c_st,
		int chain_ID,
		uint8_t *read_v,
		DA_IDX * idx,
		int l_read,
		seed_con_hash *sc_hash)
{
	int score = 0;
	chain_item * c_h = c_st + chain_ID, *combined_chian;
	uint64_t t_offset = idx->r_i_v.a[c_h->ref_ID].seq_offset;
	int ori_t_minus_q = c_h->t_st - c_h->q_st;//used to change q begin to t begin, ref_offset = t_offset + ori_t_minus_q + block_bg_t;
	int dis = ori_t_minus_q;
	int block_bg = c_h->q_st >> 8;
	int extern_overflow = 0;
	int last_score = 0;
	uint32_t q_v[100], t_v[100], last_i = 0;
	for(int block_ID = block_bg; block_ID != -1; block_ID --)
	{
		sp_HASH_table *c_sp_table = sp_table + (block_ID << 1) + c_h->direction;
		sparse_align_HASH * sp_hash = c_sp_table->sp_hash;
		int block_bg_t = block_ID << 8;
		int block_ed_t = MIN(block_bg_t + 255, c_h->q_st) - extern_overflow;
		extern_overflow = 0;
		int old_q_pos = c_h->q_st;//used to get the distance between two kmer matches
		uint64_t ref_offset = t_offset + ori_t_minus_q + block_bg_t;
		int block_len_t = block_ed_t - block_bg_t + 1;
		uint8_t ref_p[400];//get ref
		uint8_t *ref = ref_p + 30;//get 30 more ref in the right and left
		get_ref(idx->ref_bin.a, ref_p, ref_offset - 30, block_len_t + S_A_KEMR_L + 60, true);
		ref_p[0] = ref_p[ block_len_t + S_A_KEMR_L + 60] = 254;
		int ref_pos_basic = ref_offset - t_offset;
		int read_pos_basic = c_sp_table->q_st;//block_bg_q;;
		int count = sparse_align(sp_hash, ref, block_len_t, q_v, t_v);
		int MEM_last_site = c_h->q_st + 1;//used to skip already scored block
		int search_stat = 0;//0:MAX_minus_dis == 30; 1: MAX_minus_dis == 45; 2: break loop normal; 3: accept kmer found; 4:combined chain
		for(; search_stat < 2; search_stat ++)
		{
			int MAX_minus_dis = (search_stat == 0)?30:45;//try 30 at first time, when no rst, try 55
			for(int i = count - 1; i != -1; i--)//SEARCH ALL COUNTS
			{
				int c_q_pos = (read_pos_basic + q_v[i]);
				int c_dis = ref_pos_basic + t_v[i] - c_q_pos;//current t pos - q pos
				int abs_minus_dis = ABS_U(dis, c_dis);
				int minus_pos = old_q_pos - c_q_pos;//accept a kmer
				if((abs_minus_dis < 7 && minus_pos < 60) ||//condition 1
						((abs_minus_dis < MAX_minus_dis) && search_anchor(read_v, ref, c_q_pos, t_v[i], &(idx->mapQ))))//condition 2
				{
					//search combined chain
					if(minus_pos == 4 &&
						combine_chain(c_st, chain_ID, sc_hash, c_dis, true, c_q_pos, &combined_chian) == true)//search sc_hash
					{
						dis = ori_t_minus_q = c_h->t_st - c_h->q_st;
						block_ID = (c_h->q_st >> 8) + 1;
						score += sparse_align_middle_extern(combined_chian->chain_anchor_cur, sp_table, idx, read_v);//combine score
						extern_overflow = 0;
						//fprintf(stderr, "--%d\n", last_score);
						score -= last_score;
						last_score = 0;
						search_stat = 3;//3+1 = 4
						break;
					}
					old_q_pos = c_q_pos;
					last_i = i;
					//get score for new MEM
					if(MEM_last_site > c_q_pos)
					{
						//change states
						dis = c_dis;

						//get score
						int left_len  = MEM_search(read_v + c_q_pos - 1, ref + t_v[i] - 1, false, 1000);
						int right_len = MEM_search(read_v + c_q_pos + S_A_KEMR_L, ref + t_v[i] + S_A_KEMR_L, true, 3);
						int c_score = left_len + S_A_KEMR_L + right_len;
						search_stat = 2;//skip second loop,//2+1 = 3
						if(c_score >= MIN_SCORE_MEM)
						{
							if(t_v[i] < left_len)
								extern_overflow = left_len - t_v[i] + S_A_KEMR_L + 1;
							MEM_last_site = c_q_pos - left_len - 1;
							score += c_score - S_A_KEMR_L;
							last_score = c_score;
							//fprintf(stderr,"%d %d\n", c_q_pos, c_score);//DEBUG+++
						}
					}
				}
			}
		}
		if(search_stat == 3)
		{
			int total_dis_minus = c_h->t_st - c_h->q_st;
			c_h->q_st = read_pos_basic + q_v[last_i];
			c_h->t_st = ref_pos_basic + t_v[last_i];
			total_dis_minus -= (c_h->t_st - c_h->q_st);
			c_h->indel += ABS(total_dis_minus);
			//update
		}
		else if(search_stat < 3 && block_len_t > 200)//when waiting length over 200 bp, stop extension
			break;
	}
	return score;
}

int sparse_align_edge_extern_right(
		sp_HASH_table *sp_table,
		chain_item * c_st,
		int chain_ID,
		uint8_t *read_v,
		DA_IDX * idx,
		int l_read,
		seed_con_hash *sc_hash)
{
	int score = 0;
	chain_item * c_h = c_st + chain_ID, *combined_chian;
	uint64_t t_offset = idx->r_i_v.a[c_h->ref_ID].seq_offset;
	int ori_t_minus_q = c_h->t_ed - c_h->q_ed;//used to change q begin to t begin, ref_offset = t_offset + ori_t_minus_q + block_bg_t;
	int dis = ori_t_minus_q;
	int block_bg = c_h->q_ed >> 8;
	int block_ed = (l_read + 255 - S_A_KEMR_L) >> 8;
	int extern_overflow = 0;
	int last_score = 0;
	uint32_t q_v[100], t_v[100], last_i = 0;
	for(int block_ID = block_bg; block_ID != block_ed; block_ID ++)
	{
		sp_HASH_table *c_sp_table = sp_table + (block_ID << 1) + c_h->direction;
		sparse_align_HASH * sp_hash = c_sp_table->sp_hash;
		int block_bg_t = block_ID << 8, block_ed_t = block_bg_t + 255;
		block_bg_t = MAX(block_bg_t, c_h->q_ed) + extern_overflow;
		int old_q_pos = c_h->q_ed;//used to get the distance between two kmer matches
		uint64_t ref_offset = t_offset + ori_t_minus_q + block_bg_t;
		extern_overflow = 0;
		int block_len_t = block_ed_t - block_bg_t + 1;
		uint8_t ref_p[400];//get ref
		uint8_t *ref = ref_p + 30;//get 30 more ref in the right and left
		get_ref(idx->ref_bin.a, ref_p, ref_offset - 30, block_len_t + S_A_KEMR_L + 60, true);
		ref_p[0] = ref_p[ block_len_t + S_A_KEMR_L + 60] = 254;
		int ref_pos_basic = ref_offset - t_offset;
		int read_pos_basic = c_sp_table->q_st;//block_bg_q;;
		int count = sparse_align(sp_hash, ref, block_len_t, q_v, t_v);
		int MEM_last_site = c_h->q_ed - 1;
		int search_stat = 0;//0:MAX_minus_dis == 30; 1: MAX_minus_dis == 45; 2: break loop normal; 3: accept kmer found; 4:combined chain
		for(; search_stat < 2; search_stat ++)
		{
			int MAX_minus_dis = (search_stat == 0)?30:45;//try 30 at first time, when no rst, try 55
			for(int i = 0; i != count; i++)
			{
				int c_q_pos = (read_pos_basic + q_v[i]);
				int c_dis = ref_pos_basic + t_v[i] - c_q_pos;
				int abs_minus_dis = ABS_U(dis, c_dis);
				int minus_pos = c_q_pos - old_q_pos;
				if((abs_minus_dis < 7 && minus_pos < 60) ||//condition 1
						((abs_minus_dis < MAX_minus_dis) && search_anchor(read_v, ref, c_q_pos, t_v[i], &(idx->mapQ))))//condition 2
				{
					//search combined chain
					if(minus_pos == 4 &&
						combine_chain(c_st, chain_ID, sc_hash, c_dis, false, c_q_pos, &combined_chian) == true)//search sc_hash
					{
						dis = ori_t_minus_q = c_h->t_ed - c_h->q_ed;
						block_ID = (c_h->q_ed >> 8) - 1;
						score += sparse_align_middle_extern(combined_chian->chain_anchor_cur, sp_table, idx, read_v);//combine score
						extern_overflow = 0;
						//fprintf(stderr, "--%d\n", last_score);
						score -= last_score;
						last_score = 0;
						search_stat = 3;
						break;
					}
					old_q_pos = c_q_pos;
					last_i = i;
					//get score for new MEM
					if(c_q_pos > MEM_last_site)
					{
						//change states
						dis = c_dis;

						//get score
						int left_len  = MEM_search(read_v + c_q_pos - 1, ref + t_v[i] - 1, false, 3);
						int right_len = MEM_search(read_v + c_q_pos + S_A_KEMR_L, ref + t_v[i] + S_A_KEMR_L, true, 1000);
						int c_score = left_len + S_A_KEMR_L + right_len;
						search_stat = 2;//skip second loop
						if(c_score >= MIN_SCORE_MEM)
						{
							if(t_v[i] + S_A_KEMR_L + right_len > block_len_t)//???
								extern_overflow = t_v[i] + S_A_KEMR_L + right_len - block_len_t;//???
							MEM_last_site = c_q_pos + S_A_KEMR_L + right_len + 1;
							score += c_score - S_A_KEMR_L;
							last_score = c_score;
							//fprintf(stderr,"%d %d %d\n", c_q_pos, c_score, ref_pos_basic + t_v[i]);//DEBUG+++
						}
					}
				}
			}
		}
		if(search_stat == 3)
		{
			int total_dis_minus = c_h->t_ed - c_h->q_ed;
			c_h->q_ed = read_pos_basic + q_v[last_i] + S_A_KEMR_L;
			c_h->t_ed = ref_pos_basic + t_v[last_i] + S_A_KEMR_L;
			total_dis_minus -= (c_h->t_ed - c_h->q_ed);
			c_h->indel += ABS(total_dis_minus);
		}
		else if(search_stat < 3 && block_len_t > 200)//when waiting length over 200 bp, stop extension
			break;
	}
	return score;
}

void get_score_M1(SEARCH_DIR *search_dir, Classify_buff_pool *buff,
		DA_IDX * idx, uint32_t l_read, cly_r *results, seed_con_hash *sc_hash)
{
	uint32_t n_sp_table = ((l_read >> 8) + 1) << 1;//512 per block
	uint32_t n_sp_hash = n_sp_table << 9;//512 per block
	BUFF_REALLOC(buff->sp_table,      buff->m_sp_table,      n_sp_table);
	BUFF_REALLOC(buff->sp_hash,       buff->m_sp_hash,       n_sp_hash);
	sp_hash_idx(&(results->hit), buff->sp_table, buff->sp_hash, search_dir, l_read, ((idx->ek.len_e_kmer - S_A_KEMR_L) << 1));

	//Part 3:read hash
	chain_item *st_hit = results->hit.a;
	for(int i = 0; i < results->hit.n; i++)
	{
		if(st_hit[i].sum_score == 0)//skip result already combined by other hit
			continue;
		//fprintf(stderr,"\n\n\n%s %d\n\n\n\n", idx->r_i_v.a[st_hit[i].ref_ID].ref_name, st_hit[i].t_st);//DEBUG+++
		//if(504221 == st_hit[i].t_st)
		//	fprintf(stderr, " ");
		SEARCH_DIR *c_search_dir = ((search_dir->direction == st_hit[i].direction)?0:1) + search_dir;
		int score = sparse_align_middle_extern(st_hit[i].chain_anchor_cur, buff->sp_table, idx, c_search_dir->bin_read);
		score += sparse_align_edge_extern_left(buff->sp_table, st_hit, i, c_search_dir->bin_read, idx, l_read, sc_hash);
		score += sparse_align_edge_extern_right(buff->sp_table, st_hit, i, c_search_dir->bin_read, idx, l_read, sc_hash);
		st_hit[i].sum_score = score;//score -> sort score
	}
}
///------------------------------------END Get SCORE M1--------------------------------------------------

///------------------------------------Get SCORE M2--------------------------------------------------

extern uint64_t kmerMask_bit[65];

uint32_t hash_size[20] = {
		0x00001, 	0x00002, 	0x00004, 	0x00008,
		0x00010, 	0x00020, 	0x00040, 	0x00080,
		0x00100, 	0x00200,	0x00400, 	0x00800,
		0x01000, 	0x02000, 	0x04000, 	0x08000,
		0x10000,	0x20000, 	0x40000, 	0x80000};

//sa_hash = malloc(sizeof(sparse_align_HASH)*0x100000)//at most 2M blocks, 0.25M for hash table; 1.75M for max read length
int build_hash_table_M2(SEARCH_DIR *search_dir, chain_set *hit, int q_len, sparse_align_HASH **sa_hash)
{
	//INIT
	int both_dir = 0;
	for(int i = 0; i < hit->n; i++)//check direction, forward: 2; reverse: 1
	{
		both_dir |= (hit->a[i].direction == FORWARD)?0x2:0x1;
		if(both_dir == 3)
			break;
	}
	//get hash table
	int key_len = 10;//hash key length
	for(; key_len < 18; key_len++)
		if(hash_size[key_len] >= q_len)
			break;
	uint64_t MASK = kmerMask[S_A_KEMR_L];//9 kmer mask
	uint64_t KEY_MASK = kmerMask_bit[key_len];
	//for both direction
	for(int c_dir = 2; c_dir >= 1; c_dir--)
	{
		if((c_dir & both_dir) == 0)
			continue;
		//c_dir == 2, forward; i == 1 reverse
		uint32_t direction = (c_dir == 1)?REVERSE:FORWARD;
		SEARCH_DIR *c_search_dir = ((search_dir->direction == direction)?0:1) + search_dir;
		sparse_align_HASH *c_sa_hash = (c_dir == 2)?sa_hash[0]:sa_hash[1];
		//clean hash table
		int kmer_con_index = hash_size[key_len];
		for(int index = 0; index < kmer_con_index; index++)
			c_sa_hash[index].next = 0;
		uint8_t * c_q_bin = c_search_dir->bin_read;
		uint64_t kmer = bit2_nextKmer_init(c_q_bin, S_A_KEMR_L);
		for(uint32_t c_pos = 0; c_pos < q_len - S_A_KEMR_L + 1; c_pos++)
		{
			//get kmer
			kmer = bit2_nextKmerMASK(c_q_bin, kmer, S_A_KEMR_L);
			c_q_bin++;
			//get key:
			uint32_t next = kmer &KEY_MASK;
			//get block
			while(c_sa_hash[next].next != 0)
				next = c_sa_hash[next].next;
			//new a new node
			uint32_t new_node_ID = kmer_con_index++;
			c_sa_hash[new_node_ID].kmer = kmer;
			c_sa_hash[new_node_ID].next = 0;
			c_sa_hash[new_node_ID].pos = c_pos;
			c_sa_hash[next].next = new_node_ID;
		}
	}
	return 	key_len;
}

#define OVER_SEARCH_M2 50

void sdp_match_forward(uint32_t q_bg, uint32_t q_ed, uint8_t *q_str, uint8_t *t_str, uint32_t t_len, int key_len,
		sparse_align_HASH *sa_hash, spd_match_set* sms, uint32_t t_st)
{
	uint64_t KEY_MASK = kmerMask_bit[key_len];
	uint32_t t_kmer_num = t_len - S_A_KEMR_L + 1;
	uint64_t MASK = kmerMask[S_A_KEMR_L];
	//build hash table for string 1
	uint8_t *c_t_str = t_str + 4;//skip first three bit, c_t_str = t_str + i
	uint64_t kmer = bit2_nextKmer_init(c_t_str, S_A_KEMR_L);
	//sparse alignment
	for(int i = 4; i < t_kmer_num; i ++, c_t_str++)
	{
		kmer = bit2_nextKmerMASK(c_t_str, kmer, S_A_KEMR_L);
		if((i & 0x03) != 0)
			continue;
		uint32_t next = sa_hash[kmer &KEY_MASK].next;
		while(next != 0)
		{
			if(sa_hash[next].kmer == kmer)//filter 1
			{
				//get read pos
				uint32_t q_pos = sa_hash[next].pos;
				if(q_pos >= q_bg && q_bg <= q_ed)//filter 2
				{
					//search backward 3 bit
					int back_len = MEM_search(q_str + q_pos - 1, c_t_str - 1, false, 4);
					//back_len == 4:skip
					//back_len < 4; back_len = back_len;
					if(back_len < 4)//filter 4
					{
						//search forward
						uint32_t max_search = q_ed - q_pos - 1;
						max_search = MIN(max_search, t_len - i - 1) + OVER_SEARCH_M2;
						int forward_len = MEM_search(q_str + q_pos + S_A_KEMR_L, c_t_str + S_A_KEMR_L, true, max_search);
						//store if total score >= 4
						int total_len = back_len + forward_len + 1;
						if(total_len >= 4)
						{
							spd_match * spd_p;
							kv_pushp_2(spd_match, sms, spd_p);
							spd_p->len = total_len;
							spd_p->q_pos = q_pos - back_len;
							spd_p->t_pos = i - back_len + t_st;
						}
					}
				}
			}
			//last step:: reset next
			next = sa_hash[next].next;
		}
	}
}

void sdp_match_reverse(uint32_t q_bg, uint32_t q_ed, uint8_t *q_str, uint8_t *t_str, uint32_t t_len, int key_len,
		sparse_align_HASH *sa_hash, spd_match_set* sms, uint32_t t_st)
{
	uint64_t KEY_MASK = kmerMask_bit[key_len];
	uint32_t t_kmer_num = t_len - S_A_KEMR_L + 1;
	uint8_t  MOVE = ((S_A_KEMR_L << 1) - 2);
	//build hash table for string 1
	uint8_t * c_t_str = t_str + t_len - S_A_KEMR_L - 4;//skip first three bit
	uint64_t kmer = bit2_preKmer_init(c_t_str, S_A_KEMR_L);
	//sparse alignment
	for(int i = 4; i < t_kmer_num; i++, c_t_str--)
	{
		kmer = bit2_preKmerMOVE(c_t_str, kmer);
		if((i & 0x03) != 0)
			continue;
		uint32_t next = sa_hash[kmer &KEY_MASK].next;
		while(next != 0)
		{
			if(sa_hash[next].kmer == kmer)//filter 1, kmer match
			{
				//get read pos
				uint32_t q_pos = sa_hash[next].pos;
				if(q_pos >= q_bg && q_bg <= q_ed)//filter 2
				{
					//search forward 4 bit
					int forward_len = MEM_search(q_str + q_pos + S_A_KEMR_L, c_t_str + S_A_KEMR_L, true, 4);
					//forward_len == 4:skip
					//forward_len < 4; forward_len = forward_len;
					if(forward_len < 4)//filter 3
					{
						//search forward
						uint32_t max_search = q_pos;//query restriction: 0 ~ q_pos
						max_search = MIN(max_search, c_t_str - t_str) + OVER_SEARCH_M2;//target restriction:
						int back_len = MEM_search(q_str + q_pos - 1, c_t_str - 1, false, max_search);
						//store if total score >= 4
						int total_len = back_len + forward_len + 1;
						if(total_len >= 4)//store the index of end edge
						{
							spd_match * spd_p;
							kv_pushp_2(spd_match, sms, spd_p);
							spd_p->len = total_len;
							spd_p->q_pos = q_pos - back_len;//point to the left most point
							spd_p->t_pos = c_t_str - t_str - back_len + t_st;// point to the left most point
						}
					}
				}
			}
			//last step:: reset next
			next = sa_hash[next].next;
		}
	}
}

//s_st means the global (from the begin of the sequence) position of target
void sdp_match(uint32_t q_bg, uint32_t q_ed, uint8_t *q_str, uint8_t *t_str, uint32_t t_len, int key_len,
		sparse_align_HASH *sa_hash, spd_match_set* sms, uint32_t t_st, bool isForward)
{
	//global
	uint64_t KEY_MASK = kmerMask_bit[key_len];
	uint32_t t_kmer_num = t_len - S_A_KEMR_L + 1;
	if(isForward)
	{
		//get kmer from target
		uint64_t MASK = kmerMask[S_A_KEMR_L];
		//build hash table for string 1
		uint8_t *c_t_str = t_str + 4;//skip first three bit, c_t_str = t_str + i
		uint64_t kmer = bit2_nextKmer_init(c_t_str, S_A_KEMR_L);
		//sparse alignment
		for(int i = 4; i < t_kmer_num; i ++, c_t_str++)
		{
			kmer = bit2_nextKmerMASK(c_t_str, kmer, S_A_KEMR_L);
			if((i & 0x03) != 0)
				continue;
			uint32_t next = sa_hash[kmer &KEY_MASK].next;
			while(next != 0)
			{
				if(sa_hash[next].kmer == kmer)//filter 1
				{
					//get read pos
					uint32_t q_pos = sa_hash[next].pos;
					if(q_pos >= q_bg && q_bg <= q_ed)//filter 2
					{
						//search backward 3 bit
						int back_len = MEM_search(q_str + q_pos - 1, c_t_str - 1, false, 4);
						//back_len == 4:skip
						//back_len < 4; back_len = back_len;
						if(back_len < 4)//filter 4
						{
							//search forward
							uint32_t max_search = q_ed - q_pos - 1;
							max_search = MIN(max_search, t_len - i - 1) + OVER_SEARCH_M2;
							int forward_len = MEM_search(q_str + q_pos + S_A_KEMR_L, c_t_str + S_A_KEMR_L, true, max_search);
							//store if total score >= 4
							int total_len = back_len + forward_len + 1;
							if(total_len >= 4)
							{
								spd_match * spd_p;
								kv_pushp_2(spd_match, sms, spd_p);
								spd_p->len = total_len;
								spd_p->q_pos = q_pos - back_len;
								spd_p->t_pos = i - back_len + t_st;
							}
						}
					}
				}
				//last step:: reset next
				next = sa_hash[next].next;
			}
		}
	}
	else//reverse
	{
		uint8_t  MOVE = ((S_A_KEMR_L << 1) - 2);
		//build hash table for string 1
		uint8_t * c_t_str = t_str + t_len - S_A_KEMR_L - 4;//skip first three bit
		uint64_t kmer = bit2_preKmer_init(c_t_str, S_A_KEMR_L);
		//sparse alignment
		for(int i = 4; i < t_kmer_num; i++, c_t_str--)
		{
			kmer = bit2_preKmerMOVE(c_t_str, kmer);
			if((i & 0x03) != 0)
				continue;
			uint32_t next = sa_hash[kmer &KEY_MASK].next;
			while(next != 0)
			{
				if(sa_hash[next].kmer == kmer)//filter 1, kmer match
				{
					//get read pos
					uint32_t q_pos = sa_hash[next].pos;
					if(q_pos >= q_bg && q_bg <= q_ed)//filter 2
					{
						//search forward 4 bit
						int forward_len = MEM_search(q_str + q_pos + S_A_KEMR_L, c_t_str + S_A_KEMR_L, true, 4);
						//forward_len == 4:skip
						//forward_len < 4; forward_len = forward_len;
						if(forward_len < 4)//filter 3
						{
							//search forward
							uint32_t max_search = q_pos;//query restriction: 0 ~ q_pos
							max_search = MIN(max_search, c_t_str - t_str) + OVER_SEARCH_M2;//target restriction:
							int back_len = MEM_search(q_str + q_pos - 1, c_t_str - 1, false, max_search);
							//store if total score >= 4
							int total_len = back_len + forward_len + 1;
							if(total_len >= 4)//store the index of end edge
							{
								spd_match * spd_p;
								kv_pushp_2(spd_match, sms, spd_p);
								spd_p->len = total_len;
								spd_p->q_pos = q_pos - back_len;//point to the left most point
								spd_p->t_pos = c_t_str - t_str - back_len + t_st;// point to the left most point
							}
						}
					}
				}
				//last step:: reset next
				next = sa_hash[next].next;
			}
		}
	}

}

int sdp_middle_M2(Anchor * c_a, DA_IDX * idx, spd_match_set* sms, uint8_t *q_str, sparse_align_HASH *sa_hash, int key_len)
{
	int score = 10000;//basic score: 100
	uint64_t t_offset = idx->r_i_v.a[c_a->ref_ID].seq_offset;//reference global offset
	Anchor * pre_a = NULL;
	while(c_a != NULL)
	{
		//add score
		//score += c_a->a_m.mtch_len - S_A_KEMR_L + 1;
		//debug_print_score(c_a->index_in_read, c_a->a_m.mtch_len - S_A_KEMR_L + 1, 0);
		//add middle score
		pre_a = c_a->chain_anchor_pre;//for each anchor
		if(pre_a != NULL)
		{
			int pre_mch = pre_a->a_m.mtch_len;
			int pre_refoffset = pre_a->ref_offset - 3;
			int total_ref_len = c_a->ref_offset - (pre_refoffset + pre_mch) + 3;

			sms->n = 0;
			spd_match * spd_p;
			kv_pushp_2(spd_match, sms, spd_p);
			spd_p->score = score;
			spd_p->q_pos = pre_a->index_in_read;
			spd_p->t_pos = pre_refoffset;
			spd_p->len   = - S_A_KEMR_L;
			if(total_ref_len > 12)
			{
				//get reference
				uint8_t ref[2000];//get ref
				xassert(total_ref_len < 2000, "");//todo::
				uint64_t ref_offset = pre_refoffset + t_offset + pre_mch;
				get_ref(idx->ref_bin.a, ref, ref_offset, total_ref_len, true);
				sdp_match(pre_a->index_in_read + pre_mch - 4, c_a->index_in_read + 3, q_str, ref, total_ref_len, key_len, sa_hash, sms, pre_refoffset + pre_mch, true);//todo: store the first node
			}
			//push end node
			kv_pushp_2(spd_match, sms, spd_p);
			spd_p->q_pos = c_a->index_in_read;
			spd_p->t_pos = c_a->ref_offset;
			spd_p->len   = c_a->a_m.mtch_len - S_A_KEMR_L + 1;

			if(sms->n > 1)//with data
			{
				spd_match *spd_bg = sms->a + 1,  *spd_ed = sms->a + sms->n, *c_spd = spd_bg;
				for(; c_spd < spd_ed; c_spd++)
				{
					int max_score = c_spd->len;
					uint32_t max_q = c_spd->q_pos - S_A_KEMR_L + 3;
					uint32_t max_t = c_spd->t_pos - S_A_KEMR_L + 3;

					for(spd_match * c_pre_sms = c_spd - 1; c_pre_sms >= sms->a; c_pre_sms--)
					{
						//condition1~5 that can`t be chained together, condition 4 will break loop:
						if(c_pre_sms->q_pos + c_pre_sms->len > max_q)//overlap query
							continue;
						if(c_pre_sms->t_pos + c_pre_sms->len > max_t)//overlap target
							continue;
						int indel = c_pre_sms->q_pos - c_pre_sms->t_pos - (max_q - max_t);
						int ABS_indel = ABS(indel);
						if(ABS_indel > 200)
							continue;
						int new_score = c_pre_sms->score + c_spd->len - (ABS_indel >> 3);
						max_score = MAX(max_score, new_score);
					}
					score = MAX(max_score, score);
					c_spd->score = max_score;
				}
			}
		}
		else	//final one
		{
			score += c_a->a_m.mtch_len - S_A_KEMR_L + 1;
		}
		c_a = pre_a;
	}
	return score - 10000;//delete basic score: 100
}

int sdp_right_M2(DA_IDX * idx, spd_match_set* sms, uint8_t *q_str, sparse_align_HASH *sa_hash, int key_len,
		chain_item * c_st, int chain_ID, uint32_t l_read, seed_con_hash *sc_hash, int score_ori)
{
	//INIT
	score_ori += 10000;
	int total_max_score = score_ori;//global score
	int max_sms_id = 0;
	chain_item * c_h = c_st + chain_ID;//current hit
	chain_item * combined_chian;//combine chain
	sms->n = 0;
	uint8_t ref[1000];//get ref

	//store the first node
	spd_match * spd_p;
	kv_pushp_2(spd_match, sms, spd_p);
	spd_p->score = score_ori;
	spd_p->q_pos = c_h->q_ed;
	spd_p->t_pos = c_h->t_ed;
	spd_p->len   = - S_A_KEMR_L;
	uint32_t current_sms = 1;

	//t offset
	uint64_t t_offset_global = idx->r_i_v.a[c_h->ref_ID].seq_offset;//global offset for "begin" of target
	uint64_t t_length = idx->r_i_v.a[c_h->ref_ID].seq_l;
	uint32_t c_t_offset = c_h->t_ed - 3;//local offset for target
	int last_search = false;
	while(1)
	{
		//check and reload when no node in list, load 600 bp ref
		if(sms->n == current_sms)
		{
			//check:
			uint32_t next_step = t_length - c_t_offset;
			if(next_step < MIN_SCORE_MEM)
				break;
			//reload
			uint32_t max_search_ref;
			if(l_read - c_h->q_ed < 600)////restrict max search when read not very long
			{
				if(last_search == true)//ending searching when read not enough, for 2 times
					break;
				last_search = true;
				max_search_ref = l_read - c_h->q_ed + 60;
			}
			else
				max_search_ref = t_length - c_t_offset;
			max_search_ref = MIN(600, max_search_ref);//load 600 bp for most
			get_ref(idx->ref_bin.a, ref, c_t_offset + t_offset_global, max_search_ref + OVER_SEARCH_M2, true);
			uint32_t max_search_read = MIN(c_h->q_ed + 4000, l_read);//from q_end to at most 4000
			max_search_ref = MIN(max_search_ref, max_search_read - c_h->q_ed + 60);//restrict max search ref when reach end of reads
			//debug_print_ACGTstring(q_str +  c_h->q_ed, max_search_read - c_h->q_ed);//todo
			//debug_print_ACGTstring(ref, max_search_ref);//todo

			//all read
			//debug_print_ACGTstring(q_str, l_read);
			//all ref
			//uint8_t all_ref[30001];
			//get_ref(idx->ref_bin.a, all_ref, t_offset_global, 30000, true);//from c_t_offset - max_search_ref + 1 to c_t_offset
			//debug_print_ACGTstring(all_ref, 30000);


			sdp_match(c_h->q_ed, max_search_read, q_str, ref, max_search_ref, key_len, sa_hash, sms ,c_t_offset, true);
			c_t_offset += max_search_ref - S_A_KEMR_L - 3;

			//when load no new node
			if(sms->n == current_sms)
				break;
			//when the first new node`s waiting length over 600 bp, break
			if(sms->a[current_sms].t_pos > sms->a[max_sms_id].t_pos + 1000 )
				break;
		}

		//get one new node
		spd_match *c_sms = sms->a + current_sms++;
		//SDP, get max score for that new node
		int max_score = c_sms->len;
		uint32_t max_pre_q = c_sms->q_pos - S_A_KEMR_L + 3;
		uint32_t max_pre_t = c_sms->t_pos - S_A_KEMR_L + 3;
		spd_match *c_sms_ed = sms->a, *c_pre_sms = sms->a + current_sms - 2;
		for(; c_pre_sms >= c_sms_ed; c_pre_sms--)
		{
			//condition1~5 that can`t be chained together, condition 4 will break loop:
			if(c_pre_sms->q_pos + c_pre_sms->len > max_pre_q)//overlap query
				continue;
			if(c_pre_sms->t_pos + c_pre_sms->len > max_pre_t)//overlap target
				continue;
			if(c_pre_sms->q_pos + 600 < max_pre_q)//query waiting length too long
				break;
			if(c_pre_sms->t_pos + 600 < max_pre_t)//target waiting length too long
				break;//should be break
			int indel = c_pre_sms->q_pos - c_pre_sms->t_pos - (max_pre_q - max_pre_t);
			int ABS_indel = ABS(indel);
			if(ABS_indel > 200)
				continue;
			int new_score = c_pre_sms->score + c_sms->len - (ABS_indel >> 3);
			max_score = MAX(max_score, new_score);
		}
		c_sms->score = max_score;

		//try to combine
		if(c_sms->len >= 8 &&//16 bp to combine
			combine_chain(c_st, chain_ID, sc_hash, c_sms->t_pos - c_sms->q_pos, false, c_sms->q_pos, &combined_chian) == true)//search sc_hash
		{
			//using new max score, not global max score to make sure anchor nodes be chained together, "- c_sms->len" to avoid repeat calculation
			total_max_score = MAX(score_ori, max_score) - c_sms->len + sdp_middle_M2(combined_chian->chain_anchor_cur, idx, sms, q_str, sa_hash, key_len);
			//rebegin and INIT, new a new sms node
			score_ori = total_max_score;
			max_sms_id = 0;
			sms->n = 0;
			kv_pushp_2(spd_match, sms, spd_p);
			spd_p->score = total_max_score;
			spd_p->q_pos = c_h->q_ed;
			spd_p->t_pos = c_h->t_ed;
			spd_p->len   = - S_A_KEMR_L;
			current_sms = 1;
			c_t_offset = c_h->t_ed;//local offset for target
			//then continue
			continue;
		}

		//get new best score
		if(total_max_score < max_score)
		{
			total_max_score = max_score;
			max_sms_id = current_sms - 1;
		}

		//break when waiting length over 600 bp
		if(c_sms->t_pos > sms->a[max_sms_id].t_pos + 1000)
			break;
	}
	//end of extension
	c_h->q_ed = sms->a[max_sms_id].q_pos + sms->a[max_sms_id].len + S_A_KEMR_L;
	c_h->t_ed = sms->a[max_sms_id].t_pos + sms->a[max_sms_id].len + S_A_KEMR_L;
	return total_max_score - 10000;
}

int sdp_left_M2(DA_IDX * idx, spd_match_set* sms, uint8_t *q_str, sparse_align_HASH *sa_hash, int key_len,
		chain_item * c_st, int chain_ID, uint32_t l_read, seed_con_hash *sc_hash, int score_ori)
{
	//INIT
	score_ori += 10000;
	int total_max_score = score_ori;//global score
	int max_sms_id = 0;
	chain_item * c_h = c_st + chain_ID, * combined_chian;//current hit and combine chain
	sms->n = 0;
	uint8_t ref[1000];//get ref

	//store the first node
	spd_match * spd_p;
	kv_pushp_2(spd_match, sms, spd_p);
	spd_p->score = score_ori;
	spd_p->q_pos = c_h->q_st;
	spd_p->t_pos = c_h->t_st;
	uint32_t current_sms = 1;

	//t offset
	uint64_t t_offset_global = idx->r_i_v.a[c_h->ref_ID].seq_offset;//global offset for "begin" of target
	uint32_t c_t_offset = c_h->t_st + 3;//local offset for target
	int last_search = false;
	while(1)
	{
		//check and reload when no node in list, load 600 bp ref
		if(sms->n == current_sms)
		{
			//check:
			uint32_t next_step = c_t_offset;
			if(next_step < MIN_SCORE_MEM)
				break;
			//reload
			uint32_t max_search_ref;
			if(c_h->q_st < 600)////restrict max search when read not very long
			{
				if(last_search == true)//ending searching when read not enough, for 2 times
					break;
				last_search = true;
				max_search_ref = c_h->q_st + 60;
			}
			else
				max_search_ref = c_t_offset;
			max_search_ref = MIN(600, max_search_ref);//load 600 bp for most

			if(t_offset_global == 0 && c_t_offset < OVER_SEARCH_M2 + max_search_ref)//bug
				get_ref(idx->ref_bin.a, ref, c_t_offset + t_offset_global - max_search_ref, max_search_ref, true);//from c_t_offset - max_search_ref + 1 to c_t_offset
			else
				get_ref(idx->ref_bin.a, ref, c_t_offset + t_offset_global - max_search_ref - OVER_SEARCH_M2, max_search_ref + OVER_SEARCH_M2, true);//from c_t_offset - max_search_ref + 1 to c_t_offset
			int min_search_read = MAX((int)c_h->q_st - 4000, 0);//from q_end to at most 4000, must be int

			//debug_print_ACGTstring(q_str +  min_search_read, c_h->q_st - min_search_read);
			//debug_print_ACGTstring(ref, max_search_ref);

			//all read
			//debug_print_ACGTstring(q_str, l_read);
			//all ref
			//uint8_t all_ref[30001];
			//get_ref(idx->ref_bin.a, all_ref, t_offset_global, 30000, true);//from c_t_offset - max_search_ref + 1 to c_t_offset
			//debug_print_ACGTstring(all_ref, 30000);

			sdp_match(min_search_read, c_h->q_st, q_str, ref + OVER_SEARCH_M2, max_search_ref, key_len, sa_hash, sms ,c_t_offset - max_search_ref, false);
			c_t_offset = c_t_offset - max_search_ref + S_A_KEMR_L + 3;//todo

			//when load no new node
			if(sms->n == current_sms)
				break;
			//when the first new node`s waiting length over 600 bp, break
			if(sms->a[current_sms].t_pos + 1000 < sms->a[max_sms_id].t_pos)
				break;
		}
		//get one new node
		spd_match *c_sms = sms->a + current_sms++;
		//SDP, get max score for that new node
		int max_score = c_sms->len;
		uint32_t min_pre_q = c_sms->q_pos + S_A_KEMR_L + c_sms->len - 3;//next <= pre - pre_len - 9 ==> next + 9 <= pre - pre_len
		uint32_t min_pre_t = c_sms->t_pos + S_A_KEMR_L + c_sms->len - 3;
		spd_match *c_sms_ed = sms->a, *c_pre_sms = sms->a + current_sms - 2;
		for(; c_pre_sms >= c_sms_ed; c_pre_sms--)
		{
			//condition1~5 that can`t be chained together, condition 4 will break loop:
			if(c_pre_sms->q_pos < min_pre_q)//overlap query
				continue;
			if(c_pre_sms->t_pos < min_pre_t)//overlap target
				continue;
			if(min_pre_q + 600 < c_pre_sms->q_pos)//query waiting length too long, next + 600 < pre
				break;
			if(min_pre_t + 600 < c_pre_sms->t_pos)//target waiting length too long
				break;//should be break
			int indel = c_pre_sms->q_pos - c_pre_sms->t_pos - (min_pre_q - min_pre_t);
			int ABS_indel = ABS(indel);
			if(ABS_indel > 200)
				continue;
			int new_score = c_pre_sms->score + c_sms->len - (ABS_indel >> 3);
			max_score = MAX(max_score, new_score);
		}
		c_sms->score = max_score;

		//try to combine
		if(c_sms->len >= 8 &&//16 bp to combine
			combine_chain(c_st, chain_ID, sc_hash, c_sms->t_pos - c_sms->q_pos, true, c_sms->q_pos + c_sms->len, &combined_chian) == true)//search sc_hash
		{
			//using new max score, not global max score to make sure anchor nodes be chained together, "- c_sms->len" to avoid repeat calculation
			total_max_score = MAX(score_ori, max_score) - c_sms->len + sdp_middle_M2(combined_chian->chain_anchor_cur, idx, sms, q_str, sa_hash, key_len);
			//rebegin and INIT, new a new sms node
			score_ori = total_max_score;
			max_sms_id = 0;
			sms->n = 0;
			kv_pushp_2(spd_match, sms, spd_p);
			spd_p->score = total_max_score;
			spd_p->q_pos = c_h->q_st;
			spd_p->t_pos = c_h->t_st;
			current_sms = 1;
			c_t_offset = c_h->t_st;//local offset for target
			//then continue
			continue;
		}
		//get new best score
		if(total_max_score < max_score)
		{
			total_max_score = max_score;
			max_sms_id = current_sms - 1;
		}
		//break when waiting length over 1000 bp(over 1000 bp no increase in score)
		if(c_sms->t_pos + 1000 < sms->a[max_sms_id].t_pos)
			break;
	}
	//end of extension
	c_h->q_st = sms->a[max_sms_id].q_pos;
	c_h->t_st = sms->a[max_sms_id].t_pos;
	return total_max_score - 10000;
}

void get_score_M2(SEARCH_DIR *search_dir, Classify_buff_pool *buff,
		DA_IDX * idx, uint32_t l_read, cly_r *results, seed_con_hash *sc_hash)
{
	//BUILD space hash index for read

	//sa_hash = malloc(sizeof(sparse_align_HASH)*0x100000)//at most 2M blocks, 0.25M for hash table; 1.75M for max read length
	int key_len =build_hash_table_M2(search_dir, &(results->hit), l_read, buff->sa_hash);
	//Part 3:read hash
	chain_item *st_hit = results->hit.a;
	for(int i = 0; i < results->hit.n; i++)
	{
		if(st_hit[i].sum_score == 0)//skip result already combined by other hit
			continue;
		//fprintf(stderr,"\n\n\n%s %d\n\n\n\n", idx->r_i_v.a[st_hit[i].ref_ID].ref_name, st_hit[i].t_st);//DEBUG+++
		//if(504221 == st_hit[i].t_st)
		//	fprintf(stderr, " ");
		//middle
		SEARCH_DIR *c_search_dir = ((search_dir->direction == st_hit[i].direction)?0:1) + search_dir;
		sparse_align_HASH * sa_hash = (st_hit[i].direction==FORWARD)?buff->sa_hash[0]:buff->sa_hash[1];
		//todo:
		spd_match_set* sms = &(buff->sms);//todo:sms
		int score = sdp_middle_M2(st_hit[i].chain_anchor_cur,idx, sms, c_search_dir->bin_read, sa_hash, key_len);
		score = sdp_right_M2(idx, sms, c_search_dir->bin_read, sa_hash, key_len, st_hit, i, l_read, sc_hash, score);
		score = sdp_left_M2(idx, sms, c_search_dir->bin_read, sa_hash, key_len, st_hit, i, l_read, sc_hash, score);
		st_hit[i].sum_score = score;//score -> sort score
	}
}

///------------------------------------END Get SCORE M2--------------------------------------------------

static int inline chain_cmp_by_pos(const void * a_, const void * b_)
{
	chain_item *a = (chain_item *)a_, *b = (chain_item *)b_;
	if(a->ref_ID > b->ref_ID)
		return 1;
	if(a->ref_ID < b->ref_ID)
		return -1;
	//when ID equal
	if(a->t_st > b->t_st)
		return 1;
	if(a->t_st < b->t_st)
		return -1;
	if(a->sum_score < b->sum_score)
		return 1;
	if(a->sum_score > b->sum_score)
		return -1;
	return 0;
}
#define GET_SCORE_M1

//#define FILTER_MIN_LEN 170
#define FILTER_MIN_SCORE_SHROT_3G_READ 30
#define FILTER_MIN_SCORE_2G_READ 26//35 when score over this value, the results will not be filtered whatever the length; used for 2nd Generation data(Illumina)
//#define FILTER_MIN_SCORE_LEVEL_3 (FILTER_MIN_SCORE + 10)//when score over this value, the results will not be filtered whatever the length; used for long read data
//int result_filter(chain_item *c)
//{
//	if(c->q_ed - c->q_st < FILTER_MIN_LEN || c->sum_score < FILTER_MIN_SCORE)
//		return 0;
//	return 1;
//}
void delete_small_score_rst(
		DA_IDX * idx,
		cly_r *results,
		SEARCH_DIR *search_dir,
		Classify_buff_pool *buff)//chain and polish
{
	if(results->hit.n == 0)
		return;
	if(results->hit.n > 150)//delete result with only one anchor
	{
		int rst_num = 50;//delete result without top anchor
		for(;rst_num < results->hit.n && results->hit.a[rst_num].anchor_number > 1;rst_num++);
		results->hit.n = rst_num;
	}
	if(results->hit.n > 10000)
	{
		fprintf(stderr, "[%s], too many result, [%ld]\n", results->read->name.s, results->hit.n);
		results->hit.n = MIN(10000, results->hit.n);
	}
	uint32_t l_read = results->read->seq.l;

	//Part 2: //seed combine hash	//seed hash, build//store seed hash
	//seed hash, malloc
	uint32_t n_sc_hash = 256 + (results->hit.n << 1);
	BUFF_REALLOC(buff->sc_hash, buff->m_sc_hash, n_sc_hash);
	seed_con_hash *sc_hash = buff->sc_hash;
	sc_hash_idx(sc_hash, n_sc_hash, &(results->hit));

	//Part 3:read hash
	//get_score_M1(search_dir, buff, idx, l_read, results, sc_hash);
	get_score_M2(search_dir, buff, idx, l_read, results, sc_hash);

	//filter, set score to 0
	chain_item * st_c = results->hit.a, *ed_c = st_c + results->hit.n, *c_c;
	//chain seeds together
	//step1: sort
	if(results->hit.n > 1)
		qsort(results->hit.a, results->hit.n, sizeof(chain_item), chain_cmp_by_pos);
	for(c_c = st_c;c_c < ed_c - 1; c_c++)//try to combine nearby result and delete duplication results
	{
		if(c_c->sum_score == 0)
			continue;
		chain_item *next_c = c_c + 1;
		for(;next_c < ed_c; next_c++)
		{
			if(c_c->ref_ID == next_c->ref_ID)
			{
				if(c_c->direction != next_c->direction)
					continue;
				if(next_c->sum_score == 0)
					continue;
				//when totally overlap
				if(next_c->t_st - 5 < c_c->t_st && next_c->q_st - 5 < c_c->q_st && next_c->sum_score - 5 < c_c->sum_score)
				{
					next_c->sum_score = 0;
					next_c->q_ed = next_c->q_st;
					next_c->t_ed = next_c->t_st;
					continue;
				}
				//try to combine
				int dis_t = next_c->t_st - c_c->t_ed;
				int dis_q = next_c->q_st - c_c->q_ed;
				int dis_t_q = ABS(dis_t - dis_q);
				if((dis_t > -20 && dis_t < 1000 && dis_q > -20 && dis_q < 1000)//distance within threshold
						&& dis_t_q < 200)//indel size < 250bp
				{
					c_c->t_ed = MAX(c_c->t_ed, next_c->t_ed);
					c_c->q_ed = MAX(c_c->q_ed, next_c->q_ed);
					c_c->sum_score += next_c->sum_score;
					next_c->sum_score = 0;
					next_c->q_ed = next_c->q_st;
					next_c->t_ed = next_c->t_st;
				}
			}
			else
				break;
		}
	}
	buff->max_read_l = MAX(buff->max_read_l, l_read);

	if(buff->max_read_l < 510)// 2G reads
	{
		for(c_c = st_c;c_c < ed_c; c_c++)
		{
			int score = c_c->sum_score + ((c_c->q_ed - c_c->q_st) >> 5);
			if(score < FILTER_MIN_SCORE_2G_READ)
				c_c->sum_score = 0;
		}
	}
	else if(l_read < 310)//short 3G reads
	{
		for(c_c = st_c;c_c < ed_c; c_c++)
		{
			int score = c_c->sum_score + ((c_c->q_ed - c_c->q_st) >> 5);
			if(score < FILTER_MIN_SCORE_SHROT_3G_READ)
				c_c->sum_score = 0;
		}
	}
	else//long 3G reads
	{
		for(c_c = st_c;c_c < ed_c; c_c++)
		{
			int score = c_c->sum_score + ((c_c->q_ed - c_c->q_st) >> 5);
			if(score < (idx->filter_min_score_LV3) && (c_c->q_ed - c_c->q_st < idx->filter_min_length || score < idx->filter_min_score))
				c_c->sum_score = 0;
		}
	}
	if(results->hit.n > 1)
		qsort(results->hit.a, results->hit.n, sizeof(chain_item), chain_cmp_by_MEM_score);
	for(c_c = st_c;c_c < ed_c; c_c++)
		if(c_c->sum_score == 0)
			break;
	results->hit.n = c_c - st_c;
}

void detect_primary(chain_item *hit, uint32_t n_hit, uint32_t read_len)
{
	if(n_hit == 0)
		return;
	//primary list
	int primary_v[300];
	uint8_t primary_v_idx[300];
	int n_primary_v = 1;

	//set PRIMARY
	hit->pri_index = primary_v_idx[0] = 0;
	primary_v[0] = 0;
	hit->primary = PRIMARY;
	//
	chain_item *ed_hit = hit + n_hit;
	for(chain_item *c_hit = hit + 1; c_hit < ed_hit; c_hit++)
	{
		bool overlap = false;
		for(int i = 0; i < n_primary_v; i++)
		{
			int primary_st, primary_ed;
			if(hit[primary_v[i]].direction == c_hit->direction)
			{
				primary_st = hit[primary_v[i]].q_st;
				primary_ed = hit[primary_v[i]].q_ed;
			}
			else
			{
				primary_st = read_len - hit[primary_v[i]].q_ed;
				primary_ed = read_len - hit[primary_v[i]].q_st;
			}
			uint32_t overlap_st = MAX(c_hit->q_st, primary_st);
			uint32_t overlap_ed = MIN(c_hit->q_ed, primary_ed);
			if((overlap_st < overlap_ed) &&
				(((overlap_ed - overlap_st) << 1) >= (c_hit->q_ed - c_hit->q_st)))//50% overlap
					overlap = true;
			if(overlap)
			{
				c_hit->primary = SECONDARY;
				c_hit->pri_index = ++primary_v_idx[i];
				//length and scores (nearly) equal with PRIMARY
				if(c_hit->sum_score + 5 > hit[primary_v[i]].sum_score)
					c_hit->pri_index = 1;
				if(primary_v_idx[i] == 255)
					primary_v_idx[i] = 254;
				break;
			}
		}
		if(overlap == false)
		{
			c_hit->primary = SUPPLYMENTARY;
			c_hit->pri_index = primary_v_idx[n_primary_v] = 0;
			primary_v[n_primary_v++] = c_hit - hit;
		}
	}
}

/**
 * 	use buff to store bin seq; get from buff pool to avoid malloc
 */
#define MIN_READ_LEN 40//ignore reads that shorter than 40 bp
void classify_seq(kseq_t *read, DA_IDX* idx, cly_r *results, Classify_buff_pool *buff)
{
	//head
	SEARCH_DIR search_dir[2];
	uint32_t read_len = read->seq.l;
	results->anchor_v.n = 0;
	results->read = read;
	results->fast_classify = true;
	results->hit.n = 0;
	//fprintf(stderr, "%s\n", read->name.s);
	/** DEBUG code
	//get read ID
	fprintf(stderr, "%s\n", read->name.s);
	char read_name_dbg[100];
	strcpy(read_name_dbg, read->name.s);
	strcat(read_name_dbg, "!");
	char *tokens;
	tokens = strtok(read_name_dbg, ".");
	tokens = strtok(NULL, "!");
	uint32_t read_ID = strtoul(tokens, NULL, 10);
	if(read_ID > 263)///
		printf(" ");
	//else
	//	return;
	//*/
	if(read_len < MIN_READ_LEN)
		return;
	double t_island = 0;
	FUNC_GET_TIME(getIsland(read, buff, &(idx->ek), search_dir), t_island, true);//get search dir
	bool both_direction = false;//is the scores in both search directions similar?
#ifndef CONSIDER_BOTH_ORIENTATION
	both_direction = ((search_dir[0].total_score - search_dir[1].total_score) <= (search_dir[0].total_score >> 3))?true:false;
#endif
	//FAST MODE-loop1
	fast_classify(idx, search_dir, read_len, results);//fast mode for error rate < 0.2
	if(both_direction)
		fast_classify(idx, search_dir + 1, read_len, results);//fast mode for error rate < 0.2
	resolve_tree(results);
	//slow mode: loop2
	int run_slow_mode = false;
	if(results->hit.n <= 0)//when no results
		run_slow_mode = true;
	else if(results->hit.a[0].anchor_number < 5)//when not sufficient anchors
	{
		run_slow_mode = true;
		//run in 2nd Generation reads mode
		if(read_len <= 300 && results->hit.a[0].sum_score > 200)//when running high quality short reads
			run_slow_mode = false;
	}
	//SLOW MODE-loop2
	if(run_slow_mode)
	{
		results->anchor_v.n = 0; //clear chor list
		slow_classify(idx, search_dir, read_len,  results);
		resolve_tree(results); //SLOW MODE for error rate over an0.2
#ifndef CONSIDER_BOTH_ORIENTATION
		////SLOW MODE-loop3
		if(both_direction || results->hit.n <= 0 || results->hit.a[0].anchor_number < 5)
		{
			slow_classify(idx, search_dir + 1, read_len, results);
			resolve_tree(results);
		}
#endif
	}
	//delete result with too small score
	delete_small_score_rst(idx, results, search_dir, buff);
	//detect PRIMARY/SUPPLEMENTARY/SECONDARY
	detect_primary(results->hit.a, results->hit.n, read_len);
}
