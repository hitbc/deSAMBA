#include "idx.h"
#include <string.h>
#include "lib/utils.h"
#include "bwt.h"
#include "lib/desc.h"
#include <pthread.h>
#include <getopt.h>
#define B_KMER 31
uint8_t Bit[] =
{
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

uint8_t bin_Bit[] =
{//all other char are used as 'A'
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0/*A*/, 0, 1/*C*/,  0, 0, 0, 2/*G*/,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  3/*T*/, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0/*A*/, 0, 1/*C*/,  0, 0, 0, 2/*G*/,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  3/*T*/, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};

char Chars[] = {'A','C','G','T','#','$'};

uint8_t edgesNum[16] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
#define isStart(info) (((info) >> 9) & 1)
#define isEnd(info)	  (((info) >> 8) & 1)
#define setStart(marker) ((marker) |= (1 << 9))
#define setEnd(marker)   ((marker) |= (1 << 8))
#define setInEdge(info, edge)  ((info) |= (1 << ((edge) + 4)))
#define setOutEdge(info, edge) ((info) |= (1 << (edge)))
#define outEdgesNum(info) (edgesNum[((info) & 0xf)])
#define inEdgesNum(info)  (edgesNum[(((info) >> 4) & 0xf)])

#define b_search_uint64(loc, kmersValue, lowerBound, upperBound, key)	\
do{																		\
	uint64_t l_B = (lowerBound), u_B = (upperBound);					\
	while (l_B <= u_B) {												\
		loc = l_B + ((u_B - l_B) >> 1);									\
		if 		(kmersValue[loc] < key) 	l_B = loc + 1;				\
		else if (kmersValue[loc] > key) 	u_B = loc - 1;				\
		else								break;						\
	}																	\
	xassert(l_B <= u_B, "Bin search wrong, loc NOT found");				\
}while(0)

//upperBound means (number_of_list - 1); lowerBound is 0;
static inline uint64_t binSearch(uint64_t *kmersValue,  uint64_t lowerBound, uint64_t upperBound, uint64_t key)
{
	uint64_t mid;
	while (lowerBound <= upperBound) {
		mid = lowerBound + ((upperBound - lowerBound) >> 1);
		if 		(kmersValue[mid] < key) 	lowerBound = mid + 1;
		else if (kmersValue[mid] > key) 	upperBound = mid - 1;
		else								return mid;
	}
	xassert(false, "Bin search wrong, kmer NOT found");
	return -1;
}

//fore part of kmer will be use to limit the range 
void print(char *str, int len)
{
	for(uint8_t i = 0; i < len; ++i)
		fprintf(stderr,"%c",str[i]);
	fprintf(stderr,"\n");
}

uint64_t * getCounter(uint64_t  * kmerValue, uint64_t kmerNum)
{
	uint64_t counterSize = (1<<(L_PRE_IDX<<1)) + 1;//2^26 + 1 = 64M
	uint16_t move = (B_KMER - L_PRE_IDX) << 1;//(31 - 13)*2 = 36
	uint64_t * counter = (uint64_t*) xcalloc(counterSize, sizeof(uint64_t));
	//initiate ranges
	uint64_t sum = 0, temp = 0;
	for(uint64_t i=0; i<kmerNum; ++i)
		++counter[kmerValue[i]>>move];
	for(uint64_t i=0; i<counterSize;++i)
	{
		temp = counter[i];
		counter[i] = sum; sum = sum + temp;
	}
	return counter;
}

typedef struct
{
	const char *refPath;
	uint64_t * kmerValue;
	uint16_t *kmerInfo;
	uint64_V heads;
	uint64_V tails;
	uint64_t *counter;
	uint64_t len_ref;
	uint64_t MASK;
	uint64_t suffix;
}BUILD_DEB_DATA;

static void *build_deb(void *data)
{ //the path indicate a file containing all sorted k+1mers
	BUILD_DEB_DATA *b = (BUILD_DEB_DATA*) data;
	uint64_t suffix = b->suffix;
	uint64_V *heads = &(b->heads);
	uint64_V *tails = &(b->tails);
	uint64_t *kmerValue = b->kmerValue;
	uint16_t *kmerInfo = b->kmerInfo;
	uint64_t *counter = b->counter;
	uint64_t MASK = b->MASK;

	gzFile fp_z = xzopen(b->refPath,"r");
	kseq_t *trunk = kseq_init(fp_z);
	uint16_t move = (B_KMER - L_PRE_IDX) << 1;//(31 - 13)*2 = 36
	uint64_t mask = ~((uint64_t)0x3 <<((B_KMER-1)<<1));
#ifdef CONSIDER_BOTH_ORIENTATION
	uint64_t mask_rev[] = {(uint64_t)0x3 << ((B_KMER - 1)<<1),  (uint64_t)0x2 << ((B_KMER - 1)<<1), (uint64_t)0x1 << ((B_KMER - 1)<<1), 0};
#endif
	while (kseq_read(trunk) >= 0)
	{//read ref seq
		size_t seq_l = trunk->seq.l;
		b->len_ref += seq_l;
		char * seq_s = trunk->seq.s;
		//DBG code:
		if(suffix == 0)
			fprintf(stderr,"%s %ld\r", trunk->name.s, seq_l);
		for(uint64_t i = 0; i < seq_l; ++i)
		{
			if (Bit[(int8_t)seq_s[i]] != 4)
			{
				uint64_t start = i;
				while(Bit[(int8_t)seq_s[++i]]!=  4 && i < seq_l);
				if (i - start > B_KMER)
				{
					///start
					uint64_t key = char2Kmer(seq_s+start, B_KMER, Bit);
					if((key & MASK) == suffix)
					{
						uint64_t loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
						setOutEdge( kmerInfo[loc], Bit[(int8_t)seq_s[start + B_KMER]]);
						kv_push(uint64_t, *heads, loc);
					}
#ifdef CONSIDER_BOTH_ORIENTATION
					uint64_t keyR = char2KmerRC(seq_s+start, B_KMER, Bit);
					if((key & MASK) == suffix)
					{
						uint64_t loc = binSearch(kmerValue, counter[keyR>>move],counter[(keyR>>move)+1]-1, keyR);
						setInEdge(kmerInfo[loc],Bit[(int8_t)seq_s[start + B_KMER]]^0x3);
						kv_push(uint64_t, *tails, loc);
					}
#endif
					///middle
					for (uint64_t j = start + 1; j < i - B_KMER; ++j)
					{
						key = ((key & mask) << 2)| Bit[(int8_t)seq_s[j+B_KMER-1]];
						if((key & MASK) == suffix)
						{
							uint64_t loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
							setInEdge(kmerInfo[loc], Bit[(int8_t)seq_s[j-1]]);
							setOutEdge(kmerInfo[loc], Bit[(int8_t)seq_s[j+B_KMER]]);
						}
#ifdef CONSIDER_BOTH_ORIENTATION
						keyR = (keyR >> 2) | mask_rev[Bit[(int8_t)seq_s[j + B_KMER - 1]]];
						if((key & MASK) == suffix)
						{
							uint64_t loc = binSearch(kmerValue, counter[keyR>>move],counter[(keyR>>move)+1]-1, keyR);
							setOutEdge(kmerInfo[loc], (Bit[(int8_t)seq_s[j-1]]^0x3));
							setInEdge (kmerInfo[loc], Bit[(int8_t)seq_s[j + B_KMER]]^0x3);
						}
#endif
					}
					///end
					key = ((key & mask) << 2)| Bit[(int8_t)seq_s[i-1]];
					if((key & MASK) == suffix)
					{
						uint64_t loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
						setInEdge(kmerInfo[loc], Bit[(int8_t)seq_s[i - B_KMER -1]]);
						kv_push(uint64_t, *tails, loc);
					}
#ifdef CONSIDER_BOTH_ORIENTATION
					keyR = (keyR >> 2) | mask_rev[Bit[(int8_t)seq_s[i - 1]]];
					if((key & MASK) == suffix)
					{
						uint64_t loc = binSearch(kmerValue, counter[keyR>>move],counter[(keyR>>move)+1]-1, keyR);
						setOutEdge(kmerInfo[loc], (Bit[(int8_t)seq_s[i - B_KMER -1]]^0x3));
						kv_push(uint64_t, *heads, loc);
					}
#endif
				}
				else if (i - start == B_KMER)
				{
					uint64_t key = char2Kmer(seq_s + start, B_KMER, Bit);
					if((key & MASK) == suffix)
					{
						uint64_t loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
						kv_push(uint64_t, *heads, loc);
						kv_push(uint64_t, *tails, loc);
					}
#ifdef CONSIDER_BOTH_ORIENTATION
					uint64_t keyR = char2KmerRC(seq_s+start, B_KMER, Bit);
					if((key & MASK) == suffix)
					{
						uint64_t loc = binSearch(kmerValue, counter[keyR>>move],counter[(keyR>>move)+1]-1, keyR);
						kv_push(uint64_t, *heads, loc);
						kv_push(uint64_t, *tails, loc);
					}
#endif
				}
			}
		}
	}
	kseq_destroy(trunk);
	gzclose(fp_z);
	return NULL;
}

#define N_T_BDEB 16
#define MSK_BDEB 0x96//10010110
uint64_t suffix_BDEB[N_T_BDEB] =
{
		0x0,  0x2,  0x4,  0x6,
		0x10, 0x12, 0x14, 0x16,
		0x80, 0x82, 0x84, 0x86,
		0x90, 0x92, 0x94, 0x96,
};
void build_deb_mt(
		const char *refPath,
		uint64_t * kmerValue,
		uint16_t *kmerInfo,
		uint64_V *heads,
		uint64_V *tails,
		uint64_t *counter,
		uint64_t *len_ref)
{
	pthread_t tid[N_T_BDEB] = {0};
	BUILD_DEB_DATA data[N_T_BDEB] = {0};
	for (int i = 0; i < N_T_BDEB; ++i)
	{
		data[i].refPath = refPath;
		data[i].kmerValue = kmerValue;
		data[i].kmerInfo = kmerInfo;
		data[i].counter = counter;
		data[i].MASK = MSK_BDEB;
		data[i].suffix = suffix_BDEB[i];
	}
	for (int i = 0; i < N_T_BDEB; ++i)
		pthread_create(&tid[i], 0, build_deb, &data[i]);
	for (int i = 0; i < N_T_BDEB; ++i)
		pthread_join(tid[i], 0);
	*len_ref = data[0].len_ref;
	//combine heads and tails of all threads
	size_t n_total = 0;
	//for head
	for (int i = 0; i < N_T_BDEB; ++i)
		n_total += data[i].heads.n;
	heads->n = heads->m = n_total;
	heads->a = xmalloc_t(uint64_t, heads->n);
	uint64_t * heads_p = heads->a;
	size_t idx = 0;
	for (int i = 0; i < N_T_BDEB; ++i)
	{
		uint64_t *data_i = data[i].heads.a;
		for(size_t j = 0; j < data[i].heads.n; j++)
			heads_p[idx++] = data_i[j];
		free(data[i].heads.a);
	}
	//for tail
	n_total = 0;
	for (int i = 0; i < N_T_BDEB; ++i)
		n_total += data[i].tails.n;
	tails->n = tails->m = n_total;
	tails->a = xmalloc_t(uint64_t, tails->n);
	uint64_t * tails_p = tails->a;
	idx = 0;
	for (int i = 0; i < N_T_BDEB; ++i)
	{
		uint64_t *data_i = data[i].tails.a;
		for(size_t j = 0; j < data[i].tails.n; j++)
			tails_p[idx++] = data_i[j];
		free(data[i].tails.a);
	}
}

//return a 64-bit kmer as the start of the next block
void findInsertPos(
		FILE *kmer_v_fp,
		uint64_V *kmer_value,
		uint64_t key)
{
	uint64_t kmer;
	while(1)
	{
		if(fread(&kmer, 8, 1, kmer_v_fp) < 1)
			return;
		if (kmer < key)//push
			kv_push_2(uint64_t, kmer_value, kmer);
		else
		{
			fseek(kmer_v_fp, -8, SEEK_CUR);
			return;
		}
	}
	//for(uint64_t i = p_s; i < p_e; ++i)
	//	if (key <= p_kmersValue[i])
	//		return i;
	//return p_e;
}

static inline void store_hash(uint64_t *hash_index, uint64_t new_key, uint64_t *pre_key, uint64_t idx)
{
	if (new_key != *pre_key)
	{//if new a hash
		hash_index[new_key<<1] = idx;
		hash_index[(new_key<<1)+1] = idx + 1;
		*pre_key = new_key;
	} else
		++hash_index[(new_key<<1)+1];
}

#define MERGE_HASH_MOVE ((B_KMER - L_PRE_IDX) << 1)//36
void merge_kmer(
		FILE *kmer_v_fp,
		//uint64_t* kmer_v,
		uint64_t n_kmer,
		uint16_t* kmer_i,
		SP_KMER* sp_kmer,
		uint64_t n_sp_kmer,
		kstring_t *bwt_s,
		uint64_t *hash_index)
{
	uint64_t start = 0, pre = MAX_uint64_t;
	SP_KMER* bgn_sp = sp_kmer, *end_sp = sp_kmer + n_sp_kmer;
	uint64_V kmer_v_block = {0,0,0};
	for (;bgn_sp < end_sp; ++bgn_sp)
	{//for all special positon
		uint8_t sp_pos = bgn_sp->sp_pos;
		uint64_t * bgn_sp_kmer = (uint64_t*)(bgn_sp->value);
		uint64_t key = (* bgn_sp_kmer) << ((B_KMER - sp_pos)<< 1);
		//findInsertPos
		findInsertPos(kmer_v_fp, &kmer_v_block, key);//from start to end, find merge insert position
		for(uint64_t i = 0; i < kmer_v_block.n; ++i)
		{	//bwt before insert char
			kputc(Chars[(kmer_i[i + start] >> 10) & 0x7], bwt_s);//store BWT
			store_hash(hash_index, kmer_v_block.a[i] >> MERGE_HASH_MOVE, &pre, bwt_s->l - 1);//store the hash table
		}
		start += kmer_v_block.n;//reset start
		kputc(Chars[bgn_sp->last_char], bwt_s);//the insert char
		if (sp_pos >= L_PRE_IDX)
		{//store sp kmer in the end
			uint64_t new_key = (* bgn_sp_kmer) >> ((sp_pos - L_PRE_IDX) << 1);
			store_hash(hash_index, new_key, &pre, bwt_s->l - 1);
		}
		//clean kmer block
		kmer_v_block.n = 0;
	}
	findInsertPos(kmer_v_fp, &kmer_v_block, MAX_uint64_t);//from start to end, find merge insert position
	//store the last part
	for (uint64_t i = 0; i < kmer_v_block.n; ++i)
	{//store kmers left
		kputc(Chars[kmer_i[i + start]>>10], bwt_s);
		store_hash(hash_index, kmer_v_block.a[i] >> MERGE_HASH_MOVE, &pre, bwt_s->l - 1);
	}
	start += kmer_v_block.n;//reset start
	xassert(start == n_kmer, "Wrong total kmer number!\n");
}

//set s/e based on kmerInfo
int cutOffMulEdges(uint64_t *p_kmerValue, uint16_t *p_kmerInfo, uint64_t kmerNum, uint64_t *counter)
{
	uint64_t size = kmerNum;
	uint16_t move = (B_KMER - L_PRE_IDX) << 1;
	uint64_t mask = ~((uint64_t)0x3 <<((B_KMER-1)<<1));
	for (uint64_t i=0; i<size;++i)
	{
		uint8_t inCount = inEdgesNum(p_kmerInfo[i]);
		uint8_t outCount = outEdgesNum(p_kmerInfo[i]);
		//fprintf(stderr,"%lu",i);	
		if ( 1 != inCount)
		{
			setStart( p_kmerInfo[i]);
			if (inCount)
			{
				uint8_t inEdges = p_kmerInfo[i] >> 4;
				for (uint8_t j=0; j< 4; ++j,inEdges >>= 1)
				{
					if (inEdges&1)
					{
						uint64_t key = (p_kmerValue[i] >> 2)|((uint64_t)j<<((B_KMER-1)<<1));
						uint64_t loc = binSearch(p_kmerValue, counter[key>>move], counter[(key>>move)+1]-1,key );
						setEnd( p_kmerInfo[loc]);
					}
				}
			}
		}
		if ( 1 != outCount)
		{
			setEnd( p_kmerInfo[i]);
			if (outCount)
			{
				uint8_t outEdges = p_kmerInfo[i];
				for (uint8_t j=0; j< 4; ++j,outEdges >>= 1)
				{
					if (outEdges&1)
					{
						uint64_t key =((p_kmerValue[i]&mask)<<2)|j;//j is uint8_t,A mistakes 
						uint64_t loc = binSearch(p_kmerValue, counter[key>>move], counter[(key>>move)+1] - 1, key);
						setStart(p_kmerInfo[loc]);
					} 
				}	
			}
		}
	}
	return 0;
}

int handleFrstLastKmer(
		uint64_t  *p_kmerValue,
		uint16_t *p_kmerInfo,
		uint64_V *heads,
		uint64_V *tails,
		uint64_t *counter){
	uint64_t hSize = heads->n, tSize = tails->n;
	uint64_t *h = heads->a, *t = tails->a;
	uint16_t move = (B_KMER - L_PRE_IDX) << 1;

	for (uint64_t i = 0; i < hSize; ++i)
	{
		setStart(p_kmerInfo[h[i]]);
		uint8_t inCount = inEdgesNum(p_kmerInfo[h[i]]);
		if (inCount)
		{
			uint8_t inEdges = p_kmerInfo[h[i]] >> 4;
			for (uint8_t j=0; j< 4; ++j, inEdges >>= 1)
			{
				if (inEdges&1)
				{
					uint64_t key = (p_kmerValue[h[i]] >> 2)|((uint64_t)j<<((B_KMER-1)<<1));
					uint64_t loca = binSearch(p_kmerValue,counter[key>>move], counter[(key>>move)+1] - 1, key );
					setEnd(p_kmerInfo[loca]);
				} 
			}
		}
	}
	
	uint64_t mask = ~((uint64_t)0x3 <<((B_KMER-1)<<1));
	for (uint64_t i=0;i<tSize;++i)
	{
		setEnd(p_kmerInfo[t[i]]);
		uint8_t outCount = outEdgesNum(p_kmerInfo[t[i]]);
		if (outCount)
		{
			uint8_t outEdges = p_kmerInfo[t[i]];
			for (uint8_t j=0; j< 4; ++j,outEdges >>= 1)
			{
				if (outEdges&1)
				{
					uint64_t key = (p_kmerValue[t[i]] & mask) << 2 | j;
					uint64_t loca = binSearch(p_kmerValue,counter[key>>move], counter[(key>>move) + 1] - 1,key );
					setStart(p_kmerInfo[loca]);
				} 
			}
		}
	}
	return 0;
}

uint64_t countEndNum(uint16_t* p_kmerInfo,uint64_t p_kmerNum)
{
	uint64_t _count = 0;
	for (uint64_t i=0; i<p_kmerNum; ++i)
		if(isEnd(p_kmerInfo[i]))
			++_count;
	return _count;
}	

//return unitig number, set labels, it means whether or not one kmer is the end or start of a unitig
int setLabel(
		uint64_t *kmerValue,
		uint16_t *kmerInfo,
		uint64_t kmerNum,
		uint64_V *heads,
		uint64_V *tails,
		uint64_t *counter)
{//start end start&end transit
	cutOffMulEdges(kmerValue, kmerInfo, kmerNum, counter);//set s/e based on kmerInfo
	handleFrstLastKmer(kmerValue, kmerInfo, heads, tails, counter);//set s/e based on head and tail
	return countEndNum(kmerInfo, kmerNum);//get total number of end node(UNITIG)
}

void genSpKmers(uint64_t kvalue, SP_KMER* p_2k, uint8_t* p_2k_0p, uint64_t idx)
{
	uint64_t p_2k_idx = idx*(B_KMER - 1);
	uint64_t mask = (0x1ul << ((B_KMER - 1) << 1)) - 1;//000011111111...111, == MASK_60
	for(uint8_t sp_pos = B_KMER - 1; sp_pos > 0; sp_pos --, p_2k_idx++, mask >>= 2)
	{//for the pos [0~29]
		uint64_t * kmer_value = (uint64_t*)(p_2k[p_2k_idx].value);
		* kmer_value = kvalue & mask;
		p_2k[p_2k_idx].sp_pos = sp_pos;
		p_2k[p_2k_idx].last_char = (kvalue >> (sp_pos << 1)) & 0x3;
	}
	p_2k_0p[idx] = kvalue & 0x3;//for pos 30, the last pos
}

static inline void store_REF_UNITIG(
		REF_UNI * r_u,
		uint64_t ref_ID,
		uint64_t UNITIG_ID,
		uint64_t ref_offset,
		int direction)
{
	r_u->ref_ID = ref_ID;
	r_u->UNITIG_ID = UNITIG_ID;
	r_u->ref_offset = ref_offset;
	r_u->direction = direction;
}

//small to big
static int inline REF_UNITIG_cmp_by_UNITIG_ID(const void*a_,const void*b_)
{
	uint32_t UNI_ID_a = ((REF_UNI *)a_)->UNITIG_ID,	UNI_ID_b = ((REF_UNI *)b_)->UNITIG_ID;
	if(UNI_ID_a > UNI_ID_b)
		return -1;
	if(UNI_ID_a < UNI_ID_b)
		return 1;
	return 0;
}

//the path indicate a file containing all sorted k+1mers
//set ref_list part for each unitig and return sorted REF_UNITIG list
void set_ref_lists(
		char *refPath,
		uint64_t* kmer_v,
		uint16_t *kmer_i,
		uint64_V *start_kmer_,
		DA_IDX *idx,
		uint64_t *counter)
{
	UNITIG_V *unitig_v = &(idx->unitig_v);
	REF_UNI_V *r_u_v = &(idx->r_u_v);
	REF_INFO_V * r_i_v = &(idx->r_i_v);
	uint64_t *start_kmer = start_kmer_->a;
	uint64_t unitig_num = start_kmer_->n;
	//for bin_ref
	idx->ref_bin.m = idx->ref_bin.n = ((idx->ref_bin.n + 3) >> 2);
	idx->ref_bin.a = (uint8_t*) xmalloc(sizeof(uint8_t)*idx->ref_bin.n);
	uint8_t* ref_bin = idx->ref_bin.a;
	// initate preparations, start to initiate in-out degrees this
	gzFile fp_z = xzopen(refPath,"r");
	kseq_t *trunk = kseq_init(fp_z);

	uint16_t move = (B_KMER - L_PRE_IDX) << 1;//(31 - 13)*2 = 36
	int ref_ID = -1;
	uint64_t seq_offset = 0;
	REF_UNI_V ref_uni_vec = {0,0,0};
	uint8_t bit_pos = 0;
	uint64_t offset_pos = 0;
	while (kseq_read(trunk)>=0)
	{//read ref seq
		size_t seq_l = trunk->seq.l;
		char * seq_s = trunk->seq.s;
		///(1):store ref info
		REF_INFO * c_ref;
		kv_pushp(REF_INFO, *r_i_v, &c_ref);
		ref_ID++;
		c_ref->seq_l = seq_l;
		strcpy(c_ref->ref_name, trunk->name.s);
		c_ref->seq_offset = seq_offset;
		seq_offset += seq_l;
		///(2):store bin ref
		for(uint64_t i = 0; i < seq_l; i++)
		{
			switch(bit_pos)
			{
				case 0: bit_pos = 1; ref_bin[offset_pos] = (bin_Bit[(int8_t)seq_s[i]]) << 6; break;
				case 1: bit_pos = 2; ref_bin[offset_pos] |= (bin_Bit[(int8_t)seq_s[i]]) << 4; break;
				case 2: bit_pos = 3; ref_bin[offset_pos] |= (bin_Bit[(int8_t)seq_s[i]]) << 2; break;
				case 3: bit_pos = 0; ref_bin[offset_pos] |= (bin_Bit[(int8_t)seq_s[i]]); offset_pos ++; break;
			}
		}
		///(3):store ref_uni
		for(uint64_t global_index = 0; global_index < seq_l; ++global_index)
		{
			if (Bit[(int8_t)seq_s[global_index]] >= 4)
				continue;
			uint64_t start = global_index;
			while(Bit[(int8_t)seq_s[++global_index]] != 4 && global_index < seq_l);
			if (start + B_KMER > global_index)
				continue;
			while(1){
				///for the forward string
				uint64_t key_start = char2Kmer(seq_s + start, B_KMER, Bit);
				uint64_t loc_start = binSearch(kmer_v, counter[key_start>>move],counter[(key_start>>move)+1]-1, key_start);
				xassert(isStart(kmer_i[loc_start]), "Must be a start.");
				//get the length of unitig:
				//get UNITIG ID
				uint32_t unitig_ID = binSearch(start_kmer, 0, unitig_num - 1, key_start);
				uint32_t uni_len = unitig_v->a[unitig_ID].length;
				uint64_t key_end = char2Kmer(seq_s + start + uni_len - B_KMER, B_KMER, Bit);
				uint64_t loc_end = binSearch(kmer_v, counter[key_end>>move],counter[(key_end>>move)+1]-1, key_end);
				xassert(isEnd(kmer_i[loc_end]), "Must be an end.");
				REF_UNI * r_u;
				if(uni_len >= MIN_UNI_L)
				{
					kv_pushp(REF_UNI, ref_uni_vec, &r_u);
					store_REF_UNITIG(r_u, ref_ID, unitig_ID, start, FORWARD);
				}
#ifdef CONSIDER_BOTH_ORIENTATION
				///for the reverse string
				key_start = char2KmerRC(seq_s + start + uni_len - 1, B_KMER, Bit);
				loc_start = binSearch(kmer_v, counter[key_start>>move],counter[(key_start>>move)+1]-1, key_start);
				xassert(isStart(kmer_i[loc_start]), "Must be a start.");
				//(Debug code)assert end
				//get the length of unitig:
				//get UNITIG ID
				unitig_ID = binSearch(start_kmer, 0, unitig_num - 1, key_start);
				uint32_t uni_len_R = unitig_v->a[unitig_ID].length;
				xassert(uni_len_R == uni_len, "");
				key_end = char2KmerRC(seq_s + start, B_KMER, Bit);
				loc_end = binSearch(kmer_v, counter[key_end>>move],counter[(key_end>>move)+1]-1, key_end);
				xassert(isEnd(kmer_i[loc_end]), "Must be an end.");
				//store
				if(uni_len >= MIN_UNI_L)
				{
					kv_pushp(REF_UNI, ref_uni_vec, &r_u);
					store_REF_UNITIG(r_u, ref_ID, unitig_ID, start, REVERSE);
				}
#endif
				///reset start
				start += uni_len - B_KMER + 1;
				if(start + B_KMER > global_index)
				{
					xassert(start + B_KMER == global_index + 1, "");
					break;
				}
			}
		}
	}
	xassert(idx->ref_bin.n == offset_pos || idx->ref_bin.n == offset_pos + 1,"Wrong in write bin_ref");
	kseq_destroy(trunk);
	gzclose(fp_z);
	free(start_kmer_->a);
	//sort ref_uni_vec
	if(ref_uni_vec.n > 1){
		int thread_number = 16;
		if(ref_uni_vec.n < 100){
		    thread_number = 1;
		}
		FUNC_GET_TIME_P(
		ksort_stable_mt(ref_uni_vec.a, ref_uni_vec.n, sizeof(REF_UNI), REF_UNITIG_cmp_by_UNITIG_ID, thread_number),
		"sort ref_uni",true);
	}

	//set ref_list part for each unitig
	UNITIG *unitig_p = unitig_v->a;
	REF_UNI * r_u_p = ref_uni_vec.a;
	uint32_t old_unitig_ID = MAX_uint32_t, c_uni_ID;
	for(uint32_t r_u_index = 0; r_u_index < ref_uni_vec.n; r_u_index++)
	{
		//skip short unitig
		c_uni_ID = r_u_p[r_u_index].UNITIG_ID;
		if(old_unitig_ID != c_uni_ID)
		{
			if(unitig_p[c_uni_ID].ref_list != 0)
			{
				xassert(unitig_p[c_uni_ID].ref_list == r_u_index, "");
			}
			else
				unitig_p[c_uni_ID].ref_list = r_u_index;
			unitig_p[c_uni_ID + 1].ref_list = r_u_index + 1;
			old_unitig_ID = c_uni_ID;
		}
		else
			unitig_p[c_uni_ID + 1].ref_list++;
	}

	UNITIG *unitig;
	kv_pushp(UNITIG, *unitig_v, &unitig);
	unitig->length = 0;
	unitig->ref_list = ref_uni_vec.n;
	r_u_v->m = r_u_v->n = ref_uni_vec.n;
	r_u_v->a = ref_uni_vec.a;
}

typedef struct
{
	//out_put
	UNITIG_V unitig_v;
	uint64_V start_kmer;
	uint64_V sp_kmer_ori;
	//input
	uint64_t *kmerValue;
	uint16_t *kmerInfo;
	uint64_t s_kmer;//start and end of kmer index
	uint64_t e_kmer;
	uint64_t *counter;
}UNI_V_DATA;

//out_put sp_kmer, sp_KMER_0, start_kmer_v, unitig_v
static void *get_uni_v_worker(void *data)
{
	UNI_V_DATA * d = (UNI_V_DATA*)data;
	//input
	uint64_t* kmer_v = d->kmerValue;
	uint16_t* kmer_i = d->kmerInfo;
	uint64_t *counter = d->counter;
	//output
	UNITIG_V *unitig_v = &(d->unitig_v);
	//values
	uint16_t lastChar = (d->s_kmer == 0)?(0x5 << 10):(0x4 << 10);//0001 0100 0000 0000. the last char, mean '$'
	uint16_t move = (B_KMER - L_PRE_IDX) << 1;//(31 - 13)*2 = 36
	uint64_t mask = ~((uint64_t)0x3 <<((B_KMER-1)<<1));//110011....11(60 '1')
	uint8_t  bitMove = ((B_KMER-1)<<1);//60
	//loop
	for (uint64_t i = d->s_kmer; i < d->e_kmer; ++i)
	{	//travel all the kmers for begin to end
		if (isStart(kmer_i[i]))
		{///start of a UNITIG
			kmer_i[i] |= lastChar;//store last char
			uint64_t loc = i;
			uint32_t L_unitig = B_KMER;//init the length of the UNITIG
			while(!isEnd(kmer_i[loc]))
			{//middle of a UNITIG
				lastChar = (kmer_v[loc] >> bitMove) << 10;//get the ACGT value
				uint8_t next_char = 0;//get the out-edge(next char)
				while (!((kmer_i[loc]>>next_char)&1)) ++next_char;
				uint64_t key = ((kmer_v[loc]&mask)<<2)|next_char;//get new key
				loc = binSearch(kmer_v, counter[key>>move], counter[(key>>move)+1]-1, key);//get next loc
				kmer_i[loc] |= lastChar;//store new last char
				L_unitig ++;//add the L_unitig
			}//end while loop end of a UNITIG
			kv_push(uint64_t, d->sp_kmer_ori, kmer_v[loc]);//#1 < #2 .... < $; they are all less than all AGCT
			//genSpKmers(kmer_v[loc], sp_kmer, sp_kmer_0, _ind++);
			lastChar = 0x4<<10;//reset last char
			UNITIG *unitig;
			kv_pushp(UNITIG, *unitig_v, &unitig);
			unitig->length = L_unitig;
			kv_push(uint64_t, d->start_kmer, kmer_v[i]);
		}
	}
	return NULL;
}

#define N_T_UNI_V 16
int get_uni_v(char* refPath,
		uint64_t* kmer_v,
		uint64_t n_kmer,
		uint16_t* kmer_i,
		//SP_KMER* sp_kmer,
		//uint8_t* sp_kmer_0,
		DA_IDX *idx,
		uint64_t *counter,
		uint64_t n_uni,
		uint64_t *sp_kmer_ori)//sp kmer
{
	///PART ONE: get length of the unitig
	pthread_t tid[N_T_UNI_V] = {0};
	UNI_V_DATA data[N_T_UNI_V] = {0};
	uint64_t step_len = n_kmer/N_T_UNI_V;
	for (int i = 0; i < N_T_UNI_V; ++i)
	{
		//input
		data[i].kmerValue = kmer_v;
		data[i].kmerInfo = kmer_i;
		data[i].counter = counter;
		data[i].s_kmer = step_len*i;
		data[i].e_kmer = step_len*(i + 1);
	}
	data[N_T_UNI_V - 1].e_kmer = n_kmer;
	for (int i = 0; i < N_T_UNI_V; ++i)
		pthread_create(&tid[i], 0, get_uni_v_worker, &data[i]);
	for (int i = 0; i < N_T_UNI_V; ++i)
		pthread_join(tid[i], 0);
	//combine output, UNITIG_V
	{
		size_t n_total = 0;
		for (int i = 0; i < N_T_UNI_V; ++i)
			n_total += data[i].unitig_v.n;
		UNITIG_V *u_v = &(idx->unitig_v);
		u_v->n = u_v->m = n_total;
		u_v->a = xcalloc_t(UNITIG, n_total + 1000);
		UNITIG * u_v_p = u_v->a;
		size_t idx_g = 0;//global index
		for (int i = 0; i < N_T_UNI_V; ++i)
		{
			UNITIG *data_i = data[i].unitig_v.a;
			for(size_t j = 0; j < data[i].unitig_v.n; j++)
				u_v_p[idx_g++].length = data_i[j].length;
			free(data[i].unitig_v.a);
		}
	}
	//combine start_kmer
	uint64_V start_kmer;
	{
		size_t n_total = 0;
		for (int i = 0; i < N_T_UNI_V; ++i)
			n_total += data[i].start_kmer.n;
		start_kmer.m = start_kmer.n = n_total;
		start_kmer.a = xcalloc_t(uint64_t, n_total);
		uint64_t * s_k_p = start_kmer.a;
		size_t idx_g = 0;//global index
		for (int i = 0; i < N_T_UNI_V; ++i)
		{
			uint64_t *data_i = data[i].start_kmer.a;
			for(size_t j = 0; j < data[i].start_kmer.n; j++)
				s_k_p[idx_g++] = data_i[j];
			free(data[i].start_kmer.a);
		}
	}

	//combine sp kmers
	{
		size_t n_total = 0;
		for (int i = 0; i < N_T_UNI_V; ++i)
			n_total += data[i].sp_kmer_ori.n;
		xassert(n_total == n_uni, "Wrong in unitig number!\n");
		size_t idx_g = 0;//global index
		for (int i = 0; i < N_T_UNI_V; ++i)
		{
			uint64_t *c_sp_kmer_ori = data[i].sp_kmer_ori.a;
			for(size_t j = 0; j < data[i].sp_kmer_ori.n; j++)
				sp_kmer_ori[idx_g++] = c_sp_kmer_ori[j];
				//genSpKmers(c_sp_kmer_ori[j], sp_kmer, sp_kmer_0, idx_g++);
			free(data[i].sp_kmer_ori.a);
		}
	}
	//PART TWO: set ref lists for each unitig
	fprintf(stderr,"set_ref_lists\n");
	set_ref_lists(refPath, kmer_v, kmer_i, &start_kmer, idx, counter);
	return 0;
}

static inline int spkmer_cmp_l(const void* a_, const void * b_){//less than
	SP_KMER * a = (SP_KMER*) a_, * b = (SP_KMER*) b_;
	uint64_t *a_kmer_value = (uint64_t*)(a->value);
	uint64_t *b_kmer_value = (uint64_t*)(b->value);
	if (a->sp_pos < b->sp_pos ){
		uint8_t move = (b->sp_pos - a->sp_pos) << 1;
		if(*a_kmer_value <= ((*b_kmer_value)>>move))
			return 1;//compare original kmer, must be "<=", not "<"
		return -1;
	}
	if(a->sp_pos > b->sp_pos )
	{
		uint8_t move = (a->sp_pos - b->sp_pos) << 1;
		if(((*a_kmer_value) >> move) < (*b_kmer_value))
			return 1;
		return -1;
	}
	else{//equal in pos
		if((*a_kmer_value) < (*b_kmer_value))
			return 1;
		if((*a_kmer_value) > (*b_kmer_value))
			return -1;
		return 0;
	}
}//special KMER with '#' or '$'

//input : file pathes
//output: BWT/unitig list/hash table
void build_UNITIG(DA_IDX *idx, char* refPath, const char *kmerPath, kstring_t *bwt_s)
{
	///PART 1: read sorted kmers from file
	FILE *fp = xopen(kmerPath, "rb");	//read kmer file, kmer are store in 64 bit
	xread(&(idx->n_kmer), sizeof(uint64_t), 1 , fp);
	fprintf(stderr,"start reading KMER file......\n");
	uint64_t *kmer_v = xmalloc_t(uint64_t, idx->n_kmer);//120G; 15G*8
	uint16_t *kmer_i = xmalloc_t(uint16_t, idx->n_kmer);//30G+120G=150G
	xread(kmer_v, sizeof(uint64_t), idx->n_kmer, fp);
	///PART 2: start build DBG, building by reading from the start of sequence to the end
	fprintf(stderr,"Start building deBG, kmer number is[%lu]\n",idx->n_kmer);
	uint64_V heads = {0}, tails = {0};
	uint64_t *counter = getCounter(kmer_v, idx->n_kmer);//0.5G, used for faster kmer binary search
	FUNC_GET_TIME_P(
			build_deb_mt(refPath, kmer_v, kmer_i, &heads, &tails, counter, &(idx->ref_bin.n)),
	"t_build_deb", true);//get search dir
	//build_deb(refPath, kmer_v, kmer_i, idx->n_kmer, &heads, &tails, counter, &(idx->ref_bin.n));
	uint64_t n_uni = 0;
	FUNC_GET_TIME_P(
	n_uni = setLabel(kmer_v, kmer_i, idx->n_kmer, &heads, &tails, counter),//number of UNITIG
	"t_setLabel", true);//get search dir
	free(heads.a);
	free(tails.a);
	fprintf(stderr, "Number of UNITIG is [%lu]\n",n_uni);
	//PART 3: handle special kmers and build UNITIG
	fprintf(stderr,"Start building UNITIG... \n");
	uint64_t *sp_kmer_ori = xmalloc_t(uint64_t, n_uni);//1.6G
	//FUNC_GET_TIME_P(
	get_uni_v(refPath, kmer_v, idx->n_kmer, kmer_i, idx, counter, n_uni, sp_kmer_ori);//,
	//"t_get_uni_v", true);
	//free kmer_v to free 120G memory for sp_kmer sort
    free(kmer_v);//150 -120 = 30G mem

	//handle SP kmer
	uint64_t n_sp_kmer = n_uni*(B_KMER - 1);
	SP_KMER* sp_kmer   = xmalloc_t(SP_KMER, n_sp_kmer);//total (unitigNum * B_KMER) new-kmer which include '#', 200M*30*16 = 96G
	uint8_t* sp_kmer_0 = xmalloc_t(uint8_t, n_uni);//200M
	for (int i = 0; i < n_uni; ++i)
		genSpKmers(sp_kmer_ori[i], sp_kmer, sp_kmer_0, i);
	free(sp_kmer_ori);
	for (uint64_t i = 0; i < n_uni; ++i)
		kputc(Chars[sp_kmer_0[i] & 0x7], bwt_s);
	free(sp_kmer_0);
	FUNC_GET_TIME_P(
	ksort_stable_mt(sp_kmer, n_sp_kmer, sizeof(SP_KMER), spkmer_cmp_l, 15),//another 96G was allocated
	"t_ksort_stable", true);
	fprintf(stderr,"Start building BWT string.... \n");

	//merge sort and generate hash index table
	uint64_t hash_index_size = 0x1l <<((L_PRE_IDX<<1) + 1);//2^27 = 128M
	idx->bt.hash_index = xmalloc_t(uint64_t, hash_index_size);
	uint64_t* hash_index = idx->bt.hash_index;
    for(int idx = 0; idx < hash_index_size; idx++)
    	hash_index[idx] = MAX_uint64_t;
    fseek(fp, 8, SEEK_SET);
	merge_kmer(fp, idx->n_kmer, kmer_i, sp_kmer, n_sp_kmer, bwt_s, idx->bt.hash_index);
	//free memory
	free(sp_kmer);
	free(kmer_i);

    //part 4: compress bwt hash index
    uint64_t l_cps_hi = (0x1l << (L_PRE_IDX<<1)) + 1;//compressed hash index size
    uint64_t* cps_hi = xcalloc_t(uint64_t, l_cps_hi);//compressed hash index
    uint64_t pre_value = 0;
    for(int idx = 0; idx < hash_index_size; idx += 2)
    {
    	//hash_index
    	if(hash_index[idx] != MAX_uint64_t)
    	{
    		cps_hi[idx >> 1] = hash_index[idx];
    		pre_value = hash_index[idx + 1];
    	}
    	else
    		cps_hi[idx >> 1] = pre_value;
    }
    cps_hi[l_cps_hi - 1] = pre_value;
    free(idx->bt.hash_index);
    idx->bt.hash_index = cps_hi;
    fclose(fp);
}

#define SINGLE_BASE_MAX_RATIO_THEADHOLD 0.8
void set_ekmer_par(E_KMER *ek){
	uint64_t e_kmer_hash_mask = MASK_37;
	uint8_t  len_e_kmer = 20;
	switch(ek->e_kmer_size){
		case one_eighth_GigaByte:	e_kmer_hash_mask = MASK_30; len_e_kmer = 16; break;
		case quarter_GigaByte: 	 	e_kmer_hash_mask = MASK_31; len_e_kmer = 17; break;
		case half_GigaByte: 		e_kmer_hash_mask = MASK_32; len_e_kmer = 17; break;
		case one_GigaByte: 			e_kmer_hash_mask = MASK_33; len_e_kmer = 18; break;
		case two_GigaByte: 			e_kmer_hash_mask = MASK_34; len_e_kmer = 18; break;
		case four_GigaByte: 		e_kmer_hash_mask = MASK_35; len_e_kmer = 19; break;
		case eight_GigaByte: 		e_kmer_hash_mask = MASK_36; len_e_kmer = 19; break;
		case sixteen_GigaByte: 		e_kmer_hash_mask = MASK_37; len_e_kmer = 20; break;
	}
	ek->e_kmer_hash_mask = e_kmer_hash_mask;
	ek->len_e_kmer = len_e_kmer;
	ek->single_base_max = SINGLE_BASE_MAX_RATIO_THEADHOLD * len_e_kmer;
}

//#define L_E_kmer 16
extern uint64_t kmerMask[33];
void get_EXIST_kmer(E_KMER *ek, UNITIG_V *unitig_v_, uint8_V * uni_str, uint64_t kmerNum){
	//set e_kmer_size
	if	   (kmerNum < (0x1l << 31) / 9)	ek->e_kmer_size = one_eighth_GigaByte;//0.2G
	else if(kmerNum < (0x1l << 32) / 9)	ek->e_kmer_size = quarter_GigaByte;//0.4G
	else if(kmerNum < (0x1l << 33) / 9)	ek->e_kmer_size = half_GigaByte;//0.8G
	else if(kmerNum < (0x1l << 34) / 9)	ek->e_kmer_size = one_GigaByte;//1.6G
	else if(kmerNum < (0x1l << 35) / 9)	ek->e_kmer_size = two_GigaByte;//3.2G
	else if(kmerNum < (0x1l << 36) / 9)	ek->e_kmer_size = four_GigaByte;//6.4G
	else if(kmerNum < (0x1l << 37) / 9)	ek->e_kmer_size = eight_GigaByte;//12.8G
	else if(kmerNum < (0x1l << 38) / 9)	ek->e_kmer_size = sixteen_GigaByte;//25.6G//the MAX
	else{ek->e_kmer_size = sixteen_GigaByte; fprintf(stderr, "Warning: too many kmers!\n");}//the MAX
	set_ekmer_par(ek);
	fprintf(stderr,"l_e_kmer: %d",ek->len_e_kmer);

	ek->e_kmer0 = xcalloc_t(uint8_t, ek->e_kmer_size);
	ek->e_kmer1 = xcalloc_t(uint8_t, ek->e_kmer_size);
	uint8_t* e_kmer0 = ek->e_kmer0,* e_kmer1 = ek->e_kmer1;

	UNITIG * unitig_v = unitig_v_->a;
	uint64_t L_unitig_v = unitig_v_->n - 1;
	uint64_t MASK = kmerMask[ek->len_e_kmer];
	uint64_t e_kmer_hash_mask = ek->e_kmer_hash_mask;
	//for each UNITIG, store exist kmers
	uint8_t* unitig_str = uni_str->a;
	uint8_t len_e_kmer = ek->len_e_kmer;
	uint64_t s = 0, e = s + unitig_v[0].length - len_e_kmer + 1;
	for(uint64_t u = 0; u < L_unitig_v; u++)
	{
		uint64_t kmer = bit2_nextKmer_init(unitig_str + s, len_e_kmer);
		for(uint64_t i = s; i < e; i++)
		{//for each kmer
			kmer = bit2_nextKmerMASK(unitig_str + i, kmer, len_e_kmer);
			uint64_t hash1 = hash64_1(kmer) & e_kmer_hash_mask;
			e_kmer0[hash1 >> 3] = (e_kmer0[hash1 >> 3] | (0x80 >> (hash1 & 0x7)));
			uint64_t hash2 = hash64_2(kmer) & e_kmer_hash_mask;
			e_kmer1[hash2 >> 3] = (e_kmer1[hash2 >> 3] | (0x80 >> (hash2 & 0x7)));
		}
		s += (unitig_v[u].length + 1);//reset s/e
		e  = s + unitig_v[u + 1].length - len_e_kmer + 1;
	}
}

//abandoned: store uni_string in binary
void bin_uni_str(uint8_V * uni_str)
{
	uint8_t* unitig_str = uni_str->a;
	int in_word_offset = 0;
	uint64_t global_offset = 0;
	for(uint64_t i = 0; i < uni_str->n; i++)
	{
		switch(in_word_offset){
			case 0: unitig_str[global_offset]  = (unitig_str[i] & 0x3) << 6; in_word_offset = 1; break;
			case 1: unitig_str[global_offset] |= (unitig_str[i] & 0x3) << 4; in_word_offset = 2; break;
			case 2: unitig_str[global_offset] |= (unitig_str[i] & 0x3) << 2; in_word_offset = 3; break;
			case 3: unitig_str[global_offset] |= (unitig_str[i] & 0x3)     ; in_word_offset = 0; global_offset++; break;
		}
	}
	uni_str->n = ((uni_str->n + 3) >> 2);
}

void write_idx(DA_IDX *idx, const char *dirPath){
	///FIRST PART: BWT
	write_bwt(&(idx->bt), dirPath);
	//SECOND PART: E_KMER
	FILE *exk0_fp = xget_file(dirPath, ".exk0", "wb");//exist kmer
	FILE *exk1_fp = xget_file(dirPath, ".exk1", "wb");//exist kmer
	FILE *exki_fp = xget_file(dirPath, ".exki", "wb");//exist kmer
	//EXIST_kmer length: 17 or 18
	E_KMER *ek = &(idx->ek);
	fwrite(&(ek->e_kmer_size), sizeof(uint64_t), 1, exki_fp);
	//0 EXIST_kmer table
	fwrite(ek->e_kmer0, sizeof(uint8_t), ek->e_kmer_size, exk0_fp);
	//1 EXIST_kmer table
	fwrite(ek->e_kmer1, sizeof(uint8_t), ek->e_kmer_size, exk1_fp);
	fclose(exk0_fp); fclose(exk1_fp); fclose(exki_fp);
	//THIRD PART: UNITIG vector
	FILE *unv_fp = xget_file(dirPath, ".unv", "wb");//UNITIG vector
	fwrite(&(idx->unitig_v.n), sizeof(uint64_t), 1, unv_fp);
	fwrite(idx->unitig_v.a, sizeof(UNITIG), idx->unitig_v.n, unv_fp);
	fclose(unv_fp);

/*
	//UNITIG string
	FILE *uni_fp = xget_file(dirPath, ".uni", "wb");
	fwrite(&(idx->uni_str.n), sizeof(uint64_t), 1, uni_fp);
	fwrite(idx->uni_str.a, sizeof(uint8_t), idx->uni_str.n, uni_fp);
	fclose(uni_fp);
*/

	///store bin ref
	FILE *ref_b_fp = xget_file(dirPath, ".ref_b", "wb");
	fwrite(&(idx->ref_bin.n), sizeof(uint64_t), 1, ref_b_fp);
	fwrite(idx->ref_bin.a, sizeof(uint8_t), idx->ref_bin.n, ref_b_fp);
	fclose(ref_b_fp);

	//store ref info lists
	FILE *ref_i_fp = xget_file(dirPath, ".ref_i", "wb");
	fwrite(&(idx->r_i_v.n), sizeof(uint64_t), 1, ref_i_fp);
	fwrite(idx->r_i_v.a, sizeof(REF_INFO), idx->r_i_v.n, ref_i_fp);
	fclose(ref_i_fp);

	///store ref-uni list
	//conver ref-uni to ref-pos
	idx->r_p_v.n = idx->r_u_v.n; idx->r_p_v.a = xmalloc(idx->r_p_v.n * sizeof(REF_POS));
	REF_POS *bgn_rp = idx->r_p_v.a, *end_rp = bgn_rp + idx->r_p_v.n; REF_UNI * bgn_ru = idx->r_u_v.a;
	for(;bgn_rp < end_rp; bgn_rp++, bgn_ru++){
		//store global offset
		bgn_rp->ref_ID = bgn_ru->ref_ID;
		bgn_rp->direction = bgn_ru->direction;
		bgn_rp->global_offset = idx->r_i_v.a[bgn_rp->ref_ID].seq_offset + bgn_ru->ref_offset;//todo:BUG when dir == reverse
	}
	FILE *ref_p_fp = xget_file(dirPath, ".ref_p", "wb");
	fwrite(&(idx->r_p_v.n), sizeof(uint64_t), 1, ref_p_fp);
	fwrite(idx->r_p_v.a, sizeof(REF_POS), idx->r_p_v.n, ref_p_fp);
	fclose(ref_p_fp);
}

void load_idx(DA_IDX *idx, const char *dirPath)
{
	///FIRST PART: BWT
	load_bwt(&(idx->bt), dirPath);
	//SECOND PART: E_KMER
	FILE *exk0_fp = xget_file(dirPath, ".exk0", "rb");//exist kmer
	FILE *exk1_fp = xget_file(dirPath, ".exk1", "rb");//exist kmer
	FILE *exki_fp = xget_file(dirPath, ".exki", "rb");//exist kmer
	E_KMER *ek = &(idx->ek);
	xread(&(ek->e_kmer_size), sizeof(uint64_t), 1, exki_fp);
	//set e_kmer
	set_ekmer_par(ek);
	ek->e_kmer0 = (uint8_t*)xmalloc(sizeof(uint8_t) * ek->e_kmer_size);
	ek->e_kmer1 = (uint8_t*)xmalloc(sizeof(uint8_t) * ek->e_kmer_size);
	xread(ek->e_kmer0, sizeof(uint8_t), ek->e_kmer_size, exk0_fp);
	xread(ek->e_kmer1, sizeof(uint8_t), ek->e_kmer_size, exk1_fp);
	fclose(exk0_fp);
	fclose(exk1_fp);
	fclose(exki_fp);
	//THIRD PART: UNITIG vector
	FILE *unv_fp = xget_file(dirPath, ".unv", "rb");
	xread(&(idx->unitig_v.n), sizeof(uint64_t), 1, unv_fp);
	idx->unitig_v.a = (UNITIG*)xmalloc(sizeof(UNITIG)*(idx->unitig_v.n + 1));
	xread(idx->unitig_v.a, sizeof(UNITIG), idx->unitig_v.n, unv_fp);
	idx->unitig_v.a[idx->unitig_v.n].ref_list = idx->unitig_v.a[idx->unitig_v.n - 1].ref_list + 1 + idx->unitig_v.a[idx->unitig_v.n - 1].length;
	idx->bt.DOLLOR_POS = idx->unitig_v.n - 1 - 1;
	fclose(unv_fp);

/*
	//FORTH PART: UNITIG string
	FILE *uni_fp = xget_file(dirPath, ".uni", "rb");
	xread(&(idx->uni_str.n), sizeof(uint64_t), 1, uni_fp);
	idx->uni_str.a = (uint8_t*)xmalloc(sizeof(uint8_t)*idx->uni_str.n);
	xread(idx->uni_str.a, sizeof(uint8_t), idx->uni_str.n, uni_fp);
	fclose(uni_fp);
*/

	///store bin ref
	FILE *ref_b_fp = xget_file(dirPath, ".ref_b", "rb");
	xread(&(idx->ref_bin.n), sizeof(uint64_t), 1, ref_b_fp);
	idx->ref_bin.a = xmalloc_t(uint8_t, idx->ref_bin.n);
	xread(idx->ref_bin.a, sizeof(uint8_t), idx->ref_bin.n, ref_b_fp);
	fclose(ref_b_fp);

	//store ref info lists
	FILE *ref_i_fp = xget_file(dirPath, ".ref_i", "rb");
	xread(&(idx->r_i_v.n), sizeof(uint64_t), 1, ref_i_fp);
	idx->r_i_v.a = xmalloc_t(REF_INFO, idx->r_i_v.n);
	xread(idx->r_i_v.a, sizeof(REF_INFO), idx->r_i_v.n, ref_i_fp);
	fclose(ref_i_fp);

	///store ref-uni list
	FILE *ref_p_fp = xget_file(dirPath, ".ref_p", "rb");
	xread(&(idx->r_p_v.n), sizeof(uint64_t), 1, ref_p_fp);
	idx->r_p_v.a = xmalloc_t(REF_POS, idx->r_p_v.n);
	xread(idx->r_p_v.a, sizeof(REF_POS), idx->r_p_v.n, ref_p_fp);
	fclose(ref_p_fp);
}

extern uint8_t bwt_Bit[];
int bwt_cal_SA(bwt *bt, UNITIG_V *unitig_v_, uint8_V *uni_str)
{
	//malloc for SA
	bt->sa_size = (bt->len_bwt_str + NUM_BP_OCCUPY_ONE_SA - 1)/NUM_BP_OCCUPY_ONE_SA;//32 bp has an SA point
	bt->sa_taxon = xmalloc_t(SA_taxon, bt->sa_size);

	uni_str->n = uni_str->m = bt->len_bwt_str;
	uni_str->a = xmalloc_t(uint8_t, uni_str->n);//used for store unitig-string
	int64_t unitigCount = uni_str->n - 1;//the end index for the UNITIG-string, count in base pair
	uint32_t unitig_ID = unitig_v_->n - 1 - 1;//init ID and offset

	UNITIG * unitig_v = unitig_v_->a;
	SA_taxon* sa_taxon = bt->sa_taxon;
	uint8_t* unitig_str = uni_str->a;
	char* bwt_str = bt->bwt_str;

	uint32_t offset = unitig_v[unitig_ID].length - 1;
	//get the occ of '$'
	bt->DOLLOR_POS = unitig_v_->n - 1 - 1;
	uint64_t occ = bt->DOLLOR_POS;
	//store '$'
	//unitig_str[unitigCount>>1] |= 0x101 << ((unitigCount & 0x1) << 2); --unitigCount;
	unitig_str[unitigCount--] = 0x5;
	//get the char before '$'
	uint8_t c = bwt_Bit[(uint8_t)bwt_str[occ]];
	//store c
	unitig_str[unitigCount--] = c;
	//unitig_str[unitigCount>>1] |= c << ((unitigCount & 0x1) << 2); --unitigCount;
	//store SA
	if((occ & SA_MASK) == 0)
	{//SA point
		sa_taxon[occ >> SA_OFF].unitig_ID = unitig_ID;
		sa_taxon[occ >> SA_OFF].offset = offset;
	}
	offset --;
	while(1)
	{//LFC loop
		//get new occ
		occ = build_LFC(bt, occ,c);
		c = bwt_Bit[(uint8_t)bwt_str[occ]];
		if(c == 4)
		{
			if(unitig_ID == 0)
        		    break;
			unitig_ID--;
			xassert(offset == MAX_uint32_t, "Not end of offset!\n");
			offset = unitig_v[unitig_ID].length;
		}//'#'
		if(c == 5)
		{
			xassert(offset == MAX_uint32_t, "Not end of offset!\n");
			break;
		}//mean '$'
		//store c
		//unitig_str[unitigCount>>1] |= c << ((unitigCount & 0x1) << 2); --unit11igCount;
		unitig_str[unitigCount--] = c;
		//store SA
		if((occ & SA_MASK) == 0)
		{//SA point
			sa_taxon[occ >> SA_OFF].unitig_ID = unitig_ID;
			sa_taxon[occ >> SA_OFF].offset = offset;
		}
		offset --;
	}
	//for the last '$'
	if((occ & SA_MASK) == 0)
	{
		unitig_ID = unitig_v_->n - 1 - 1;
		offset = unitig_v[unitig_ID].length;
		sa_taxon[occ >> SA_OFF].unitig_ID = unitig_ID;
		sa_taxon[occ >> SA_OFF].offset = offset;
	}
	xassert(unitigCount == MAX_uint64_t, "DEBUG: Wrong unitigCount");
	return 0;
}

static int ind_usage(){
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:   %s\n", PACKAGE_NAME);
	fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact:   %s\n\n", CONTACT);
	fprintf(stderr, "  Usage:     %s  index  <Options> [SortedKmer] [Reference] [IndexDir]\n", PACKAGE_NAME);
	fprintf(stderr, "  Basic:     \n");
	fprintf(stderr, "    [SortedKmer]  FILE   sorted kmers file \"kmer.srt\" genereted by \"kmersort\"\n");
	fprintf(stderr, "    [Reference]   FILE   one fasta REF file, multiple files need to be combined\n");
	fprintf(stderr, "    [IndexDir]    FOLDER the directory to store %s index\n", PACKAGE_NAME);
	fprintf(stderr, "  Options:\n");
	fprintf(stderr, "    -h            help\n");
	fprintf(stderr, "\n");
	return 0;
}

int build_index_main(int argc, char *argv[]){
	int c;
	while((c = getopt(argc, argv, "k:h")) >= 0)
	{
		if(c == 'h') return ind_usage();
	}

	if(optind + 3 > argc)
	{
		fprintf(stderr,"Warning: [SortedKmer] / [Reference] or [IndexDir] missed.\n");
		return ind_usage();
	}
	char *kmerPath = argv[optind++];
	char *refPath  = argv[optind++];
	char *dirPath  = argv[optind++];

	DA_IDX idx = {0};
	kstring_t bwt_s = {0};//used to store BWT string
	build_UNITIG(&idx, refPath, kmerPath, &bwt_s);
	FUNC_GET_TIME_P(build_BWT(&(idx.bt), bwt_s.s, bwt_s.l, dirPath),"build_BWT", true);
	//store SA part
	FUNC_GET_TIME_P(bwt_cal_SA(&(idx.bt), &(idx.unitig_v), &(idx.uni_str)),
			"bwt_cal_SA", true);
	FUNC_GET_TIME_P(get_EXIST_kmer(&(idx.ek), &(idx.unitig_v), &(idx.uni_str), idx.n_kmer),//e_kmer
			"get_EXIST_kmer", true);
	write_idx(&idx, dirPath);
	return 0;
}
