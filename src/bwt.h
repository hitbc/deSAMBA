#ifndef _BWT_H
#define _BWT_H
#include <stdint.h>

//defines for SA
#define NUM_BP_OCCUPY_ONE_SA 8
#define SA_MASK 0x7
#define SA_OFF 3

typedef struct{
	uint32_t unitig_ID;
	uint32_t offset;
}SA_taxon;

typedef struct{
	//USED ONLY IN BUILDING
	char*			bwt_str;//the BWT string when build BWT, 8bit for each bp
	uint64_t 		len_bwt_str;//length of bwt

	uint64_t* 		occCheck;//check point list when building
	uint64_t 		occCount;//number of check point when building

	uint8_t* 		bin_bwt;//BWT string when build, 4 bit for each base pair
	uint64_t 		len_bwt_bin;//length of bin BWT when build

	///USED ONLY IN CLASSIFYING
	//FM-INDEX
	uint8_t* 		bwt_occ;//BWT string and FM index when search, 4 bit for each base pair
	uint64_t		byteLen;//length of BWT string and FM index when search

	//USED ONLY IN BOTH
	uint64_t		rank[6];//rank, {'A', 'C', 'G', 'T', '#'}

	uint8_t* 		AGCTCounter[5];//300k counter

	SA_taxon* 		sa_taxon;
	uint64_t 		sa_size;

	uint64_t*		hash_index;

	uint64_t 		DOLLOR_POS;//the position of "$" in the bwt string, END of original string, equal (number of UNITIG - 1)
}bwt;

void build_BWT(bwt *bt, char *bwt_string, uint64_t bwt_len, char *dirPath);
uint64_t occ(bwt* bt, uint64_t r, uint8_t *c);
uint64_t build_LFC(bwt *bt, uint64_t r, uint8_t c);
int write_bwt(bwt *bt, const char *dirPath);
int load_bwt(bwt * bt, const char *dirPath);
#endif
