#include "bwt.h"

#include "lib/utils.h"
#include <string.h>
#include <stdlib.h>

#include "lib/desc.h"

#define L_PRE_IDX 13

uint8_t bwt_Bit[] =
{	//A->0 ; C->1; G->2; T->3; 4->4; 5->5;
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4/*#*/,5/*$*/,4, 4, 4,  		4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4,
	4, 0,/*A*/ 4, 	1/*C*/, 4, 4, 4, 2/*G*/,4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	3/*T*/, 4, 4, 4,  		4, 4, 4, 4,  4, 4, 4, 4,
	4, 0,/*A*/ 4, 	1/*C*/, 4, 4, 4, 2/*G*/,4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	3/*T*/, 4, 4, 4,  		4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  	4, 4, 4, 4,  			4, 4, 4, 4,  4, 4, 4, 4
};
///-------------search function----------------------------------------------------

// to find out how many character 'c' are occurrence before the 'c' in pos r
//for each 256 bp bwt data, one bwt_occ block are stored
//the first 40 byte store the 64-bit *5 "base" part, the last 128 byte
//store the 4-bit * 256 "data" part
//for the base part, it store the total numbers of ACGT# before the 256 bp part
//for the data part, each base pair(one of ACGT#) occupy 4-bit,
//to count the total number of c before the c in position r, you should
//add the base of c of that block and the cs in that 256 block
//before the c of r, in this implement, 4 bps are count one time to speed up
//(by check the AGCTCounter table)
uint16_t occ_mask[4] = {0xFFFF, 0xFFF0, 0xFF00, 0xF000};//constricted by operating system only for linux
uint64_t occ(bwt *bt, uint64_t r, uint8_t *c)
{
	uint64_t pos = r >> 8;
	uint64_t p_occ = pos * 168;
	uint64_t count = 0;
	uint16_t *s_dc = (uint16_t *)(bt->bwt_occ + p_occ + 40);//start of double char
	uint16_t *e_dc = s_dc + ((r & 0xff) >> 2);
	if(*c == MAX_UINT8_t)
	{
		*c = ((*e_dc) >> ((r & 0x3) << 2)) & 0xf;//current chat
		xassert(*c <= 5, "In [bwt::occ]:c > 5\n");
		if(*c == 5)
			return bt->DOLLOR_POS;
	}
	else if(*c > 5)
		fprintf(stderr, "In [bwt::occ]:c > 5\n");
	uint64_t base = *(uint64_t *)(bt->bwt_occ + p_occ + ((*c)<<3));
	uint8_t* p_AGCTCounter = bt->AGCTCounter[*c];
	for (; s_dc < e_dc; s_dc++)//[0~64 times loop, mean 32]
		count += p_AGCTCounter[*s_dc];//here won't surpass the range
	count += p_AGCTCounter[(*s_dc) | occ_mask[r & 0x3]];
	return base + count;
}

//load bwt and taxonID
int load_bwt(bwt * bt, const char *dirPath)
{
	FILE *bwt_fp = xget_file(dirPath, ".bwt", "rb");
	FILE *sa_fp = xget_file(dirPath, ".sa", "rb");
	FILE *acg_fp = xget_file(dirPath, ".acg", "rb");

	///part one: FM index length
	xread(&(bt->byteLen), 8, 1, bwt_fp);//8
	//part two FM index
	bt->bwt_occ  = xmalloc_t(uint8_t, bt->byteLen);//21÷32×kmerNum
	xread(bt->bwt_occ, sizeof(uint8_t), bt->byteLen, bwt_fp);
	//part three: rank
	xread(bt->rank, sizeof(uint64_t), 5, bwt_fp);//40
	bt->rank[5] = bt->rank[0] - 1;//set rank for '&', it is before 'A'
	///part four: hash index: 1G
	uint64_t hash_index_size = (0x1l << (L_PRE_IDX<<1)) + 1;
	bt->hash_index = xmalloc_t(uint64_t, hash_index_size);//1G
	xread(bt->hash_index,sizeof(uint64_t),hash_index_size,bwt_fp);
	///part five: ACGT counter size
	uint64_t AGCTCounterSize;
	xread(&AGCTCounterSize, sizeof(uint64_t), 1, acg_fp);
	///part six: AGCTCounter
	for(int i = 0; i < 5; i++){
		bt->AGCTCounter[i] = xmalloc_t(uint8_t, AGCTCounterSize);
		xread(bt->AGCTCounter[i], sizeof(uint8_t), AGCTCounterSize, acg_fp);
	}
	///part seven: SA part size
	xread(&(bt->sa_size), sizeof(uint64_t), 1, sa_fp);//1÷32×8×kmerNum
	bt->sa_taxon = (SA_taxon*)xmalloc(sizeof(SA_taxon)*bt->sa_size);
	//part eight: SA part
	xread(bt->sa_taxon, sizeof(SA_taxon), bt->sa_size, sa_fp);
	//part ten: unitig file
	fclose(bwt_fp);
	fclose(sa_fp);
	fclose(acg_fp);
	return 0;
}

///-------------------------build BWT---------------------------------------------//

#define BP_PER_BLOCK 256//base pair per block
//calculate and store check point("base" part) and rank
void bwt_cal_check_point(bwt* bt)
{
	//256 is interval
	uint64_t *rank = bt->rank;
	uint64_t len_bwt_str = bt->len_bwt_str;
	bt->occCount = ((len_bwt_str + BP_PER_BLOCK - 1)/BP_PER_BLOCK)*5;//total number of check point
	bt->occCheck = (uint64_t*)xmalloc(bt->occCount * sizeof(uint64_t));//the check point list
	uint64_t* occCheck = bt->occCheck;
	char* bwt_str = bt->bwt_str;
	uint64_t ind = 0;//index of check point
	uint64_t tempOcc[] = {0,0,0,0,0,0};//store current check point{A C G T # $}
	uint8_t kago = 0;//store check point for each 256-bp block
	for(uint64_t i=0; i < len_bwt_str; ++i)
	{
		if (kago  == 0)
		{
			for(int j = 0; j < 5; j++)
				occCheck[ind + j] = tempOcc[j];
			ind += 5;
		}
		++tempOcc[bwt_Bit[(int8_t)bwt_str[i]]];
		++kago;
	}
	rank[0] = tempOcc[4] + tempOcc[5];//A after '#' and '$'
	rank[1] = rank[0]    + tempOcc[0];//C
	rank[2] = rank[1]    + tempOcc[1];//G
	rank[3] = rank[2]    + tempOcc[2];//T
	rank[4] = 0;					  //#
}

//BWT str include 6 kind character : A C G T 4 (means '#') 5 (means '$')
//store those chars in 4-bit: A->0 ; C->1; G->2; T->3; 4->4; 5->5;
void bwt_str2bwt_occ(bwt *bt)
{
	bt->len_bwt_bin = ((bt->len_bwt_str + 1) >> 1);
	bt->bin_bwt = (uint8_t*)xmalloc(bt->len_bwt_bin * sizeof(uint8_t));
	uint64_t leftOneChar = bt->len_bwt_str & 0x1;
	uint64_t i;
	uint8_t* bin_bwt = bt->bin_bwt;
	char* bwt_str = bt->bwt_str;
	uint64_t len_bwt_str = bt->len_bwt_str;
	for(i=0; i< len_bwt_str - leftOneChar; i += 2)
		bin_bwt[i>>1] = (bwt_Bit[(int8_t)bwt_str[i+1]]<<4)|(bwt_Bit[(int8_t)bwt_str[i]]);//used for little endian mode
	if (leftOneChar)
		bin_bwt[i>>1] = 0xF0|bwt_Bit[(int8_t)bwt_str[len_bwt_str-1]];//used for little endian mode
}

uint8_t countZero(uint64_t v)
{
	uint8_t counter = 0;
	for (uint8_t i=0; i<4;++i)
	{
		if ((v & 0xf) == 0)
			++counter;
		v >>= 4;
	}
	return counter;
}

void bwt_cal_AGCTCounter(bwt* bt)
{
	uint64_t AGCTCounterSize = 1 << 16;//64k*5 ~ 320k
	uint64_t mask[5] = {0x0000,0x1111,0x2222,0x3333,0x4444};
	uint8_t* p_AGCTCounter;
	for(uint8_t j=0; j < 5; ++j)
	{
		bt->AGCTCounter[j] = xmalloc_t(uint8_t, AGCTCounterSize);
		p_AGCTCounter = bt->AGCTCounter[j];
		uint64_t c_mask = mask[j];
		for (uint64_t i=0; i <= 0xFFFF; ++i)//64k
			p_AGCTCounter[i] = countZero(i^c_mask);
	}
}

//the LFC function when build BWT
uint64_t build_LFC(bwt *bt, uint64_t r, uint8_t c)
{
	static uint16_t mask[] = {0xFFFF, 0xFFF0, 0xFF00, 0xF000};//constricted by operating system only for linux
	uint64_t pos = r >> 8;//the block ID
	uint64_t base = bt->occCheck[pos*5 + c];
	uint64_t extra = r & 0x3;//the last 4 bp
	uint64_t start = pos << 7;//start of a block
	uint64_t end = (r - extra) >> 1;//the begin of the last 4 bp
	uint64_t count = 0;
	uint8_t *p_AGCTCounter = bt->AGCTCounter[c];
	//change every 4 bp block to count
	uint8_t* bin_bwt = bt->bin_bwt;
	for (uint64_t i = start; i < end; i += 2)
	{//[0~64 times loop, mean 32]
		count += p_AGCTCounter[(*(uint16_t *)(bin_bwt+i))];//here won't surpass the range
	}
	count += p_AGCTCounter[(*(uint16_t *)(bin_bwt + end) | mask[extra])];//for the last 4 bp, used for little endian mode
	return bt->rank[c] + base + count;//get the char in the first line
}

#define BWT_BYTE_PER_BLOCK (BP_PER_BLOCK/2)
int write_bwt(bwt *bt, const char *dirPath)
{
	xmkdir(dirPath);
	FILE *bwt_fp = xget_file(dirPath, ".bwt", "wb");
	FILE *sa_fp = xget_file(dirPath, ".sa", "wb");
	FILE *acg_fp = xget_file(dirPath, ".acg", "wb");//ACGT counter file

	//here combine occ and bwt_str,
	uint32_t block_size = BWT_BYTE_PER_BLOCK + 40;
	uint32_t bufferSize = block_size << 8;//(128 + 40) * 256 ~ 42K

	uint8_t *buffer = (uint8_t*)xmalloc(bufferSize * sizeof(uint8_t));// just ensure that buffer is long enough acutally it won't surpass buffersize + 2
	uint64_t block_num = (bt->len_bwt_str + BP_PER_BLOCK - 1)/BP_PER_BLOCK;//the total number of blocks
	///part one: FM index length
	uint64_t byteLen = block_num * block_size;//length of FM index, in byte
	fwrite(&byteLen, 8, 1, bwt_fp);
	uint64_t L_bwt_occ = ((bt->len_bwt_str + 1) >> 1);
	uint64_t bufferPoint = 0, occPoint = 0, bwtPoint = 0;
	//part two FM index
	for(uint64_t i = 0; i < block_num; i++)
	{
		//store the base part
		memcpy((uint64_t *)(buffer + bufferPoint),bt->occCheck + occPoint, sizeof(uint64_t)*5);
		occPoint += 5;
		bufferPoint += 40;
		//store the BWT part
		if(i < block_num - 1)
			memcpy(buffer + bufferPoint, bt->bin_bwt + bwtPoint, BWT_BYTE_PER_BLOCK); //for the first (n - 1) blocks
		else
			memcpy(buffer + bufferPoint, bt->bin_bwt + bwtPoint, L_bwt_occ - bwtPoint);//for the last block
		bwtPoint += BWT_BYTE_PER_BLOCK;
		bufferPoint += BWT_BYTE_PER_BLOCK;
		//store in hard disk
		if(bufferPoint >= bufferSize || (i == block_num - 1))
		{
			fwrite(buffer, sizeof(uint8_t), bufferPoint, bwt_fp);
			bufferPoint = 0;
		}
	}
	//part three: rank
	fwrite(bt->rank, sizeof(uint64_t),5,bwt_fp);
	///part four: hash index: 1G
	uint64_t hash_index_size = (0x1l << (L_PRE_IDX<<1)) + 1;
	fwrite(bt->hash_index,sizeof(uint64_t),hash_index_size,bwt_fp);
	///part five: ACGT counter size
	uint64_t AGCTCounterSize = 1 << 16;
	fwrite(&AGCTCounterSize, sizeof(uint64_t), 1, acg_fp);
	///part six: AGCTCounter
	for(int i = 0; i < 5; i++)
		fwrite(bt->AGCTCounter[i], sizeof(uint8_t), AGCTCounterSize, acg_fp);
	///part seven: SA part size
	fprintf(stderr,"SA texst2\n");
	fwrite(&(bt->sa_size), sizeof(uint64_t), 1, sa_fp);
	//part eight: SA part
	fwrite(bt->sa_taxon, sizeof(SA_taxon), bt->sa_size, sa_fp);

	fclose(bwt_fp);
	fclose(sa_fp);
	fclose(acg_fp);
	free(buffer);//free todo: end of program, no need to free
	return 0;
}

void build_BWT(bwt *bt, char *bwt_string, uint64_t bwt_len, char *dirPath){
	fprintf(stderr,"Start building BWT index....\n");
	bt->bwt_str = bwt_string;
	bt->len_bwt_str = bwt_len;
	bwt_cal_check_point(bt);//store check point, (FM index, "base" part), checked
	bwt_str2bwt_occ(bt);//bin BWT string
	bwt_cal_AGCTCounter(bt);//ACGT counter
}

