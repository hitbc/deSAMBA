//jellyfish count -m 31 -s 10000000 -t 4  83263.fasta

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>

#include "lib/desc.h"
#include "lib/utils.h"

#define L_BUFF (1<<20)//1M
#define L_B 12 //BUCKET LENGTH, need less than KMER_LENGTH_PlusOne
#define B_CAP (1<<(L_B<<1))//BUCKET CAPACITY, 16M
#define n_thread 32 //should not surpass 2^8-1
#define STKSIZ 40

typedef struct
{
	size_t 		h_size;//header_size
	FILE * 		file;
	uint64_t 	key_bits;
	uint64_t 	key_len;
	uint64_t 	val_len;
	uint64_t 	key_ct;
}JReader;

void jreader_open(JReader * jr, char* filename_str)
{
	jr->file = xopen(filename_str, "rb");

	//jelly fish header
	char JF_FILE_T[10] = "JFLISTDN";//file type
	char JF_FILE_H[10];//file header
	xread(JF_FILE_H, 1, strlen(JF_FILE_T), jr->file);
	xassert(strncmp(JF_FILE_H, JF_FILE_T, strlen(JF_FILE_T)) == 0," not in proper format\n");
	//jelly fish information
	xread(&(jr->key_bits), sizeof(uint64_t), 1, jr->file);
	xread(&(jr->val_len), sizeof(uint64_t), 1, jr->file);
	fseek(jr->file, 48, SEEK_SET);
	xread(&(jr->key_ct), sizeof(uint64_t), 1, jr->file);
	xassert(jr->val_len == 4, "can only handle 4 byte DB values\n");

	jr->key_len = jr->key_bits / 8 + !! (jr->key_bits % 8);//"!!"means: 0->0; others->1
	jr->h_size = 72 + 2 * (4 + 8 * jr->key_bits);
}

//return the upper bound
uint64_t b_Search(uint64_t *data, uint64_t key, int64_t low, int64_t up)
{
 	 int64_t mid = 0;
 	 while(low <= up)
 	 {
		mid = (low + up)>>1;
		if(key < data[mid])
			up  =  mid - 1;
		else if(key > data[mid])
			low =  mid + 1;
		else
			return mid;
  	}
 	return low;
}

int cmp_uint64(const void *a, const void *b)
{
	uint64_t va = *(uint64_t *) a, vb= *(uint64_t *) b;
	if(va <  vb)
		return -1;
	else if(va == vb)
		return 0;
	return 1;
}

typedef struct
{
	uint64_t up_seg;
	uint64_t low_seg;
	uint64_t *count;
	uint64_t *hashKmer;
	uint64_t totalKmerNum;
}MT_SORT;

void *multiThreadSort(void *data)
{
	MT_SORT * d = (MT_SORT*)data;
	uint64_t low = d->low_seg;
	uint64_t up  = d->up_seg;
	uint64_t *hashKmer = d->hashKmer;
	uint64_t *count = d->count;
	uint64_t i;
	for(i = low; i < up; i++)
		qsort(hashKmer + count[i], count[i + 1]-count[i], sizeof(uint64_t), cmp_uint64);
	if(i == B_CAP - 1)
		qsort(hashKmer + count[i], d->totalKmerNum - count[i],   sizeof(uint64_t), cmp_uint64);
	return NULL;
}

void sort_kmers(uint8_t l_kmer, char *jf_fn ,char *out_file)
{
	//STEP1:open jelly fish file
	JReader jr = {0};
	jreader_open(&jr, jf_fn);
	uint64_t val_len = jr.val_len;
	uint64_t key_ct = jr.key_ct;
	uint64_t l_jf_block = jr.key_len;//length of jelly fish block
	uint64_t l_p = l_jf_block + val_len;//pair size
	char pair[l_p];
	fseek(jr.file, jr.h_size, SEEK_SET);
	fprintf(stderr,"kmer number:%ld\n",key_ct);

	///STEP2:KMER-counting, create a copy of the offsets array
	///for use as insertion positions, and calculate "count" and store kmer in hard disk
	fprintf(stderr,"start KMER sorting\n");
	uint64_t *count = xcalloc_t(uint64_t, B_CAP);//16M*8 = 128M
	uint64_t *w_b = xcalloc_t(uint64_t, L_BUFF);//write Buff, 1M*8 = 8M
	uint64_t idx_b = 0;//buff index
	const char *kmerPath = out_file;
	FILE *fpKmer = xopen(kmerPath ,"wb");
	uint64_t MOVE = (l_kmer - L_B) << 1; //32 defines kmer_plus_one
	uint64_t idx_key = key_ct;
	while(idx_key--)
	{
		xread(pair, 1, l_p, jr.file);
		uint64_t kmer = 0;
		memcpy(&kmer, pair, l_jf_block);
		++count[kmer >> MOVE];//store kmer
		w_b[idx_b++] = kmer;
#ifdef CONSIDER_BOTH_ORIENTATION//store RC kmer;
		uint64_t rcKmer = getRcKmer(kmer, l_kmer);
		++count[rcKmer >> MOVE];
		w_b[idx_b++] = rcKmer;
		xassert(rcKmer < kmer, "'-C' option should be set in Jellyfish.");
#endif
		if(idx_b >= L_BUFF)
		{
			fwrite(w_b,sizeof(uint64_t),L_BUFF,fpKmer);
			idx_b=0;
		}
	}
	fwrite(w_b, sizeof(uint64_t), idx_b, fpKmer);
	fclose(jr.file);
	fclose(fpKmer);
	free(w_b);

	///STEP3:get count and segCount, count point to the low bound of each bucket,
	//while segCount used for distributing works for each thread
	for(uint64_t i = 1; i < B_CAP; i++)
		count[i] += count[i-1];
	xassert(key_ct == count[B_CAP-1], "Wrong in get kmer number!");
	uint64_t segCount[n_thread];
	segCount[0] = 0;
	for(uint64_t i = 1; i < n_thread; i++)
	{
		uint64_t locateNum = i * (key_ct / n_thread);
		uint64_t locate = b_Search(count, locateNum, 0, B_CAP - 1);
		segCount[i] = locate;
	}

	//STEP4:KMER distributing
	fprintf(stderr,"start KMER distributing\n");
	uint64_t *r_b = xcalloc_t(uint64_t, L_BUFF);//read buff
	uint64_t *hashKmer = xcalloc_t(uint64_t, key_ct);
	fpKmer   = xopen(kmerPath ,"rb");
	uint64_t true_r;//true read number
	while((true_r = fread(r_b,sizeof(uint64_t),L_BUFF,fpKmer))>0)
	{
		for(uint64_t i = 0; i < true_r; ++i)
		{
			uint64_t hash_pos = --(count[r_b[i] >> MOVE]);//because it is upper bound(cannot reach)
			hashKmer[hash_pos] = r_b[i];
		}
	}
	fclose(fpKmer);

	//STEP5:KMER MT sorting
	fprintf(stderr,"start KMER sorting\n");
	pthread_t t_ID[n_thread];
	MT_SORT data[n_thread];
	for(uint64_t i = 0; i < n_thread; i++)
	{
		data[i].count = count;
		data[i].hashKmer = hashKmer;
		data[i].low_seg = segCount[i];
		data[i].up_seg = (i < n_thread-1)?segCount[i + 1]:B_CAP - 1;
		data[i].totalKmerNum = key_ct;
	}
	for(uint64_t i=0;i<n_thread;i++)
		xassert(pthread_create( &t_ID[i], NULL, multiThreadSort, data + i) == 0,
				"Error - pthread_create()");
	for(uint64_t i=0;i<n_thread;i++)
		pthread_join(t_ID[i], NULL);

	///store data to disk
	fpKmer = xopen(kmerPath,"wb");
	fwrite(&key_ct, sizeof(uint64_t), 1, fpKmer);
	fwrite(hashKmer,sizeof(uint64_t),key_ct,fpKmer);
	fclose(fpKmer);
	free(r_b);
	free(hashKmer);
	free(count);
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: %s\n", PACKAGE_NAME);
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: %s\n\n", CONTACT);
	fprintf(stderr, "  Usage:   %s kmersort <options> [JF Kmer file]\n\n", PACKAGE_NAME);
	fprintf(stderr, "  Basic:\n");
	fprintf(stderr, "    [JF Kmer file] FILE   Unsorted kmer file generated by Jellyfish\n");
	fprintf(stderr, "  Options:\n");
	fprintf(stderr, "       -k INT        length of kmer in base pair[31]\n");
	fprintf(stderr, "       -o FILE       output sorted kmer file[\"kmer.srt\"]\n");
	fprintf(stderr, "       -h            help\n");
	fprintf(stderr, "\n");
	return 0;
}

int sort_kmers_main(int argc, char *argv[])
{
	uint8_t l_kmer = 31;
	char *out_file = "kmer.srt";
	int opin;
	while((opin = getopt(argc, argv, "k:o:h"))>=0)
	{
		switch(opin)
		{
			case 'o': out_file = optarg; break;
		    case 'k': l_kmer = (uint8_t)atoi(optarg);break;
		    case 'h': return usage(); break;
		    default:  fprintf(stderr,"inappropriate parameters\n");	return usage();
		}
	}
	if(optind + 1 > argc)
	{
		fprintf(stderr,"[JF Kmer file] missed.\n");
		return usage();
	}
	char *jf_fn = argv[optind];//jelly fish file name

	struct timeval start;
	gettimeofday(&start, NULL);
	sort_kmers(l_kmer, jf_fn, out_file);
	fprintf(stderr,"sort read time is [%f], CPU: [%.3f] sec\n", realduration(start), cputime());
	return 0;

}

