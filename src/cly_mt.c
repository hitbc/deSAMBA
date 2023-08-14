/*
 * classify_main.c
 *
 *  Created on: 2018-5-14
 *      Author: fenghe
 */
#include <string.h>
#include "lib/desc.h"
#include <getopt.h>
#include "lib/utils.h"
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "lib/kthread.h"
#include "zlib.h"

#include "bwt.h"
#include "cly.h"
#define N_NEEDED 5000
#define MAX_read_size 10000000//10M
uint64_t total_sequences = 0;

typedef struct
{
	kseq_t	*	seqs;
	cly_r   *	results;
	long int 	readNum;
	void 	*	share_data_pointer;
}CLASSIFY_THREAD_DATA;

typedef struct
{
	//shared
	DA_IDX 				*idx;
	kstream_t 			*_fp;
	MAP_opt 			 *o;
	//for each pipeline thread
	Classify_buff_pool 	*buff;
	CLASSIFY_THREAD_DATA *data;
}CLASSIFY_SHARE_DATA;

int read_reads(kstream_t *_fp, kseq_t *_seqs, int n_needed)
{
	kseq_t *temp = _seqs;
	int i, rst = 0, total_length = 0;
	for( i = 0; i < n_needed &&  total_length < MAX_read_size; ++i)
	{
		temp[i].f = _fp;
		rst = kseq_read(temp+i);
		if(rst < 0)
			break;
		total_length += temp[i].seq.l;
	}
	total_sequences += i;
	return i;
}

static char primary_string[3][4] = { "PRI", "SEC",	"SUP"};

static inline void print_hit(chain_item *c, REF_INFO *r_i, int rst_cnt, FILE* outfile)
{
	fprintf(outfile,
			"%3d "
			"%s "
			"%s "
			"%20s "

			"ts:%-10d "
			"te:%-10d "
			"qs:%-10d "
			"qe:%-10d "

			//debug value:
			//"%6d "
			//"%5d "
			"%-5d\t"
			"%d\t"
			//"%d\t"
			//"%d\t"
			//"%s\t"
			//"ln:%-10d\t"
			//"%d\t"
			"\n",

			rst_cnt,
			primary_string[c->primary - 1],
			(c->direction)?"F":"R",
			r_i[c->ref_ID].ref_name,

			c->t_st,
			c->t_ed,
			c->q_st,
			c->q_ed,
			//debug value:
			//c->ref_ID,
			//c->chain_id,
			c->sum_score,
			c->indel
			//(c->t_ed - c->t_st) - (c->q_ed - c->q_st)
			//c->anchor_number,
			//(c->with_top_anchor)?"WT":"NT",
			//c->q_ed - c->q_st,
			//c->pri_index
	);
}

void print_anchor(Anchor* anchor_b, REF_INFO *r_i, FILE* outfile)
{
	//basic part
	fprintf(outfile,
			"%c " 			//dir
			"%20s\t"			//ref ID
			"%5d\t"			//index in read
			"%10u\t",		//ref-offset
			(anchor_b->direction)?'F':'R',
			r_i[anchor_b->ref_ID].ref_name,
			anchor_b->index_in_read,
			anchor_b->ref_offset);

	//debug_part
	if(DEBUG)
		fprintf(outfile,
			"%4d\t"			//seed ID
			"%4d\t"			//chain ID
			"%c "			//useless
			"%d\t"
			"%d\t"
			"%d\t"
			"%d\t"
			"%c ",
			anchor_b->seed_ID,
			anchor_b->chain_id,
			(anchor_b->anchor_useless)?'*':' ',
			anchor_b->a_m.left_len,
			anchor_b->a_m.left_ED,
			anchor_b->a_m.rigt_len,
			anchor_b->a_m.rigt_ED,
			(anchor_b->duplicate)?'D':'N'
	);
	//am part
	fprintf(outfile,
			"%d\t"
			"%d\t"
			"\n",
			anchor_b->a_m.score,
			anchor_b->a_m.mtch_len);//,
}

static int inline cmp_anchor(const void * a_, const void * b_)
{
	Anchor *a = (Anchor *)a_, *b = (Anchor *)b_;
	if(a->chain_id != b->chain_id)
		return a->chain_id - b->chain_id;
	return a->index_in_read - b->index_in_read;
}

#define SHOW_ANCHOR_WHEN_NO_RST 0
void output_one_result_des(DA_IDX *idx, cly_r *p_rst, MAP_opt *o)
{
	REF_INFO *r_i = idx->r_i_v.a;
	//NAME:
	fprintf(o->outfile,
			"%s\t"
			"%s\t"
			"%s\t"
			"%ld\t"
			"n_rst:[%ld]\t"
			"n_anc:[%ld]\t"
			"\n",
			p_rst->read->name.s,
			(p_rst->hit.n)?"CLASSIFY":"UNCLASSIFY",
			(p_rst->fast_classify)?"FAST":"SLOW",
			p_rst->read->seq.l,
			p_rst->hit.n,
			p_rst->anchor_v.n
	);
	//RESULT
	chain_item * c_s = p_rst->hit.a, *c_e = c_s + p_rst->hit.n;
	int rst_cnt = 0;
	for(chain_item *c = c_s; c < c_e; c++)//primary + sup
		if(c->pri_index == 0)
			print_hit(c, r_i, rst_cnt++, o->outfile);
	for(chain_item *c = c_s; c < c_e; c++)//secondary
		if(c->pri_index > 0 && c->pri_index <= o->max_sec_N)
			print_hit(c, r_i, rst_cnt++, o->outfile);
	//ANCHORs
	//sort anchor
	if(o->show_anchor == true || (rst_cnt == 0 && SHOW_ANCHOR_WHEN_NO_RST))//when result == 0, show anchors
	{
		qsort(p_rst->anchor_v.a, p_rst->anchor_v.n, sizeof(Anchor), cmp_anchor);
		Anchor* anchor_b = p_rst->anchor_v.a, *anchor_e = anchor_b + p_rst->anchor_v.n;
		for(; anchor_b < anchor_e; anchor_b++)
		{
			if(anchor_b->a_m.score == 0)
				continue;
			print_anchor(anchor_b, r_i, o->outfile);
		}
	}
	fprintf(o->outfile, "\n");
}

void output_one_result_full(DA_IDX *idx, cly_r *p_rst, MAP_opt *o)
{
	REF_INFO *r_i = idx->r_i_v.a;
	//NAME:
	fprintf(o->outfile,
			"%s\t"
			"%s\t"
			"%s\t"
			"%ld\t"
			"n_rst:[%ld]\t"
			"n_anc:[%ld]\t"
			"\n",
			p_rst->read->name.s,
			(p_rst->hit.n)?"CLASSIFY":"UNCLASSIFY",
			(p_rst->fast_classify)?"FAST":"SLOW",
			p_rst->read->seq.l,
			p_rst->hit.n,
			p_rst->anchor_v.n
	);
	//RESULT
	chain_item * c_s = p_rst->hit.a, *c_e = c_s + p_rst->hit.n;
	int rst_cnt = 0;
	for(chain_item *c = c_s; c < c_e; c++)//primary + sup
		if(c->pri_index == 0)
			print_hit(c, r_i, rst_cnt++, o->outfile);
	for(chain_item *c = c_s; c < c_e; c++)//secondary
		if(c->pri_index > 0)
			print_hit(c, r_i, rst_cnt++, o->outfile);
	//ANCHORs
	if((o->show_anchor == true && rst_cnt == 0))//when rst_cnt == 0; always show anchors
	{
		qsort(p_rst->anchor_v.a, p_rst->anchor_v.n, sizeof(Anchor), cmp_anchor);	//sort anchor
		Anchor* anchor_b = p_rst->anchor_v.a, *anchor_e = anchor_b + p_rst->anchor_v.n;
		for(; anchor_b < anchor_e; anchor_b++)
		{
			if(anchor_b->a_m.score == 0)
				continue;
			print_anchor(anchor_b, r_i, o->outfile);
		}
	}
	fprintf(o->outfile, "\n");
}

void output_one_result_sam(DA_IDX *idx, cly_r *p_rst, int output_seq,  MAP_opt *o)
{
	REF_INFO *r_i = idx->r_i_v.a;
	char star[2] = "*", *seq_s = (output_seq)?p_rst->read->seq.s:star, *qual_s = (output_seq)?p_rst->read->qual.s:star;
	//unmapped
	if(p_rst->hit.n == 0)
	{		//NAME:
		fprintf(o->outfile, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t"
				"%s\t"
				"%s\t"
				"\n",
				p_rst->read->name.s,
				seq_s,
				qual_s
		);
		return;
	}
	//mapped
	//primary
	uint32_t read_l = p_rst->read->seq.l;
	chain_item * c_s = p_rst->hit.a, *c_e = c_s + p_rst->hit.n;
	int flag = c_s->direction?0:0x10;//direction
	int mapQ_PRI = 0;
	if(p_rst->hit.n == 1 ||(c_s->sum_score - c_s[1].sum_score > 5))
		mapQ_PRI = 30;
	else
		mapQ_PRI = (c_s->sum_score - c_s[1].sum_score)<<2;
	fprintf(o->outfile,
			"%s\t"
			"%d\t"
			"%s\t"
			"%d\t"
			"%d\t"
			"%dS%dM%dS\t"//CIGAR
			"*\t0\t0\t"
			"%s\t"
			"%s\t"
			"AS:i:%d\t"//number of mapping 9-mers, used as score
			//"di:i:%d\t"//sum number of big deletion and insertions(only an es30timated value)
			"\n",
			p_rst->read->name.s,
			flag,
			r_i[c_s->ref_ID].ref_name,
			c_s->t_st,
			mapQ_PRI,
			c_s->q_st, c_s->q_ed - c_s->q_st, read_l - c_s->q_ed,
			seq_s,
			qual_s,
			c_s->sum_score
			//c_s->indel
	);
	//supplementary and secondary
	int rst_cnt = 0;
	for(int loop = 0; loop <= 1; loop++)
	{
		for(chain_item *c = c_s + 1; c < c_e; c++)//sup + sec
		{
			int show_rst = false;
			int flag = c->direction?0:0x10;//direction
			int mapQ = 0;
			if(loop == 0 && c->pri_index == 0)//supplementary
			{
				show_rst = true;
				flag += 0x800;//supplementary alignment
				mapQ = MIN(30, mapQ_PRI);
			}
			else if(loop == 1 && c->pri_index > 0 && c->pri_index <= o->max_sec_N)//secondary
			{
				show_rst = true;
				flag += 0x100;//secondary alignment
			}
			if(show_rst == true)
			{
				rst_cnt++;
				fprintf(o->outfile,
					"%s\t"
					"%d\t"
					"%s\t"
					"%d\t"
					"%d\t"
					"%d%c%dM%d%c\t"//CIGAR
					"*\t0\t0\t"
					"*\t"
					"*\t"
					"AS:i:%d\t"//number of mapping 9-mers, used as score
					//"di:i:%d\t"//sum number of big deletion and insertions(only an estimated value)
					"\n",
					p_rst->read->name.s,
					flag,
					r_i[c->ref_ID].ref_name,
					c->t_st,
					mapQ,
					c->q_st, (loop == 0)?'H':'S', c->q_ed - c->q_st, read_l - c->q_ed, (loop == 0)?'H':'S',
					c->sum_score
					//c->indel
				);
			}
		}
	}
}

#define OUTPUT_MODE_SAM 	 1
#define OUTPUT_MODE_SAM_FULL 2
#define OUTPUT_MODE_DES 	 3
#define OUTPUT_MODE_DES_FULL 4
void output_results(DA_IDX *idx, cly_r *results, long int n_results, MAP_opt *o)
{
	cly_r *p_rst = results, *p_e_rst = p_rst + n_results;
	if(o->out_format == OUTPUT_MODE_SAM)
		for(; p_rst < p_e_rst; ++p_rst)//out put SAM results ,but with out read sequence
			output_one_result_sam(idx, p_rst, false, o);//the last 0 occupy the MAPQ part
	else if(o->out_format == OUTPUT_MODE_SAM_FULL)
		for(; p_rst < p_e_rst; ++p_rst)//out put SAM results
			output_one_result_sam(idx, p_rst, true, o);
	else if(o->out_format == OUTPUT_MODE_DES)
		for(; p_rst < p_e_rst; ++p_rst)//out put DES results
			output_one_result_des(idx, p_rst, o);//the last 0 occupy the MAPQ part
	else if(o->out_format == OUTPUT_MODE_DES_FULL)
		for(; p_rst < p_e_rst; ++p_rst)//out put FULL results
			output_one_result_full(idx, p_rst, o);//the last 0 occupy the MAPQ part
}

extern void classify_seq(kseq_t *read, DA_IDX* idx, cly_r *results, Classify_buff_pool *buff);

static void inline worker_for(void *_data, long data_index, int thread_index)
{ // kt_for() callback
	CLASSIFY_THREAD_DATA *d = (CLASSIFY_THREAD_DATA *) _data;
	CLASSIFY_SHARE_DATA  *s = (CLASSIFY_SHARE_DATA  *) (d->share_data_pointer);
	classify_seq(d->seqs + data_index,
			s->idx, d->results + data_index,
			s->buff + thread_index);
}

void *classify_pipeline(void *shared, int step, int tid, void *_data)
{
	CLASSIFY_SHARE_DATA * s = (CLASSIFY_SHARE_DATA*) shared;
	//step0: read read data from files; step1: process; step2: output result
	if 		(step == 0)
	{
		if((s->data[tid].readNum = read_reads(s->_fp, s->data[tid].seqs, N_NEEDED)))
			return (void *)1;
	}
	else if (step == 1)
	{
		kt_for(s->o->thread_num, worker_for, s->data + tid, s->data[tid].readNum);
		return (void *)1;
	}
	else if (step == 2)
	{
		output_results(s->idx, s->data[tid].results, s->data[tid].readNum, s->o);
		return (void *)1;
	}
	return 0;
}

//function usage
/**
 *  int Q_MEM[Q_MEM_MAX];
 *  int Q_LV[MAX_LV_WRONG][MAX_LV_R_LEN];
 *  double P_E = 0.15;
 *  uint64_t L_REF = (0x1<<30);
 *  L_REF *= 10;
 *
 *  calculate_MAPQ_TABLE(Q_MEM, Q_LV, P_E, L_REF);
 *  int (*Q_CONFLICT)[Q_CONF_MAX][Q_CONF_MAX][3] = (int(*)[Q_CONF_MAX][Q_CONF_MAX][3])malloc(sizeof(int)*Q_CONF_MAX*Q_CONF_MAX*Q_CONF_MAX*3);
 *  calculate_ANCHOR_CONFLICT_MAPQ_TABLE(Q_CONFLICT);
 *
 */
void calculate_MAPQ_TABLE(
		int * Q_MEM,
		int (* Q_LV)[MAX_LV_R_LEN],
		double P_E, uint64_t L_REF)
{
	//step0: scores
	double REF_SIZE_PUNALTY = -10 * log(L_REF) 			/ log(10);
	double MATCH_SCORE 		= -10 * log(0.25/(1 - P_E)) / log(10);
	double MISMATCH_PUNALTY = -10 * log(0.75/(P_E)) 	/ log(10);
	//step1: get Q_MEM
	for(int i = 0; i < Q_MEM_MAX; i++)
		Q_MEM[i] = REF_SIZE_PUNALTY + i * MATCH_SCORE + 0.5;

	//step2: get Q_LV
	for(int j = 0; j < MAX_LV_R_LEN; j++)
	{
		for(int i = 0; i < MAX_LV_WRONG; i++)
		{
			Q_LV[i][j] = (j-i) * MATCH_SCORE + i * MISMATCH_PUNALTY + 0.5;
			if(j < 5)
				Q_LV[i][j] += 15;//gain score when ref_l is too small
			Q_LV[i][j] = MAX(Q_LV[i][j], -8);
		}
	}
}

void report_stats(struct timeval start)
{
	  double seconds = realduration(start);
	  fprintf(stderr, "%ld sequences processed in %.3fs (%.1f Kseq/m).\n",
			  total_sequences,
			  seconds,
			  total_sequences / 1.0e3 /(seconds / 60));
}

static void classify_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:   %s\n", PACKAGE_NAME);
	fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact:   %s\n\n", CONTACT);
	fprintf(stderr, "  Usage:     %s  classify  [Options] <IndexDir> [ReadFiles.fa][...]>\n", PACKAGE_NAME);
	fprintf(stderr, "  Basic:   \n");
	fprintf(stderr, "    <IndexDir>      FOLDER   the directory contains %s index\n", PACKAGE_NAME);
	fprintf(stderr, "    [ReadFiles.fa]  FILES    reads files, FASTQ(A) format, separated by space\n");
	fprintf(stderr, "  Options:\n");
	fprintf(stderr, "    -h,             help\n");
	//fprintf(stderr, "    -u,             run in \"Strain mode\": this mode will increase strain\n");
	//fprintf(stderr, "                    level classification sensitivity and accuracy, but the speed\n");
	//fprintf(stderr, "                    may be slower.\n");
	fprintf(stderr, "    -t, INT         number of threads[4]\n");
	fprintf(stderr, "    -l, INT         minimum matching length, ignored for NGS reads [170]\n");
	fprintf(stderr, "    -r, INT         max Output number of secondary alignments[5]\n");
	fprintf(stderr, "    -o, FILE        output results into file [stdout]\n");
	fprintf(stderr, "    -s, INT         MIN score[64]\n");
	fprintf(stderr, "    -f, STR         output format, one of:\n");
	fprintf(stderr, "                    - SAM: SAM-like results without SEQ and QUAL and header, default\n");
	fprintf(stderr, "                    - SAM_FULL: SAM-like results with SEQ and QUAL\n");
	fprintf(stderr, "                    - DES: smallest format\n");
	fprintf(stderr, "                    - DES_FULL: all results are showed, ignore '-r' opinion\n");
	fprintf(stderr, "\n");
	//fprintf(stderr, "Options:   -e,           Error rate[0.85]\n");
	//fprintf(stderr, "Options:   -a,           Output all results\n");
}

#define PIPELINE_T_NUM 3//one for reading; one for classifying; one for writing
#define STEP_NUM PIPELINE_T_NUM
kvec_T(kstring_t, kstring_V)

int classify_main(int argc, char *argv[])
{//only use one thread, random means random seed
	double P_E = 0.15;//ERROR rate
	int c = 0;
	MAP_opt o = {170, 4, 5, OUTPUT_MODE_SAM, false, stdout, 64};//, false};
	while((c = getopt(argc, argv, "ht:l:r:f:o:s:")) >= 0)
	{
		     if(c == 'h') {classify_usage(); return 0;}
		else if(c == 't') o.thread_num = (int)atoi(optarg);
		//else if(c == 'u') o.strain_mode = true;
		else if(c == 'l') o.L_min_matching = (int)atoi(optarg);
		else if(c == 'r') o.max_sec_N = (int)atoi(optarg);
		else if(c == 'o') o.outfile = xopen(optarg, "w");
		else if(c == 's') o.min_score = (int)atoi(optarg);
		else if(c == 'f')
		{
				 if(strcmp(optarg, "SAM") == 0) 		o.out_format = OUTPUT_MODE_SAM;
			else if(strcmp(optarg, "SAM_FULL") == 0) o.out_format = OUTPUT_MODE_SAM_FULL;
			else if(strcmp(optarg, "DES") == 0) 		o.out_format = OUTPUT_MODE_DES;
			else if(strcmp(optarg, "DES_FULL") == 0) o.out_format = OUTPUT_MODE_DES_FULL;
		}
	}
	if (optind + 2 > argc)
	{
		classify_usage();
		return 0;
	}
	char * index_dir = argv[optind++];
	kstring_V reads_files = {0};
	while (optind < argc)
	{
		kstring_t *reads_file;
		kv_pushp(kstring_t, reads_files, &reads_file);
		kstring_initp(reads_file);
		kputs(argv[optind++], reads_file);
	}
	fprintf(stderr,"loading index\t");	//load index
	DA_IDX idx = {0};
	load_idx(&idx, index_dir);
	idx.filter_min_length = o.L_min_matching;
	idx.filter_min_score = o.min_score;
	idx.filter_min_score_LV3 = o.min_score + 10;
	//idx.strain_mode = o.strain_mode;
	idx.mapQ.Q_MEM = xmalloc_t(int, Q_MEM_MAX);
	idx.mapQ.Q_LV  = (int (*)[MAX_LV_R_LEN])xmalloc(MAX_LV_WRONG*MAX_LV_R_LEN*sizeof(int));
	calculate_MAPQ_TABLE(idx.mapQ.Q_MEM, idx.mapQ.Q_LV, P_E, idx.ref_bin.n*4);

	struct timeval start;
	gettimeofday(&start, NULL);
	fprintf(stderr,"Start classify\n");
	double cpu_time = cputime();

	CLASSIFY_THREAD_DATA data[PIPELINE_T_NUM];
	o.thread_num = MAX(o.thread_num,1);
	Classify_buff_pool *buff = xcalloc_t(Classify_buff_pool, o.thread_num);
	//init: buff
	for(int i = 0; i < o.thread_num; i++)
	{
		buff[i].sa_hash[0] = xmalloc(sizeof(sparse_align_HASH)*0x100000);//1M for each
		buff[i].sa_hash[1] = xmalloc(sizeof(sparse_align_HASH)*0x100000);//1M for each
	}
	CLASSIFY_SHARE_DATA share = {&idx, NULL, &o, buff, data};

	for(int i = 0; i < PIPELINE_T_NUM; i++)
	{
		data[i].seqs	= xcalloc_t(kseq_t, N_NEEDED);
		data[i].results	= xcalloc_t(cly_r, N_NEEDED);
		data[i].share_data_pointer = &share;
	}
	for (uint32_t i = 0; i < reads_files.n; ++i)
	{
		gzFile fp = xzopen(reads_files.a[i].s, "r");
		share._fp = ks_init(fp);
		fprintf(stderr, "Processing file: [%s].\n", reads_files.a[i].s);
		kt_pipeline(PIPELINE_T_NUM, classify_pipeline, &share, STEP_NUM);
		gzclose(fp);
	}
	report_stats(start);
	fprintf(stderr, "Classify CPU: %.3f sec\n", cputime() - cpu_time);
	return 0;
}


