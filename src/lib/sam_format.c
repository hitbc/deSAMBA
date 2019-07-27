/**

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdint.h>

#include "utils.h"
#include "desc.h"


// * result format used by minimap2:SAM
// * use function  mm_write_sam_hdr and mm_write_sam2 to write head line and sam line

#define MM_F_NO_SELF     0x001
#define MM_F_AVA         0x002
#define MM_F_CIGAR       0x004
#define MM_F_OUT_SAM     0x008
#define MM_F_NO_QUAL     0x010
#define MM_F_OUT_CG      0x020
#define MM_F_OUT_CS      0x040
#define MM_F_SPLICE      0x080 // splice mode
#define MM_F_SPLICE_FOR  0x100 // match GT-AG
#define MM_F_SPLICE_REV  0x200 // match CT-AC, the reverse complement of GT-AG
#define MM_F_NO_LJOIN    0x400
#define MM_F_OUT_CS_LONG 0x800
#define MM_F_SR          0x1000
#define MM_F_FRAG_MODE   0x2000
#define MM_F_NO_PRINT_2ND  0x4000
#define MM_F_2_IO_THREADS  0x8000
#define MM_F_LONG_CIGAR    0x10000
#define MM_F_INDEPEND_SEG  0x20000
#define MM_F_SPLICE_FLANK  0x40000
#define MM_F_SOFTCLIP      0x80000
#define MM_F_FOR_ONLY      0x100000
#define MM_F_REV_ONLY      0x200000

#define MM_I_HPC          0x1
#define MM_I_NO_SEQ       0x2
#define MM_I_NO_NAME      0x4

#define MM_IDX_MAGIC   "MMI\2"

#define MM_MAX_SEG       255

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static char mm_rg_id[256];

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

unsigned char seq_comp_table[256] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};


int mm_verbose = 1;


typedef struct		 //reference information
{
	char *name;      // name of the db sequence
	uint64_t offset; // offset in mm_idx_t::S
	uint32_t len;    // length
} mm_idx_seq_t;

typedef struct		// minimap2 index
{
	uint32_t n_seq;            // number of reference sequences
	mm_idx_seq_t *seq;         // sequence name, length and offset
	uint32_t *S;               // 4-bit packed reference sequence
} mm_idx_t;

// minimap2 alignment result
typedef struct
{
	uint32_t capacity;          // the capacity of cigar[]??
	int32_t  dp_score;			// DP score;
	int32_t  dp_max;			// score of the max-scoring segment;
	int32_t  dp_max2;  			// score of the best alternate mappings
	uint32_t n_ambi:30, 		// number of ambiguous bases;
			 trans_strand:2; 	// transcript strand: 0 for unknown, 1 for +, 2 for -
	uint32_t n_cigar;           // number of cigar operations in cigar[]
	uint32_t cigar[];
} mm_extra_t;

typedef struct
{
	int32_t id;             // ID for internal uses (see also parent below)
	int32_t cnt;            // number of minimizers; if on the reverse strand
	int32_t rid;            // reference index; if this is an alignment from inversion rescue
	int32_t score;          // DP alignment score
	int32_t qs, qe; 		// query start and end;
	int32_t rs, re;			// reference start and end
	int32_t parent;			// parent==id if primary;
	int32_t subsc;   		// best alternate mapping score
	int32_t as;             // offset in the a[] array (for internal uses only)
	int32_t mlen;     		// seeded exact match length;
	int32_t blen;			// seeded alignment block length
	int32_t n_sub;          // number of sub-optimal mappings
	int32_t score0;         // initial chaining score (before chain merging/spliting)
	uint32_t mapq:8, split:2, rev:1, inv:1, sam_pri:1, proper_frag:1, pe_thru:1, seg_split:1, dummy:16;
	uint32_t hash;
	float div;
	mm_extra_t *p;
} mm_reg1_t;

//minimap read sequence
typedef struct
{
	int l_seq, rid;
	char *name, *seq, *qual;
} mm_bseq1_t;

static inline int mm_qname_len(const char *s)
{
	int l;
	l = strlen(s);
	//delete /1 in q name
	return l >= 3 && s[l-1] >= '0' && s[l-1] <= '9' && s[l-2] == '/'? l - 2 : l;
}

//change string "\t" to true '\t' , and change "\\" to true '\\'
static char *mm_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p)
	{
		if (*p == '\\')
		{
			++p;
			if (*p == 't')
				*q++ = '\t';
			else if (*p == '\\')
				*q++ = '\\';
		}
		else *q++ = *p;
	}
	*q = '\0';
	return s;
}

//check the header @RG line and print to str
static void sam_write_rg_line(kstring_t *str, const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(mm_rg_id, 0, 256);
	if (s == 0)
		return;
	if (strstr(s, "@RG") != s)
	{
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR] the read group line is not started with @RG\n");
		goto err_set_rg;
	}
	if (strstr(s, "\t") != NULL)
	{
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n");
		goto err_set_rg;
	}
	rg_line = strdup(s);
	mm_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0)
	{
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR] no ID within the read group line\n");
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256)
	{
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR] @RG:ID is longer than 255 characters\n");
		goto err_set_rg;
	}
	for (q = p, r = mm_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	sprintf_lite(str, "%s\n", rg_line);

err_set_rg:
	free(rg_line);
}

// *
// * Writing header line,
// *
// * LINE like:
// * @SQ	SN:tid|642493|ref|NC_012665.1	LN:2743
// * @SQ	SN:tid|1814280|ref|NC_031932.1	LN:42702
// * @PG	ID:minimap2	PN:minimap2	VN:2.10-r764-dirty	CL:minimap2 -a -t 8 ./miniMap2_index/meta_ref.mmi SRR/ERR1050068.fastq
// *

void sam_write_hdr(const mm_idx_t *idx, const char *rg, const char *ver, int argc, char *argv[])
{
	kstring_t str = {0,0,0};
	if (idx)
	{
		uint32_t i;
		for (i = 0; i < idx->n_seq; ++i)
			printf(
					"@SQ\t"
					"SN:%s\t"
					"LN:%d\n",
					idx->seq[i].name,
					idx->seq[i].len);
	}
	if (rg)
		sam_write_rg_line(&str, rg);
	sprintf_lite(&str,
			"@PG\t"
			"ID:%s\t"
			"PN:%s",
			PACKAGE_NAME,
			PACKAGE_NAME);//program name
	if (ver)
		sprintf_lite(&str,
				"\tVN:%s",
				PACKAGE_VERSION);//version number
	if (argc > 1)
	{
		int i;
		sprintf_lite(&str,
				"\tCL:%s",
				PACKAGE_NAME);//comment line
		for (i = 1; i < argc; ++i)
			sprintf_lite(&str,
					" "
					"%s",
					argv[i]);
	}
	sprintf_lite(&str,
			"\n");
	fputs(str.s, stdout);//print the string to screen
	free(str.s);
}

//get reference string from index, it is stored in 4-bit format
#define mm_seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)//todo:::change to 2 bit format
int mm_idx_getseq(
		const mm_idx_t *mi,
		uint32_t rid,//reference id
		uint32_t st,
		uint32_t en,
		uint8_t *seq)//output
{
	uint64_t i, st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len)
		return -1;
	if (en > mi->seq[rid].len)
		en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st;
	en1 = mi->seq[rid].offset + en;
	for (i = st1; i < en1; ++i)
		seq[i - st1] = mm_seq4_get(mi->S, i);
	return en - st;
}

static void write_cs(
		kstring_t *s,//output
		const mm_idx_t *mi,//ref
		const mm_bseq1_t *t,
		const mm_reg1_t *r,
		int no_iden)
{
	int i, q_off, t_off;
	uint8_t *qseq, *tseq;
	char *tmp;
	if (r->p == 0)
		return;
	sprintf_lite(s,
			"\tcs:Z:");
	qseq = (uint8_t*)malloc(r->qe - r->qs);
	tseq = (uint8_t*)malloc(r->re - r->rs);
	tmp = (char*)malloc(r->re - r->rs > r->qe - r->qs? r->re - r->rs + 1 : r->qe - r->qs + 1);
	mm_idx_getseq(mi, r->rid, r->rs, r->re, tseq);
	if (!r->rev)
	{
		for (i = r->qs; i < r->qe; ++i)
			qseq[i - r->qs] = seq_nt4_table[(uint8_t)t->seq[i]];
	}
	else
	{
		for (i = r->qs; i < r->qe; ++i)
		{
			uint8_t c = seq_nt4_table[(uint8_t)t->seq[i]];
			qseq[r->qe - i - 1] = c >= 4? 4 : 3 - c;
		}
	}
	for (i = q_off = t_off = 0; i < r->p->n_cigar; ++i)
	{
		int j, op = r->p->cigar[i]&0xf, len = r->p->cigar[i]>>4;
		assert(op >= 0 && op <= 3);
		if (op == 0)
		{
			int l_tmp = 0;
			for (j = 0; j < len; ++j)
			{
				if (qseq[q_off + j] != tseq[t_off + j])
				{
					if (l_tmp > 0)
					{
						if (!no_iden)
						{
							tmp[l_tmp] = 0;
							sprintf_lite(s,
									"=%s",
									tmp);
						}
						else
							sprintf_lite(s,
									":%d",
									l_tmp);
						l_tmp = 0;
					}
					sprintf_lite(s,
							"*%c"
							"%c",
							"acgtn"[tseq[t_off + j]],
							"acgtn"[qseq[q_off + j]]);
				}
				else
					tmp[l_tmp++] = "ACGTN"[qseq[q_off + j]];
			}
			if (l_tmp > 0)
			{
				if (!no_iden)
				{
					tmp[l_tmp] = 0;
					sprintf_lite(s, "=%s", tmp);
				}
				else
					sprintf_lite(s, ":%d", l_tmp);
			}
			q_off += len, t_off += len;
		}
		else if (op == 1)
		{
			for (j = 0, tmp[len] = 0; j < len; ++j)
				tmp[j] = "acgtn"[qseq[q_off + j]];
			sprintf_lite(s, "+%s", tmp);
			q_off += len;
		}
		else if (op == 2)
		{
			for (j = 0, tmp[len] = 0; j < len; ++j)
				tmp[j] = "acgtn"[tseq[t_off + j]];
			sprintf_lite(s, "-%s", tmp);
			t_off += len;
		}
		else
		{
			assert(len >= 2);
			sprintf_lite(s,
					"~"
					"%c"
					"%c"
					"%d"
					"%c"
					"%c",
					"acgtn"[tseq[t_off]],
					"acgtn"[tseq[t_off+1]],
					len,
					"acgtn"[tseq[t_off+len-2]],
					"acgtn"[tseq[t_off+len-1]]);
			t_off += len;
		}
	}
	assert(t_off == r->re - r->rs && q_off == r->qe - r->qs);
	free(qseq);
	free(tseq);
	free(tmp);
}

static inline void write_tags(kstring_t *s, const mm_reg1_t *r)
{
	int type;
	if (r->id == r->parent)
		type = r->inv? 'I' : 'P';
	else
		type = r->inv? 'i' : 'S';
	if (r->p)
	{
		sprintf_lite(s,
				"\tNM:i:%d"
				"\tms:i:%d"
				"\tAS:i:%d"
				"\tnn:i:%d",
				r->blen - r->mlen + r->p->n_ambi,
				r->p->dp_max,
				r->p->dp_score,
				r->p->n_ambi);
		if (r->p->trans_strand == 1 || r->p->trans_strand == 2)
			sprintf_lite(s,
					"\t"
					"ts:A:%c",
					"?+-?"[r->p->trans_strand]);
	}
	sprintf_lite(s,
			"\t"
			"tp:A:%c\t"
			"cm:i:%d\t"
			"s1:i:%d",
			type,
			r->cnt,
			r->score);
	if (r->parent == r->id)
		sprintf_lite(s,
				"\t"
				"s2:i:%d",
				r->subsc);
	if (r->div >= 0.0f && r->div <= 1.0f)
	{
		char buf[8];
		if (r->div == 0.0f)
			buf[0] = '0', buf[1] = 0;
		else
			sprintf(buf, "%.4f", r->div);
		sprintf_lite(s,
				"\t"
				"dv:f:%s",
				buf);
	}
	if (r->split)
		sprintf_lite(s,
				"\t"
				"zd:i:%d",
				r->split);
}

void mm_write_paf(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, void *km, int opt_flag)
{
	s->l = 0;
	sprintf_lite(s,
			"%s\t"
			"%d\t"
			"%d\t"
			"%d\t"
			"%c\t",
			t->name,
			t->l_seq,
			r->qs,
			r->qe,
			"+-"[r->rev]);
	if (mi->seq[r->rid].name)
		sprintf_lite(s, "%s", mi->seq[r->rid].name);
	else
		sprintf_lite(s, "%d", r->rid);
	sprintf_lite(s,
			"\t"
			"%d\t"
			"%d\t"
			"%d",
			mi->seq[r->rid].len,
			r->rs,
			r->re);
	sprintf_lite(s,
			"\t"
			"%d\t"
			"%d",
			r->mlen,
			r->blen);
	sprintf_lite(s,
			"\t%d",
			r->mapq);
	write_tags(s, r);
	if (r->p && (opt_flag & MM_F_OUT_CG))
	{
		uint32_t k;
		sprintf_lite(s,
				"\tcg:Z:");
		for (k = 0; k < r->p->n_cigar; ++k)
			sprintf_lite(s,
					"%d"
					"%c",
					r->p->cigar[k]>>4,
					"MIDN"[r->p->cigar[k]&0xf]);
	}
	if (r->p && (opt_flag & MM_F_OUT_CS))
		write_cs(s, mi, t, r, !(opt_flag&MM_F_OUT_CS_LONG));
}

static void sam_write_sq(kstring_t *s, char *seq, int l, int rev, int comp)
{
	extern unsigned char seq_comp_table[256];
	if (rev)
	{
		int i;
		str_enlarge(s, l);
		for (i = 0; i < l; ++i)
		{
			int c = seq[l - 1 - i];
			s->s[s->l + i] = c < 128 && comp? seq_comp_table[c] : c;
		}
		s->l += l;
	}
	else kstrcpy(s, seq, seq + l);
}

static inline const mm_reg1_t *get_sam_pri(int n_regs, const mm_reg1_t *regs)
{
	int i;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].sam_pri)
			return &regs[i];
	assert(n_regs == 0);
	return NULL;
}

static void write_sam_cigar(
		kstring_t *s,
		int sam_flag,
		int in_tag,
		int qlen,
		const mm_reg1_t *r,
		int opt_flag)
{
	if (r->p == 0)
	{
		sprintf_lite(s, "*");
	}
	else
	{
		uint32_t k, clip_len[2];
		clip_len[0] = r->rev? qlen - r->qe : r->qs;
		clip_len[1] = r->rev? r->qs : qlen - r->qe;
		if (in_tag)
		{
			int clip_char = (sam_flag&0x800) && !(opt_flag&MM_F_SOFTCLIP)? 5 : 4;
			sprintf_lite(s,
					"\t"
					"CG:B:I");
			if (clip_len[0])
				sprintf_lite(s,
						","
						"%u",
						clip_len[0]<<4|clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				sprintf_lite(s,
						","
						"%u",
						r->p->cigar[k]);
			if (clip_len[1])
				sprintf_lite(s,
						","
						"%u",
						clip_len[1]<<4|clip_char);
		}
		else
		{
			int clip_char = (sam_flag&0x800) && !(opt_flag&MM_F_SOFTCLIP)? 'H' : 'S';
			if (clip_len[0])
				sprintf_lite(s,
						"%d"
						"%c",
						clip_len[0],
						clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				sprintf_lite(s,
						"%d"
						"%c",
						r->p->cigar[k]>>4,
						"MIDN"[r->p->cigar[k]&0xf]);
			if (clip_len[1])
				sprintf_lite(s,
						"%d"
						"%c",
						clip_len[1],
						clip_char);
		}
	}
}


// **[PART1 ] ERR1050068.1
// * [PART2 ] 0
// * [PART3 ] tid|186538|ref|NC_002549.1
// * [PART3 ] 14271
// * [PART3 ] 60
// * [PART3 ] 27S11M2D3M1I42M2D1M1D3M3I18M1D1M3I6M21D4M29S
// * [PART3 ] *
//// * [PART3 ] 0
// * [PART3 ] 0
// * [PART3 ] ACTGTCTTGTAAATGCGGGATCGGTCAAACAATATCATGCAGCTAGTCCTAGAGCAGGTGCTTGGAGGGTTATCCTTCTTGAA
// * [PART3 ] ###################################################################################
// * [PART3 ] NM:i:359
 // * [PART3 ] ms:i:792
// * [PART3 ] AS:i:792
// * [PART3 ] nn:i:0
// * [PART3 ] tp:A:P
// * [PART3 ] cm:i:15
// * [PART3 ] s1:i:125
// * [PART3 ] s2:i:0
//// * [PART3 ] dv:f:0.1823

void mm_write_sam2(
		kstring_t *s, //output
		const mm_idx_t *mi,//index, to get ref seq info
		const mm_bseq1_t *t,//read seq
		int seg_idx,
		int reg_idx,
		int n_seg,
		const int *n_regss,
		const mm_reg1_t *const* regss,//result
		int opt_flag)
{
	const int max_bam_cigar_op = 65535;
	int flag, n_regs = n_regss[seg_idx], cigar_in_tag = 0;
	int this_rid = -1, this_pos = -1, this_rev = 0;
	const mm_reg1_t *regs = regss[seg_idx], *r_prev = NULL, *r_next;
	const mm_reg1_t *r = n_regs > 0 && reg_idx < n_regs && reg_idx >= 0? &regs[reg_idx] : NULL;

	// find the primary of the previous and the next segments, if they are mapped
	if (n_seg > 1)
	{
		int i, next_sid = (seg_idx + 1) % n_seg;
		r_next = get_sam_pri(n_regss[next_sid], regss[next_sid]);
		if (n_seg > 2)
		{
			for (i = 1; i <= n_seg - 1; ++i)
			{
				int prev_sid = (seg_idx + n_seg - i) % n_seg;
				if (n_regss[prev_sid] > 0)
				{
					r_prev = get_sam_pri(n_regss[prev_sid], regss[prev_sid]);
					break;
				}
			}
		}
		else
			r_prev = r_next;
	}
	else
		r_prev = r_next = NULL;

	// write QNAME[PART1]
	s->l = 0;
	sprintf_lite(s, "%s", t->name);
	if (n_seg > 1)
		s->l = mm_qname_len(t->name); // trim the suffix like /1 or /2

	// write flag
	flag = n_seg > 1? 0x1 : 0x0;
	if (r == 0)
	{
		flag |= 0x4;
	}
	else
	{
		if (r->rev)
			flag |= 0x10;
		if (r->parent != r->id)
			flag |= 0x100;
		else if (!r->sam_pri)
			flag |= 0x800;
	}
	if (n_seg > 1)
	{
		if (r && r->proper_frag)
			flag |= 0x2; // TODO: this doesn't work when there are more than 2 segments
		if (seg_idx == 0)
			flag |= 0x40;
		else if (seg_idx == n_seg - 1)
			flag |= 0x80;
		if (r_next == NULL)
			flag |= 0x8;
		else if (r_next->rev)
			flag |= 0x20;
	}
	sprintf_lite(s, "\t%d", flag);

	// write coordinate, MAPQ and CIGAR
	if (r == 0)//no result
	{
		if (r_prev)//with previous segment
		{
			this_rid = r_prev->rid, this_pos = r_prev->rs;
			sprintf_lite(s,
					"\t"
					"%s\t"
					"%d\t0\t*", mi->seq[this_rid].name, this_pos+1);
		}
		else
			sprintf_lite(s,
					"\t"
					"*\t"
					"0\t"
					"0\t"
					"*");
	}
	else
	{
		this_rid = r->rid, this_pos = r->rs, this_rev = r->rev;
		//ref name; position; mapQ
		sprintf_lite(s,
				"\t"
				"%s\t"
				"%d\t"
				"%d\t",
				mi->seq[r->rid].name,
				r->rs+1,
				r->mapq);
		if ((opt_flag & MM_F_LONG_CIGAR) && r->p && r->p->n_cigar > max_bam_cigar_op - 2)
		{
			int n_cigar = r->p->n_cigar;
			if (r->qs != 0)
				++n_cigar;
			if (r->qe != t->l_seq)
				++n_cigar;
			if (n_cigar > max_bam_cigar_op)
				cigar_in_tag = 1;
		}
		if (cigar_in_tag)
		{
			if (flag & 0x100)
				sprintf_lite(s, "0S"); // secondary alignment
			else if (flag & 0x800)
				sprintf_lite(s, "%dS", r->re - r->rs); // supplementary alignment
			else
				sprintf_lite(s, "%dS", t->l_seq);
		}
		else
			write_sam_cigar(s, flag, 0, t->l_seq, r, opt_flag);
	}

	// write mate positions
	if (n_seg > 1)
	{
		int tlen = 0;
		if (this_rid >= 0 && r_next)
		{
			if (this_rid == r_next->rid)
			{
				int this_pos5 = r && r->rev? r->re - 1 : this_pos;
				int next_pos5 = r_next->rev? r_next->re - 1 : r_next->rs;
				tlen = next_pos5 - this_pos5;
				sprintf_lite(s, "\t=\t");
			}
			else
				sprintf_lite(s,
						"\t"
						"%s\t",
						mi->seq[r_next->rid].name);
			sprintf_lite(s,
					"%d\t",
					r_next->rs + 1);
		}
		else if (r_next)
		{ // && this_rid < 0
			sprintf_lite(s,
					"\t"
					"%s\t"
					"%d\t",
					mi->seq[r_next->rid].name,
					r_next->rs + 1);
		}
		else if (this_rid >= 0)
		{ // && r_next == NULL
			int this_pos5 = this_rev? r->re - 1 : this_pos; // this_rev is only true when r != NULL
			tlen = this_pos - this_pos5; // next_pos5 will be this_pos
			sprintf_lite(s,
					"\t=\t%d\t",
					this_pos + 1); // next segment will take r's coordinate
		}
		else
			sprintf_lite(s,
					"\t"
					"*\t"
					"0\t"); // neither has coordinates
		if (tlen > 0)
			++tlen;
		else if (tlen < 0)
			--tlen;
		sprintf_lite(s, "%d\t", tlen);
	}
	else
		sprintf_lite(s,
				"\t"
				"*\t"
				"0\t"
				"0\t");

	// write SEQ and QUAL
	if (r == 0)
	{
		sam_write_sq(s, t->seq, t->l_seq, 0, 0);
		sprintf_lite(s, "\t");
		if (t->qual)
			sam_write_sq(s, t->qual, t->l_seq, 0, 0);
		else
			sprintf_lite(s, "*");
	}
	else
	{
		if ((flag & 0x900) == 0 || (opt_flag & MM_F_SOFTCLIP))
		{
			sam_write_sq(s, t->seq, t->l_seq, r->rev, r->rev);
			sprintf_lite(s, "\t");
			if (t->qual)
				sam_write_sq(s, t->qual, t->l_seq, r->rev, 0);
			else
				sprintf_lite(s, "*");
		}
		else if (flag & 0x100)
		{
			sprintf_lite(s, "*\t*");
		}
		else
		{
			sam_write_sq(s, t->seq + r->qs, r->qe - r->qs, r->rev, r->rev);
			sprintf_lite(s, "\t");
			if (t->qual)
				sam_write_sq(s, t->qual + r->qs, r->qe - r->qs, r->rev, 0);
			else
				sprintf_lite(s, "*");
		}
	}

	// write tags
	if (mm_rg_id[0])
		sprintf_lite(s,
				"\t"
				"RG:Z:%s",
				mm_rg_id);
	if (n_seg > 2)
		sprintf_lite(s,
				"\t"
				"FI:i:%d",
				seg_idx);
	if (r)
	{
		write_tags(s, r);
		if (r->parent == r->id && r->p && n_regs > 1 && regs && r >= regs && r - regs < n_regs)
		{ // supplementary aln may exist
			int i, n_sa = 0; // n_sa: number of SA fields
			for (i = 0; i < n_regs; ++i)
				if (i != r - regs && regs[i].parent == regs[i].id && regs[i].p)
					++n_sa;
			if (n_sa > 0)
			{
				sprintf_lite(s, "\tSA:Z:");
				for (i = 0; i < n_regs; ++i)
				{
					const mm_reg1_t *q = &regs[i];
					int l_M, l_I = 0, l_D = 0, clip5 = 0, clip3 = 0;
					if (r == q || q->parent != q->id || q->p == 0)
						continue;
					if (q->qe - q->qs < q->re - q->rs)
						l_M = q->qe - q->qs, l_D = (q->re - q->rs) - l_M;
					else
						l_M = q->re - q->rs, l_I = (q->qe - q->qs) - l_M;
					clip5 = q->rev? t->l_seq - q->qe : q->qs;
					clip3 = q->rev? q->qs : t->l_seq - q->qe;
					sprintf_lite(s, "%s,%d,%c,", mi->seq[q->rid].name, q->rs+1, "+-"[q->rev]);
					if (clip5)
						sprintf_lite(s, "%dS", clip5);
					if (l_M)
						sprintf_lite(s, "%dM", l_M);
					if (l_I)
						sprintf_lite(s, "%dI", l_I);
					if (l_D)
						sprintf_lite(s, "%dD", l_D);
					if (clip3)
						sprintf_lite(s, "%dS", clip3);
					sprintf_lite(s,
							","
							"%d"
							","
							"%d"
							";",
							q->mapq,
							q->blen - q->mlen + q->p->n_ambi);
				}
			}
		}
		if (r->p && (opt_flag & MM_F_OUT_CS))
			write_cs(s, mi, t, r, !(opt_flag&MM_F_OUT_CS_LONG));
		if (cigar_in_tag)
			write_sam_cigar(s, flag, 1, t->l_seq, r, opt_flag);
	}
	s->s[s->l] = 0; // we always have room for an extra byte (see str_enlarge)
}

void mm_write_sam(
		kstring_t *s,
		const mm_idx_t *mi,
		const mm_bseq1_t *t,
		const mm_reg1_t *r,
		int n_regs,
		const mm_reg1_t *regs)
{
	int i;
	for (i = 0; i < n_regs; ++i)
		if (r == &regs[i])
			break;
	mm_write_sam2(s, mi, t, 0, i, 1, &n_regs, &regs, 0);
}
*/
