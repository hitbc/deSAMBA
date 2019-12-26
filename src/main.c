#include "stdio.h"
#include "string.h"
#include <getopt.h>
#include "lib/desc.h"
#include "lib/utils.h"

int sort_kmers_main(int argc, char *argv[]);
int build_index_main(int argc, char *argv[]);
int classify_main(int argc, char *argv[]);
int simDataTest(int argc, char *argv[]);

//int hash_cmp_main(int argc, char *argv[]);
//int test_ekmer_main(int argc, char *argv[]);
//int test_fastq_main(int argc, char *argv[]);

int main_test_cpp();

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: %s\n", PACKAGE_NAME);
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: %s\n\n", CONTACT);
	fprintf(stderr, "  Usage  : %s <command> \n", PACKAGE_NAME);
	fprintf(stderr, "  Command: \n");
	fprintf(stderr, "		index   :	index the reference sequences\n");
	fprintf(stderr, "		classify:	classify meta-genomic reads\n");
	fprintf(stderr, "		analysis:	analysis results\n");
	fprintf(stderr, "		kmersort:	sort the kmer-result of Jellyfish\n");
	fprintf(stderr, "\n");
	return 1;
}


int main(int argc, char *argv[]){
	int SHOW_TITLE = 1;
	//if(SHOW_TITLE)
	//	fprintf(stderr, "Program: [%s]\tVersion: [%s]\n", PACKAGE_NAME, PACKAGE_VERSION);
		 if (argc <= 1) 						usage();
	else if (strcmp(argv[1],"kmersort") == 0) 	sort_kmers_main(argc - 1, argv + 1);
	else if (strcmp(argv[1],"index")    == 0) 	build_index_main(argc - 1, argv + 1);
	else if (strcmp(argv[1],"classify") == 0) 	classify_main(argc - 1, argv + 1);
	else if (strcmp(argv[1],"analysis") == 0) 	{simDataTest(argc - 1, argv + 1); SHOW_TITLE = 0;}
	else if (strcmp(argv[1],"test_cpp") == 0) 	main_test_cpp();
	else 		 								usage();
	//else if (strcmp(argv[1],"hash_cmp") == 0) 	hash_cmp_main(argc - 1, argv + 1);
	//else if (strcmp(argv[1],"e_kmer") == 0) 	test_ekmer_main(argc - 1, argv + 1);
	//else if (strcmp(argv[1],"test_fastq") == 0) test_fastq_main(argc - 1, argv + 1);

	if(SHOW_TITLE)
		fprintf(stderr, "Normal end program, MAX MEM:[%f]Gbp.\n\n", (float)peakrss()/1024/1024);
	return 0;
}
