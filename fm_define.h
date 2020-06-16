/*************************************************************************
	> File Name: fem_define.h
	> Author: Xu Kai
	> Created Time: Wed 11 Dec 2019 08:20:34 PM CST
 ************************************************************************/
#ifndef __FM_DEFINE__
#define __FM_DEFINE__

#define index_t uint64_t
#define sindex_t uint32_t
#define half_word_t  uint16_t

#define SEED_STEP 4
#define INVALID_SEED 0xffffffffffffffff
#define SINVALID_SEED 0xffffffff
#define INVALID_SFLAG 0xffffffff
#define MAX_CIGAR 256
#define CIGARS_PER_BLOCK 256
#define CIGAR_LEN 16
#define CIGAR_AVG_LEN 16

typedef struct ref_info{
    //in
    sindex_t *hash_count;
    sindex_t *hash_offset;
    sindex_t *hash_positions;
    sindex_t *base_2bits;
    sindex_t *base_3bits;
    half_word_t *pos_lseeds;
    sindex_t *reads_seq_len[2];
    char *reads[2];
    index_t reads_name_length;
    index_t reads_seq_length;
    index_t seed_length;
    index_t errors;
    index_t ref_length;

    index_t cigars_per_read;
    index_t cigar_avg_len;
    index_t reads_num;

    //out
    char *cigar_buffer[2];
    sindex_t *cigar_num[2];
    sindex_t *ref_pos_buffer[2];

    //support 
    index_t *cigar_ind;
    index_t *unmap_ind;
    index_t cur_block;
    index_t cur_read;
    index_t cigar_size;
    index_t cigar_len;
    index_t lock_read;
}ref_info_t;

#endif

	

