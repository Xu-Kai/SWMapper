#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<fcntl.h>
//#include "htslib/kseq.h"
#include "kseq.h"
#include <athread.h>
#include<zlib.h>
#include<mpi.h>
#include<getopt.h>
#include "fm_define.h"
//#define CPE_PROFILE
//#define MPE_PROFILE 
#ifdef MPE_PROFILE
#include "swlu.h"
#endif

#ifdef CPE_PROFILE
#define MPE
#define LWPF_UNITS U(MAPPING)
#include "lwpf2/lwpf2.h"
#endif


#define MAX_HASH_INDEX (8*10000*10000UL)
#define MB (1024*1024)
#define CIGARS_PER_READ 1024 
//#define CIGARS_PER_READ 128 


#define LOCATIONS (CIGARS_PER_BLOCK + 2)

#define PRO_VERSION "1.0"
#define PRO_NAME "XXX"
#define SAM_VERSION "1.4"

const char FWD_MAPPING = '+';
const char REV_MAPPING = '-';

// SAM tags

// Header tags.
const char *HD_TAG = "@HD";
const char *VERSION_TAG = "VN:";
const char *SO_TAG = "SO:";
const char *SEQ_TAG = "@SQ";
const char *SEQNAME_TAG = "SN:";
const char *SEQLEN_TAG = "LN:";
const char *PROG_TAG = "@PG	ID:XXX";
const char *PROGNAME_TAG = "PN:XXX";
const char *CMD_LINE_TAG = "CL:";

// Result tags.
const char CIGAR = 'M';
const char *N_TAG = "N:";
const char *NM_TAG = "NM:i:";
const char *XM_TAG = "XM:i:";
const char *XA_Z_TAG = "XA:Z:";
const char *XA_I_TAG = "XA:i:";
const char *MD_Z_TAG = "MD:Z:";
const char UNAVAILABLE_1 = '*';
const char UNAVAILABLE_2 = '0';
const uint32_t UNAVAILABLE_3 = 255;

// Formatting.
const char SEPARATOR = '\t';
#define SEP  "\t"
const char COMMA = ',';
const char SEMI_COLON = ';';

// SAM FLAGS
#define	READ_UNMAPPED  4
#define	READ_REVERSE_MAPPED  16
#define	READ_PAIR_UNMAPPED1  77
#define	READ_PAIR_UNMAPPED2  141
#define	READ_PAIR_FORWARD_MAPPED1  99
#define	READ_PAIR_REVERSE_MAPPED1  83
#define	READ_PAIR_FORWARD_MAPPED2  147
#define	READ_PAIR_REVERSE_MAPPED2  163
#define	READ_PAIR_FIRST_FORWARD_MAPPED  0
#define	READ_PAIR_FIRST_REVERSE_MAPPED  16
#define	READ_PAIR_FIRST_UNMAPPED  4
#define	READ_PAIR_SECOND_FORWARD_MAPPED  0
#define	READ_PAIR_SECOND_REVERSE_MAPPED  16
#define	READ_PAIR_SECOND_UNMAPPED  4
#define	NOT_PRIMARY_ALIGN  256

#define index_t uint64_t
#define sindex_t uint32_t
#define half_word_t  uint16_t
/*
#define SEED_STEP 4
#define INVALID_SEED 0xffffffffffffffff
#define SINVALID_SEED 0xffffffff
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
*/





index_t total_malloc_size;
gzFile* sam_ptr;
void check_mallc(void *ptr){
    if(ptr == 0){
        printf("malloc wrong\n");
    }
}
void* __wrap_malloc(size_t size){
    void *ptr = (void*)__real_malloc(size);
    if(ptr == 0)
        printf("malloc wrong!!!!!!!");
    total_malloc_size += size;
    return ptr;
}

//extern SLAVE_FUNC(find_positions)();
char prefix[200];
char out_prefix[200];
sindex_t seed_length;
sindex_t *hash_count;
sindex_t *hash_offset;
sindex_t *positions;
sindex_t *base_2bits;
sindex_t *base_3bits;
half_word_t *pos_lseeds;

index_t ref_length;

extern void slave_find_positions(ref_info_t *h);
//typedef struct hash_table_c{
//    sindex_t loc;
//    sindex_t read_num;
//}hash_table_t;

index_t max_match_size;
sindex_t ref_len[500];
sindex_t ref_name_len[500];
char ref_name[500][25];
sindex_t ref_num;
//hash_table_t *loc_table;
index_t cigar_block_size;
index_t cigar_unmap_size;
index_t cigar_block_size_pre;
index_t cigar_unmap_size_pre;
index_t mapped_reads;

KSEQ_INIT(int, read)

    kseq_t *kseq_open(const char *path){
        int fd = open(path, O_RDONLY);
        if(fd == -1){
            printf("open reads file wrong!\n"); 
            MPI_Abort(MPI_COMM_WORLD, 3);
        }
        return kseq_init(fd);
    }

index_t get_file_size(FILE *fp){
    fseek(fp, 0L, SEEK_END);
    index_t sz =  ftell(fp);
    fseek(fp, 0, SEEK_SET);
    return sz;
}
void load_ref_info(){
    index_t file_size;
    FILE *fptr;
//#define REFDIR "hg384/"  
//#define REFDIR "hg384/"  
    //#define REFDIR "chr1/"  
    char file_name[200];
    ref_num = 0;
    //fptr = fopen(REFDIR"ref_info","r");
    sprintf(file_name, "%s%s", prefix, ".ref_info");
    fptr = fopen(file_name,"r");
    char *tmp = NULL;
    int len = 0;
    int rlen = 0;
    ref_len[0] = 0;
    while((rlen = getline(&tmp, &len, fptr)) != -1){
        tmp[rlen - 1] = 0;
        ref_len[ref_num + 1] = (sindex_t)atoi(tmp);

        rlen = getline(&tmp, &len, fptr);
        tmp[rlen - 1] = 0;
        ref_name_len[ref_num] = (sindex_t)atoi(tmp);

        rlen = getline(&tmp, &len, fptr);
        tmp[rlen - 1] = 0;
        memcpy(ref_name[ref_num], tmp, rlen);
        ref_num ++;
    }
    if(tmp)
        free(tmp);
    fclose(fptr);
    for(len = 0; len < ref_num; len ++){
        ref_len[len + 1] += ref_len[len];
        //printf("%u\n", ref_len[len]);
        //printf("%u\n", ref_name_len[len]);
        //printf("%s\n", ref_name[len]);
    }
    printf("ref num %u\n", ref_num);




    sprintf(file_name, "%s%s", prefix, ".hash_count");
    //printf("file name %s\n", file_name);
    //fptr = fopen(REFDIR"hash_count", "rb");

    fptr = fopen(file_name,"rb");
    file_size = get_file_size(fptr);
    hash_count = (sindex_t*)malloc(file_size);
    fread(hash_count, sizeof(sindex_t), file_size/sizeof(sindex_t), fptr);
    fclose(fptr);

    hash_offset = (sindex_t*)malloc(file_size + 128);
    index_t i = 0;
    hash_offset[0] = 0;
    for(i = 0; i< (1<<(seed_length*2)); i ++){
        hash_offset[i + 1] = hash_offset[i] + hash_count[i];
    }
    //printf("file name %s\n", file_name);
    //sprintf(file_name, "%s%s", prefix, ".hash_offset");
    //fptr = fopen(file_name,"r");
    //fptr = fopen(REFDIR"hash_offset", "rb");
    //file_size = get_file_size(fptr);
    //hash_offset = (sindex_t*)malloc(file_size);

    //printf("hash offset %lx %lu\n", hash_offset, file_size); 
    //fread(hash_offset, sizeof(sindex_t), file_size/sizeof(sindex_t), fptr);
    //fclose(fptr);

    sprintf(file_name, "%s%s", prefix, ".positions");
    fptr = fopen(file_name,"rb");
    //fptr = fopen(REFDIR"hash_positions", "rb");
    file_size = get_file_size(fptr);
    positions = (sindex_t*)malloc(file_size + 128);
    //printf("base 2bits %lx %lu\n", positions, file_size); 
    fread(positions, sizeof(sindex_t), file_size/sizeof(sindex_t), fptr);
    fclose(fptr);

    //fptr = fopen(REFDIR"pos_lseeds", "rb");
    sprintf(file_name, "%s%s", prefix, ".pos_lseeds");
    fptr = fopen(file_name,"rb");
    file_size = get_file_size(fptr);
    pos_lseeds = (half_word_t*)malloc(file_size + 128);
    fread(pos_lseeds, sizeof(half_word_t), file_size/sizeof(half_word_t), fptr);
    fclose(fptr);

    sprintf(file_name, "%s%s", prefix, ".base_2bits");
    fptr = fopen(file_name,"rb");
    //fptr = fopen(REFDIR"base_2bits", "rb");
    file_size = get_file_size(fptr);
    ref_length = file_size * 16;

    base_2bits = (sindex_t*)malloc(((file_size + 64) / sizeof(sindex_t))*sizeof(sindex_t));
    //printf("base 2bits %lx %lu\n", base_2bits, file_size); 
    //if(base_2bits == 0){
    //    printf("wrong !!!!!\n");
    //}

    fread(base_2bits, sizeof(sindex_t), file_size/sizeof(sindex_t), fptr);
    fclose(fptr);


    sprintf(file_name, "%s%s", prefix, ".base_3bits");
    fptr = fopen(file_name,"rb");
    //fptr = fopen(REFDIR"base_3bits", "rb");
    file_size = get_file_size(fptr);
    base_3bits = (sindex_t*)malloc(file_size);
    //printf("base 3bits %lx %lu\n", base_3bits, file_size); 
    fread(base_3bits, sizeof(sindex_t), file_size/sizeof(sindex_t), fptr);
    fclose(fptr);
    //index_t  i;
    //for( i = 0; i < 1 << (seed_length*2); i ++){
    //    hash_offset[i] -= hash_count[i];
    //}


}
void sam_header(){
    // Writing HD TAG.
    gzprintf(sam_ptr, "%s%c%s%s%c%s%s\n", HD_TAG, SEPARATOR, VERSION_TAG, 
            SAM_VERSION, SEPARATOR, SO_TAG, "unsorted");

    sindex_t i;
    sindex_t pre_len = 0;
    for( i = 0; i < ref_num; i ++){
        //printf("%s\n", ref_name[i]);
        gzprintf(sam_ptr, "%s%c%s%s%c%s%u\n", SEQ_TAG, SEPARATOR, SEQNAME_TAG, ref_name[i], 
                SEPARATOR, SEQLEN_TAG,  (ref_len[i+1] - pre_len));
        pre_len = ref_len[i + 1];
    }
    //// Writing SQ TAGS.
    //unsigned currLen = 0;
    //unsigned prevLen = 0;
    //for (unsigned i = 0; i < indexer.chromosomeLengths.size(); i++) {
    //	currLen = indexer.chromosomeLengths[i];
    //	resultStream << SEQ_TAG << SEPARATOR << SEQNAME_TAG << indexer.chromosomeNames[i]
    //			<< SEPARATOR << SEQLEN_TAG << (currLen - prevLen) << endl;
    //	prevLen = currLen;
    //}
    // Writing PROG TAG.
    gzprintf(sam_ptr, "%s%c%s%s%c%s%c%s%c%s%c\n", PROG_TAG, SEPARATOR, VERSION_TAG, PRO_VERSION, SEPARATOR,
            PROGNAME_TAG, SEPARATOR, CMD_LINE_TAG, '\"', "CMDLINE", '\"');

}
void reverse_read(char *reverse_seq, char *reverse_qual, char *read_seq, char *read_qual, sindex_t read_len){
    index_t i;
    index_t j = read_len - 1;
    for(i = 0; i < read_len; i ++){
        reverse_qual[j - i] = read_qual[i];
        switch(read_seq[i]){
            case 'A':
                reverse_seq[j-i] = 'T';
                break;
            case 'C':
                reverse_seq[j-i] = 'G';
                break;
            case 'G':
                reverse_seq[j-i] = 'C';
                break;
            case 'T':
                reverse_seq[j-i] = 'A';
                break;
            default:
                reverse_seq[j-i] = read_seq[i];
        }
    }
    reverse_seq[read_len] = 0;
    reverse_qual[read_len] = 0;
}
sindex_t get_ref_name(sindex_t ref_loc, sindex_t pre_ind){

    if( ref_loc >= ref_len[pre_ind] && ref_loc < ref_len[pre_ind + 1])
        return pre_ind;
    if(ref_loc >= ref_len[pre_ind + 1] && ref_loc < ref_len[pre_ind + 2])
        return pre_ind + 1;
    sindex_t left = 0;
    sindex_t right = ref_num;
    sindex_t mid = (0 + ref_num) / 2;
    while( right > left){
        if(ref_loc < ref_len[mid]){
            right = mid - 1; 
        }else if(ref_loc > ref_len[mid]){
            left = mid + 1; 
        }else{
            return mid;
        }
        mid = (left + right) / 2; 
    }
    if(ref_loc < ref_len[mid]){
        return mid - 1;
    }else if (ref_loc >= ref_len[mid]){
        return mid;
    }else{
        return mid;
    }

}
//int probe_loc(sindex_t ref_loc, index_t i){
//    sindex_t ind = ref_loc % 99991; 
//    if(loc_table[ind].loc == ref_loc && loc_table[ind].read_num == i){
//        return 0;
//    }else{
//        loc_table[ind].loc = ref_loc;
//        loc_table[ind].read_num = i;
//        return 1;
//    }
//
//}
#define append_char(curc, c) *curc  = c; curc ++;
#define append_str(curc, s, slen) memcpy(curc, s, slen); curc += slen;
#define append_int(curc, n) curc=append_int_(curc, n);
#define append_const(curc, n) append_str(curc, #n, strlen(#n)); 
char* append_int_(char* curc, int n){
    if(n == 0){
        append_char(curc, '0')
            return curc;
    }
    char c2i[20];
    int  i = 19;
    while(n){
        c2i[i--] = '0' + n%10;
        n /= 10;
    }

    memcpy(curc, c2i + i + 1, 19  - i); 
    curc  += 19 - i;
    return curc;
}
void sam_format_2(char *reads_name, char *reads_seq, char *reads_qual, char *cigar_buffer,
        sindex_t *ref_pos, sindex_t *reads_name_len, sindex_t *reads_seq_len, 
        sindex_t *cigar_num2, index_t reads_size,
        sindex_t read_name_max_len, sindex_t read_seq_max_len){
    index_t i, j, k;

    //uint32_t map_quality = UNAVAILABLE_3;
    //char ref_next = UNAVAILABLE_1; //unavailable
    //unsigned relative_pos = 0;
    //char pos_next = UNAVAILABLE_2;
    //char template_len = UNAVAILABLE_2;

    char reverse_reads[256];
    char reverse_quality[256];

    //unsigned int sam_flags = 0;

    char buffer[1024];
    char *curc = buffer;
    for( i = 0; i < cigar_unmap_size_pre; i ++){
        // sindex_t read_ind = cur_read / 2;
        sindex_t read_ind = cigar_num2[i]/2; 
        //printf("reads num %lu\n", i);
        //sam_flags += READ_UNMAPPED;

        curc = buffer;
        //gzprintf(sam_ptr, "%s%c", reads_name + read_ind * read_name_max_len, SEPARATOR);
        append_str(curc, reads_name + read_ind*read_name_max_len, reads_name_len[read_ind]);

        //gzprintf(sam_ptr, "%u%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%c%s\n",
        //        sam_flags, SEPARATOR, UNAVAILABLE_1, SEPARATOR, UNAVAILABLE_2, 
        //        SEPARATOR, UNAVAILABLE_2, SEPARATOR, UNAVAILABLE_1, SEPARATOR, 
        //        UNAVAILABLE_1, SEPARATOR, UNAVAILABLE_2, SEPARATOR, UNAVAILABLE_2, 
        //        SEPARATOR , reads_seq + read_ind * read_seq_max_len , SEPARATOR, 
        //        reads_qual + read_ind * read_seq_max_len);
        append_char(curc, '\t');
        append_char(curc, '4');
        append_char(curc, '\t');
        append_char(curc, '*');
        append_char(curc, '\t');
        append_char(curc, '0');
        append_char(curc, '\t');
        append_char(curc, '0');
        append_char(curc, '\t');
        append_char(curc, '*');
        append_char(curc, '\t');
        append_char(curc, '*');
        append_char(curc, '\t');
        append_char(curc, '0');
        append_char(curc, '\t');
        append_char(curc, '0');
        append_char(curc, '\t');
        gzwrite(sam_ptr, buffer, curc - buffer);
        gzwrite(sam_ptr, reads_seq + read_ind * read_seq_max_len, reads_seq_len[read_ind]);
        gzwrite(sam_ptr, "\t", 1);
        gzwrite(sam_ptr, reads_qual + read_ind * read_seq_max_len, reads_seq_len[read_ind]);
        gzwrite(sam_ptr, "\n", 1);
    } 
    index_t pre_size = 0;
    for( i = 0; i < cigar_block_size_pre; i ++){
        sindex_t cur_read = ref_pos[i*LOCATIONS]; 
        sindex_t read_ind = cur_read/2; 
        sindex_t cigar_num = ref_pos[i*LOCATIONS + 1]; 
        sindex_t *ref_pos_ptr = ref_pos + i * LOCATIONS + 2;
        //if(cigar_num > max_match_size)
        char *cigar = cigar_buffer + i * CIGARS_PER_BLOCK * CIGAR_AVG_LEN;

        index_t reverse_flag = cur_read & 1; 
        index_t main_flag = cigar_num >> 31;
        cigar_num  = cigar_num & 0xfffffff;

        max_match_size += cigar_num;

        curc = buffer;
        append_str(curc, reads_name + read_ind*read_name_max_len, reads_name_len[read_ind]);
        char *name_ptr = curc;
        if(reverse_flag == 0){
            sindex_t ref_ind = 0;
            for(j = 0; j < cigar_num; j ++){

                //if(main_flag == 0 && j != 0) sam_flags = NOT_PRIMARY_ALIGN;

                sindex_t ref_loc = ref_pos_ptr[j];
                ref_ind = get_ref_name(ref_loc, ref_ind);
                curc = name_ptr;
                //gzprintf(sam_ptr, "%s"SEP"%u"SEP"%s"SEP"%u"SEP,reads_name + read_ind*read_name_max_len, 
                //        sam_flags, ref_name[ref_ind], 
                //        ref_loc - ref_len[ref_ind]);
                //gzprintf(sam_ptr, "%u\t", ref_loc); 

                append_char(curc, '\t');
                if(main_flag == 1 && j == 0) {
                    append_char(curc, '0');
                }
                else{
                    append_const(curc, 256);
                }

                append_char(curc, '\t');
                append_str(curc, ref_name[ref_ind], ref_name_len[ref_ind]);
                append_char(curc, '\t');
                append_int(curc, ref_loc - ref_len[ref_ind]);
                append_char(curc, '\t');


                //gzprintf(sam_ptr, "%u"SEP, map_quality);
                append_const(curc, 255);
                append_char(curc, '\t');


                sindex_t distance  = cigar[0];
                sindex_t cigar_len = cigar[1];
                //gzprintf(sam_ptr, "%s"SEP, cigar + 2);
                append_str(curc, cigar + 2, cigar_len);
                append_char(curc, '\t');
                //gzprintf(sam_ptr, "%c%c%c%c", ref_next, SEPARATOR, pos_next, SEPARATOR);
                append_char(curc, '*');
                append_char(curc, '\t');
                append_char(curc, '0');
                append_char(curc, '\t');

                if(main_flag && j == 0){
                    //gzprintf(sam_ptr, "%c%c%s%c%s%c", template_len, SEPARATOR, 
                    //        reads_seq + read_ind*read_seq_max_len, SEPARATOR,
                    //        reads_qual + read_ind*read_seq_max_len, 
                    //        SEPARATOR);
                    append_char(curc, '0');
                    append_char(curc, '\t');
                    append_str(curc, reads_seq + read_ind * read_seq_max_len, reads_seq_len[read_ind]);
                    append_char(curc, '\t');
                    append_str(curc, reads_qual + read_ind * read_seq_max_len, reads_seq_len[read_ind]);
                    append_char(curc, '\t');
                }
                else{	
                    //gzprintf(sam_ptr, "%c%c%c%c%c%c", template_len, SEPARATOR, '*', SEPARATOR, 
                    //        '*', SEPARATOR);

                    append_char(curc, '0');
                    append_char(curc, '\t');
                    append_char(curc, '*');
                    append_char(curc, '\t');
                    append_char(curc, '*');
                    append_char(curc, '\t');
                }
                //gzprintf(sam_ptr, "%s%u%c%s%u\n", NM_TAG, distance, SEPARATOR, XM_TAG,
                //        distance);
                append_char(curc, 'N');
                append_char(curc, 'M');
                append_char(curc, ':');
                append_char(curc, 'i');
                append_char(curc, ':');
                append_char(curc, '0' + distance);
                append_char(curc, '\t');
                append_char(curc, 'X');
                append_char(curc, 'M');
                append_char(curc, ':');
                append_char(curc, 'i');
                append_char(curc, ':');
                append_char(curc, '0' + distance);
                append_char(curc, '\n');
                gzwrite(sam_ptr, buffer, curc - buffer); 
                cigar = cigar + 2 + cigar_len;
            }
        }else{


            cigar = cigar_buffer + i*CIGARS_PER_BLOCK * CIGAR_AVG_LEN;

            char reverse_seq[256];
            char reverse_qual[256];
            if(main_flag &&  cigar_num > 0){
                reverse_read(reverse_seq, reverse_qual, reads_seq +  read_ind* read_seq_max_len, reads_qual +  read_ind * read_seq_max_len, reads_seq_len[read_ind]);
            }
            sindex_t ref_ind = 0;
            for(j = 0; j < cigar_num; j ++){
                //if(main_flag && j == 0) sam_flags = READ_REVERSE_MAPPED;
                //else sam_flags = READ_REVERSE_MAPPED + NOT_PRIMARY_ALIGN;

                curc = name_ptr;

                sindex_t ref_loc = ref_pos_ptr[j];

                ref_ind = get_ref_name(ref_loc, ref_ind);

                //gzprintf(sam_ptr, "%s%c%u%c%s%c%u%c",reads_name + read_ind*read_name_max_len, SEPARATOR, 
                //        sam_flags, SEPARATOR, ref_name[ref_ind], SEPARATOR, 
                //        ref_loc - ref_len[ref_ind], SEPARATOR);

                append_char(curc, '\t');
                if(main_flag == 1 && j == 0) {
                    append_const(curc, 16);
                }
                else{
                    append_const(curc, 272);
                }
                append_char(curc, '\t');
                append_str(curc, ref_name[ref_ind], ref_name_len[ref_ind]);
                append_char(curc, '\t');
                append_int(curc, ref_loc - ref_len[ref_ind]);
                append_char(curc, '\t');

                //gzprintf(sam_ptr, "%u\t", ref_loc); 
                //gzprintf(sam_ptr, "%u%c", map_quality, SEPARATOR);

                append_const(curc, 255);
                append_char(curc, '\t');


                //cigar[2+cigar_len] = 0;
                sindex_t distance  = cigar[0];
                sindex_t cigar_len = cigar[1];
                //gzprintf(sam_ptr, "%s%c", cigar + 2, SEPARATOR);
                append_str(curc, cigar + 2, cigar_len);
                append_char(curc, '\t');
                //gzprintf(sam_ptr, "%c%c%c%c", ref_next, SEPARATOR, pos_next, SEPARATOR);

                append_char(curc, '*');
                append_char(curc, '\t');
                append_char(curc, '0');
                append_char(curc, '\t');
                if(main_flag && j == 0){
                    //gzprintf(sam_ptr, "%c%c%s%c%s%c", template_len, SEPARATOR, reverse_seq,
                    //       SEPARATOR, reverse_qual, SEPARATOR);
                    append_char(curc, '0');
                    append_char(curc, '\t');
                    append_str(curc, reverse_seq, reads_seq_len[read_ind]);
                    append_char(curc, '\t');
                    append_str(curc, reverse_qual, reads_seq_len[read_ind]);
                    append_char(curc, '\t');
                }
                else{	
                    //gzprintf(sam_ptr, "%c%c%c%c%c%c", template_len, SEPARATOR, '*', SEPARATOR, 
                    //        '*', SEPARATOR);
                    append_char(curc, '0');
                    append_char(curc, '\t');
                    append_char(curc, '*');
                    append_char(curc, '\t');
                    append_char(curc, '*');
                    append_char(curc, '\t');
                }
                //                    gzprintf(sam_ptr, "%s%u%c%s%u\n", NM_TAG, distance, SEPARATOR, XM_TAG,
                //                            distance);
                append_char(curc, 'N');
                append_char(curc, 'M');
                append_char(curc, ':');
                append_char(curc, 'i');
                append_char(curc, ':');
                append_char(curc, '0' + distance);
                append_char(curc, '\t');
                append_char(curc, 'X');
                append_char(curc, 'M');
                append_char(curc, ':');
                append_char(curc, 'i');
                append_char(curc, ':');
                append_char(curc, '0' + distance);
                append_char(curc, '\n');
                gzwrite(sam_ptr, buffer, curc - buffer); 
                cigar = cigar + 2 + cigar_len;
            }
        }
    }
}
void printf_help(){
    printf("-l\tmax length (less than 200)\n");
    printf("-o\toutput file prefix\n");
    printf("-e\terror size\n");
    printf("-p\tindex prefix\n");
    printf("-s\t[opt, default:12] seed length\n");
    printf("-i\t[opt, default:12] single data input\n");
    exit(0);
}
int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

#ifdef MPE_PROFILE
    GPTLinitialize();
    fem_swlu_prof_init();
    fem_swlu_prof_init();
    fem_swlu_prof_start();
#endif
    athread_init();
#ifdef CPE_PROFILE
    perf_config_t conf;
    conf.pcrc = PCRC_ALL;
    conf.pcr0 = PC0_CYCLE;
    conf.pcr1 = PC1_CYCLE;
    conf.pcr2 = PC2_CNT_GLD;
    lwpf_init(&conf);
#endif


    char *optstring  = "h::l::s::p:i::o::e:f::";
    //char *optstring  = "h:::p:";
    int opt;

    char reads_path[512];
    char result_name[512];
    char reads_file_name[512];
    int file_flag = 0;
    index_t reads_name_length = 20;
    index_t reads_seq_length = 120;
    index_t reads_num = 10000;
    index_t errors = 4;
    seed_length = 12;
    while((opt = getopt(argc,argv,optstring)) != -1){
        switch(opt){
            case ('l'):
                reads_seq_length = (index_t)atoi(optarg);
                ///memcpy(ref_path, optarg, strlen(optarg));
                break;
            case('p'):
                memcpy(prefix, optarg, strlen(optarg));
                break;
            case('s'):
                seed_length = (index_t)atoi(optarg);
                break;
            case('h'):
                printf_help();
                break;
            case('o'):
                memcpy(out_prefix, optarg, strlen(optarg));
                break;
            case('e'):
                errors = (index_t)(atoi(optarg));
                break;
            case('i'):
                memcpy(reads_path, optarg, strlen(optarg));
                break;
            case('f'):
                file_flag = 1;
                memcpy(reads_file_name, optarg, strlen(optarg));
                break;

        }
    }


    total_malloc_size = 0;
    //sam_ptr = gzopen("data4.sam.gz", "w3");
    //index_t reads_num = 2000000;
    index_t cigars_per_read = CIGARS_PER_READ;
    index_t cigar_avg_len = CIGAR_AVG_LEN;

    //index_t reads_name_max_len = reads_name_length;
    //index_t reads_seq_max_len = reads_seq_length;

    char *reads_name_buffer[2];
    char *reads_seq_buffer[2];
    char *reads_qual_buffer[2];

    sindex_t *reads_nl_buffer[2];
    sindex_t *reads_sl_buffer[2];

    reads_name_buffer[0] = (char*)malloc(reads_num * reads_name_length);
    reads_name_buffer[1] = (char*)malloc(reads_num * reads_name_length);
    reads_seq_buffer[0] = (char*)malloc(reads_num * reads_seq_length);
    reads_seq_buffer[1] = (char*)malloc(reads_num * reads_seq_length);
    reads_qual_buffer[0] = (char*)malloc(reads_num * reads_seq_length);
    reads_qual_buffer[1] = (char*)malloc(reads_num * reads_seq_length);

    reads_nl_buffer[0] = (sindex_t*)malloc(reads_num * sizeof(sindex_t));
    reads_nl_buffer[1] = (sindex_t*)malloc(reads_num * sizeof(sindex_t));
    reads_sl_buffer[0] = (sindex_t*)malloc(reads_num * sizeof(sindex_t));
    reads_sl_buffer[1] = (sindex_t*)malloc(reads_num * sizeof(sindex_t));

    char *cigar_buffer[2];
    int *tmp[2];
#define CIGAR
#ifdef CIGAR
    tmp[0] = (int*)malloc((2*reads_num * cigars_per_read * cigar_avg_len + 12));
    tmp[1] = (int*)malloc((2*reads_num * cigars_per_read * cigar_avg_len + 12));
#endif
    //cigar_buffer[0] = (char*)malloc(2*reads_num * cigars_per_read * cigar_avg_len);
    //cigar_buffer[1] = (char*)malloc(2*reads_num * cigars_per_read * cigar_avg_len);

    cigar_buffer[0] = (char*)tmp[0];
    cigar_buffer[1] = (char*)tmp[1];

    sindex_t *cigar_num[2];

    cigar_num[0] = (sindex_t*)malloc(2*reads_num * sizeof(sindex_t));
    cigar_num[1] = (sindex_t*)malloc(2*reads_num * sizeof(sindex_t));

    memset(cigar_num[0], 0, 2 * reads_num * sizeof(sindex_t));
    memset(cigar_num[1], 0, 2 * reads_num * sizeof(sindex_t));

#ifdef CIGAR
    memset(cigar_buffer[0], 0, 2 * reads_num * cigars_per_read * cigar_avg_len + 12 );
    memset(cigar_buffer[1], 0, 2 * reads_num * cigars_per_read * cigar_avg_len + 12);
#endif
    sindex_t *ref_pos_buffer[2];
#ifdef CIGAR
    ref_pos_buffer[0] = (sindex_t*)malloc(2*reads_num * cigars_per_read * sizeof(sindex_t));
    ref_pos_buffer[1] = (sindex_t*)malloc(2*reads_num * cigars_per_read * sizeof(sindex_t));
#endif

#ifdef CIGAR
    //loc_table = (hash_table_t*)malloc(100000 * sizeof(loc_table));
    //memset(loc_table, 0, 100000*sizeof(loc_table));
#endif



    //open reads file
    int my_id; 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    //printf("reads name: %s\n", reads_path); 
    //printf("result name: %s\n", result_name); 




    //kseq_t *kreads = kseq_open(argv[1]);
    //sam_ptr = gzopen(argv[2], "w1");

    sprintf(result_name, "%s.sam.gz", out_prefix);


    if(file_flag){
      FILE *lf = fopen(reads_file_name, "r");
      int i;
      for(i = 0; i <= my_id; i ++){
        fscanf(lf, "%s\n", reads_path);
      //  printf("reads name 0: %s\n", reads_path);
      }
      sprintf(result_name, "%s_%04d.sam.gz", out_prefix, my_id);
      //printf("reads name: %s\n", reads_path);
    }
    kseq_t *kreads = kseq_open(reads_path);
    sam_ptr = gzopen(result_name, "w1");

    if(sam_ptr == NULL){
        printf("open results file wrong!\n");
    }


    int l = 0;
    //seed_length = 12;

#ifdef MPE_PROFILE
    GPTLstart("load ref");
#endif
    load_ref_info();
#ifdef MPE_PROFILE
    GPTLstop("load ref");
#endif
#ifdef MPE_PROFILE
    GPTLstart("Total time");
#endif
    sam_header(); 

    //init slave parameter
    ref_info_t ref_info_h;
    ref_info_h.hash_count = hash_count;
    ref_info_h.hash_offset = hash_offset;
    ref_info_h.hash_positions = positions;
    ref_info_h.base_2bits = base_2bits;
    ref_info_h.base_3bits = base_3bits;
    ref_info_h.pos_lseeds = pos_lseeds;
    ref_info_h.reads_seq_len[0] = reads_sl_buffer[0];
    ref_info_h.reads_seq_len[1] = reads_sl_buffer[1];
    ref_info_h.reads[0] = reads_seq_buffer[0];
    ref_info_h.reads[1] = reads_seq_buffer[1];
    ref_info_h.reads_name_length = reads_name_length;
    ref_info_h.reads_seq_length = reads_seq_length;
    ref_info_h.seed_length = seed_length;
    ref_info_h.errors = errors;
    ref_info_h.ref_length = ref_length;
    ref_info_h.cigars_per_read = cigars_per_read;
    ref_info_h.cigar_avg_len = cigar_avg_len;
    ref_info_h.cigar_buffer[0] = cigar_buffer[0];
    ref_info_h.cigar_buffer[1] = cigar_buffer[1];
    ref_info_h.cigar_num[0] = cigar_num[0];
    ref_info_h.cigar_num[1] = cigar_num[1];
    ref_info_h.ref_pos_buffer[0] = ref_pos_buffer[0];
    ref_info_h.ref_pos_buffer[1] = ref_pos_buffer[1];
    //support
    ref_info_h.cur_block = 0;
    printf("error is %lu \n", ref_info_h.errors);

    index_t cigar_ind = 0;
    index_t unmap_ind = 0;
    ref_info_h.cigar_ind = &cigar_ind;
    ref_info_h.unmap_ind = &unmap_ind;


    max_match_size = 0;


    index_t cur_block = 0; 
    index_t pre_block = 0;
    index_t pre_reads_size = 0;
    index_t cur_reads_size = 0;


    char *reads_name = reads_name_buffer[cur_block]; 
    char *reads_seq = reads_seq_buffer[cur_block]; 
    char *reads_qual = reads_qual_buffer[cur_block];

    sindex_t *reads_nl = reads_nl_buffer[cur_block]; 
    sindex_t *reads_sl = reads_sl_buffer[cur_block]; 

    char *cigar = cigar_buffer[cur_block]; 

    int read_flag  = 0;
    mapped_reads = 0;


    kreads->name.s = reads_name;
    kreads->seq.s = reads_seq;
    kreads->qual.s = reads_qual;

    kreads->name.m = reads_name_length;
    kreads->seq.m = reads_seq_length;
    kreads->qual.m = reads_seq_length;
#ifdef MPE_PROFILE
    GPTLstart("load reads");
#endif
    while(read_flag >= 0 && cur_reads_size < reads_num){

        read_flag = kseq_read(kreads); 
        if(read_flag == -1)
            break;
        //printf("%s\n%s\n%s\n", kreads->name.s, kreads->seq.s, kreads->qual.s);
        kreads->name.s = kreads->name.s + reads_name_length;
        kreads->seq.s = kreads->seq.s + reads_seq_length;
        kreads->qual.s = kreads->qual.s + reads_seq_length;

        reads_nl[cur_reads_size] = kreads->name.l; 
        reads_sl[cur_reads_size] = kreads->seq.l; 
        //reads_nl_buffer[cur_block][cur_reads_size] = kreads->name.l; 
        //reads_sl_buffer[cur_block][cur_reads_size] = kreads->seq.l; 
        cur_reads_size ++;
    }
#ifdef MPE_PROFILE
    GPTLstop("load reads");
#endif


    ref_info_h.cur_block = cur_block; 
    ref_info_h.reads_num = cur_reads_size; 
    ref_info_h.lock_read = 0;
#ifdef MPE_PROFILE
    GPTLstart("athread compute");
#endif
    athread_spawn(find_positions, &ref_info_h);
#define ASYN
#ifndef ASYN
    athread_join();
#ifdef MPE_PROFILE
    GPTLstop("athread compute");
#endif
    cigar_block_size = *ref_info_h.cigar_ind;
    cigar_unmap_size = *ref_info_h.unmap_ind;
    //printf("cur_block_size == %lu\n", cigar_block_size);
#endif

    int loop_count = 0;
    while(read_flag >= 0){
        //printf("while loop\n");
        printf("my_id: %d loop count %d\n", my_id, loop_count++);
        pre_block = cur_block;
        cur_block = (cur_block + 1) % 2;
        pre_reads_size = cur_reads_size;
        cur_reads_size = 0;
#ifndef ASYN
        cigar_block_size_pre =   cigar_block_size;
        cigar_unmap_size_pre =  cigar_unmap_size;
#endif
        reads_name = reads_name_buffer[cur_block]; 
        reads_seq = reads_seq_buffer[cur_block]; 
        reads_qual = reads_qual_buffer[cur_block];

        reads_nl = reads_nl_buffer[cur_block]; 
        reads_sl = reads_sl_buffer[cur_block]; 

        //cigar = cigar_buffer[cur_block]; 
        //read next block reads

        kreads->name.s = reads_name;
        kreads->seq.s = reads_seq;
        kreads->qual.s = reads_qual;

        kreads->name.m = reads_name_length;
        kreads->seq.m = reads_seq_length;
        kreads->qual.m = reads_seq_length;


#ifdef MPE_PROFILE
        GPTLstart("load reads");
#endif
        while(read_flag >= 0 && cur_reads_size < reads_num){
            read_flag = kseq_read(kreads);
            if(read_flag == -1)
                break;
            kreads->name.s = kreads->name.s + reads_name_length;
            kreads->seq.s = kreads->seq.s + reads_seq_length;
            kreads->qual.s = kreads->qual.s + reads_seq_length;

            reads_nl[cur_reads_size] = kreads->name.l; 
            reads_sl[cur_reads_size] = kreads->seq.l; 

            cur_reads_size ++;
        }
#ifdef MPE_PROFILE
        GPTLstop("load reads");
#endif
        //printf("wait job join\n");
#ifdef ASYN
        athread_join();
#ifdef MPE_PROFILE
        GPTLstop("athread compute");
#endif
        cigar_block_size = *ref_info_h.cigar_ind;
        cigar_unmap_size = *ref_info_h.unmap_ind;
        cigar_block_size_pre =   cigar_block_size;
        cigar_unmap_size_pre =  cigar_unmap_size;
        printf("my_id: %d cur_block_size == %lu\n", my_id, cigar_block_size);
#endif

        ref_info_h.cur_block = cur_block; 
        ref_info_h.reads_num = cur_reads_size; 
        ref_info_h.lock_read = 0;
        *ref_info_h.cigar_ind = 0;
        *ref_info_h.unmap_ind = 0;

#ifdef MPE_PROFILE
        GPTLstart("athread compute");
#endif
        athread_spawn(find_positions, &ref_info_h);
#ifndef ASYN
        athread_join();
#ifdef MPE_PROFILE
        GPTLstop("athread compute");
#endif
        cigar_block_size = *ref_info_h.cigar_ind;
        cigar_unmap_size = *ref_info_h.unmap_ind;
        printf("my_id: %d cur_block_size == %lu\n", my_id, cigar_block_size);
#endif

#ifdef MPE_PROFILE
        GPTLstart("write results");
#endif
#define IOPUT
#ifdef IOPUT
        //write SAM to File
        sam_format_2(reads_name_buffer[pre_block], reads_seq_buffer[pre_block], reads_qual_buffer[pre_block], cigar_buffer[pre_block],
                ref_pos_buffer[pre_block], reads_nl_buffer[pre_block], reads_sl_buffer[pre_block], 
                cigar_num[pre_block], pre_reads_size,
                reads_name_length, reads_seq_length);
#endif 
#ifdef MPE_PROFILE
        GPTLstop("write results");
#endif

    }

    pre_block = cur_block;
    pre_reads_size = cur_reads_size;

#ifdef ASYN
    athread_join();
#ifdef MPE_PROFILE
    GPTLstop("athread compute");
#endif
    cigar_block_size = *ref_info_h.cigar_ind;
    cigar_unmap_size = *ref_info_h.cigar_ind;
    printf("my_id: %d cur_block_size == %lu\n", my_id, cigar_block_size);
#endif
    cigar_block_size_pre =   cigar_block_size;
    cigar_unmap_size_pre =  cigar_unmap_size;
#ifdef MPE_PROFILE
    GPTLstart("write results");
#endif
#ifdef IOPUT
    sam_format_2(reads_name_buffer[pre_block], reads_seq_buffer[pre_block], reads_qual_buffer[pre_block], cigar_buffer[pre_block],
            ref_pos_buffer[pre_block], reads_nl_buffer[pre_block], reads_sl_buffer[pre_block], 
            cigar_num[pre_block], pre_reads_size,
            reads_name_length, reads_seq_length);
#endif
#ifdef MPE_PROFILE
    GPTLstop("write results");
    gzclose(sam_ptr);




    GPTLstop("Total time");
#endif
    printf("max match size %lu\n", max_match_size);
    //printf("mapped reads %lu\n", mapped_reads);
    double msize = (double)total_malloc_size / 1024 / 1024 / 1024;
    printf("memory usage :  %.2f GB\n", msize);


#ifdef MPE_PROFILE
    fem_swlu_prof_stop();
    fem_swlu_prof_print();
#endif
#ifdef CPE_PROFILE
       lwpf_report_summary(stdout, &conf);
#endif
#ifdef MPE_PROFILE
    GPTLpr(my_id);
    GPTLpr_summary(MPI_COMM_WORLD);
#endif



    MPI_Finalize();
    return 0;
}
