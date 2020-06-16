#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<fcntl.h>
#include "htslib/kseq.h"
#include <athread.h>
#include<zlib.h>
#include<mpi.h>


#define MPE



#define index_t uint64_t
#define sindex_t uint32_t
#define half_word_t uint16_t

#define SEED_STEP 4
#define INVALID_SEED 0xffffffffffffffff
#define SINVALID_SEED 0xffffffff
#define MAX_HASH_INDEX (8*10000*10000UL)
#define MB (1024*1024)
#define CIGARS_PER_READ 1536 
#define CIGAR_AVG_LEN 16
#define CIGARS_PER_BLOCK 256
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
    return ptr;
}

extern SLAVE_FUN(find_positions)();

sindex_t seed_length;
sindex_t *hash_count;
sindex_t *hash_offset;
sindex_t *positions;
sindex_t *base_2bits;
sindex_t *base_3bits;
half_word_t *pos_lseeds; 

index_t ref_length;

typedef struct ref_info{
    //in
    sindex_t *hash_count;
    sindex_t *hash_offset;
    sindex_t *hash_positions;
    sindex_t *base_2bits;
    sindex_t *base_3bits;
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

typedef struct hash_table_c{
    sindex_t loc;
    sindex_t read_num;
}hash_table_t;

index_t max_match_size;
sindex_t ref_len[500];
sindex_t ref_name_len[500];
char ref_name[500][25];
sindex_t ref_num;
hash_table_t *loc_table;
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

#define MAXERROR 10000
#define Word unsigned int

#define FIND  2899070901
int myersEditDistance(char* pattern, char *text, int slocation, int readLen, int distThresh)
{
    int band_down = 2 * distThresh;
    int band_length = band_down + 1; 
    int refLen = readLen + band_down; 

  Word Peq[128];
  Word tmp_Peq[128];

  Peq['A'] = Peq['T'] = Peq['G'] = Peq['C'] = (Word)0;

  Word tmp = (Word)1;
  int i;
  for(i =0; i < band_length; i++){
    Peq[pattern[i]]=Peq[pattern[i]] | tmp;
    tmp = tmp << 1;
  }
  //printf("%x\n", Peq['A']);
  //printf("%x\n", Peq['C']);
  //printf("%x\n", Peq['G']);
  //printf("%x\n", Peq['T']);

  tmp_Peq['A'] = Peq['A'];
  tmp_Peq['C'] = Peq['C'];
  tmp_Peq['G'] = Peq['G'];
  tmp_Peq['T'] = Peq['T'];
  int symbol;
  for(symbol = 0; symbol < 128; symbol++) Peq[symbol]=(Word)0;

  Peq['A']=tmp_Peq['A'];
  Peq['C']=tmp_Peq['C'];
  Peq['T']=tmp_Peq['T'];
  Peq['G']=tmp_Peq['G'];

  Word Mask=(Word)1<<(band_length-1);
  Word VP=0;
  Word VN=0;
  Word X=0;
  Word D0=0;
  Word HN=0;
  Word HP=0;
  //int i;
  //for(i = 0; i < refLen; i ++){
  //  cout <<  pattern[i];
  //}
  //cout << endl;

  int err=0;


  Word err_mask=(Word)1;
  int i_bd= band_down;
  int last_high= band_length - readLen + refLen - band_down - 1;
  for(i = 0; i < readLen - 1; i++){
    printf("%x\n", Peq[text[i]]);
    printf("%c\n", text[i]);

    X=Peq[text[i]]|VN;

    D0=((VP+(X&VP))^VP)|X;
    HN=VP&D0;
    HP=VN|~(VP|D0);

    X=D0>>1;
    VN=X&HP;
    VP=HN|~(X|HP);

    printf("VP %x\n", VP); 
    printf("VN %x\n", VN); 

    //printf("%x\n", Peq[text[i]]);
    //printf("%x\n", HP);
    //printf("%x\n", HN);
    //printf("%x\n", VP);
    //printf("%x\n", VN);
    //printf("%x\n", D0);

    //printf("\n");

    //printf("VP %x\n", VP); 
    //printf("VN %x\n", VN); 
    if(!(D0&err_mask)){
      ++err;
      //if((err-last_high)>distThresh) return MAXERROR;
    }

    Peq['A']=Peq['A']>>1;
    Peq['C']=Peq['C']>>1;
    Peq['G']=Peq['G']>>1;
    Peq['T']=Peq['T']>>1;

    ++i_bd;
    Peq[pattern[i_bd]]=Peq[pattern[i_bd]] | Mask;

  //printf("%x\n", Peq['A']);
  //printf("%x\n", Peq['C']);
  //printf("%x\n", Peq['G']);
  //printf("%x\n", Peq['T']);
  }

  X=Peq[text[readLen - 1]]|VN;
  D0=((VP+(X&VP))^VP)|X;
  HN=VP&D0;
  HP=VN|~(VP|D0);
  X=D0>>1;
  VN=X&HP;
  VP=HN|~(X|HP);

    //printf("%x\n", Peq[text[readLen - 1]]);
    //printf("%x\n", HP);
    //printf("%x\n", HN);
    //printf("%x\n", VP);
    //printf("%x\n", VN);
    //printf("%x\n", D0);

    //printf("\n");

 
//    printf("VP %x\n", VP); 
//    printf("VN %x\n", VN); 
  if(!(D0&err_mask)){
    ++err;
    //if((err-last_high)>distThresh) return MAXERROR;
  }

//  cout << "inital " << err << endl;

 // cout << err - last_high << endl;

  int site=refLen-last_high-1;

  int location = -1;
  int error = MAXERROR;
  if((err <= distThresh)&&(err < error)){
    error = err;
    location = site;
  }

    //printf("VP %x\n", VP); 
    //printf("VN %x\n", VN); 
  for(i = 0; i < last_high; i++){
    err=err+((VP>>i)&(Word)1);
    err=err-((VN>>i)&(Word)1);

    if((err<=distThresh)&&(err<error)){
      error = err;
      location = site+i+1;
    }
  }
    printf("myers err %d %d\n", error, (slocation - err) - (location - readLen));
//  cout << "final " << error << endl;
  return error;
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
#define REFDIR "hg384/"  
    //#define REFDIR "chr1/"  
    ref_num = 0;
    fptr = fopen(REFDIR"ref_info","r");
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




    fptr = fopen(REFDIR"hash_count", "rb");
    file_size = get_file_size(fptr);
    hash_count = (sindex_t*)malloc(file_size);
    fread(hash_count, sizeof(sindex_t), file_size/sizeof(sindex_t), fptr);
    fclose(fptr);


    fptr = fopen(REFDIR"hash_offset", "rb");
    file_size = get_file_size(fptr);
    hash_offset = (sindex_t*)malloc(file_size);

    printf("hash offset %lx %lu\n", hash_offset, file_size); 
    fread(hash_offset, sizeof(sindex_t), file_size/sizeof(sindex_t), fptr);
    fclose(fptr);

    fptr = fopen(REFDIR"hash_positions", "rb");
    file_size = get_file_size(fptr) + 128;
    positions = (sindex_t*)malloc(file_size);

    printf("base 2bits %lx %lu\n", positions, file_size); 
    fread(positions, sizeof(sindex_t), file_size/sizeof(sindex_t), fptr);
    fclose(fptr);

    fptr = fopen(REFDIR"pos_lseeds", "rb");
    file_size = get_file_size(fptr) + 128;
    pos_lseeds = (half_word_t*)malloc(file_size);

    printf("base 2bits %lx %lu\n", positions, file_size); 
    fread(pos_lseeds, sizeof(half_word_t), file_size/sizeof(half_word_t), fptr);
    fclose(fptr);



    fptr = fopen(REFDIR"base_2bits", "rb");
    file_size = get_file_size(fptr);
    ref_length = file_size * 16;

    base_2bits = (sindex_t*)malloc(((file_size + 64) / sizeof(sindex_t))*sizeof(sindex_t));
    printf("base 2bits %lx %lu\n", base_2bits, file_size); 
    if(base_2bits == 0){
        printf("wrong !!!!!\n");
    }

    fread(base_2bits, sizeof(sindex_t), file_size/sizeof(sindex_t), fptr);
    fclose(fptr);


    fptr = fopen(REFDIR"base_3bits", "rb");

    file_size = get_file_size(fptr);
    base_3bits = (sindex_t*)malloc(file_size);

    printf("base 3bits %lx %lu\n", base_3bits, file_size); 
    fread(base_3bits, sizeof(sindex_t), file_size/sizeof(sindex_t), fptr);
    fclose(fptr);
    index_t  i;
    for( i = 0; i < 1 << (seed_length*2); i ++){
        hash_offset[i] -= hash_count[i];
    }


}
void reverse_read(char *seq, char *reverse_seq, sindex_t s_len){
    int i = 0;
    for ( i = 0; i < s_len; i ++){
        switch(seq[i]){
            case('A'):
                reverse_seq[s_len - i - 1] = 'T';
                break;
            case('C'):
                reverse_seq[s_len - i - 1] = 'G';
                break;
            case('G'):
                reverse_seq[s_len - i - 1] = 'C';
                break;
            case('T'):
                reverse_seq[s_len - i - 1] = 'A';
                break;
            default:
                reverse_seq[s_len - i - 1] = 'A';
                break;
        }
    }
}
sindex_t get_seed(char *read, sindex_t l_off, sindex_t s_off, sindex_t seed_len){
    sindex_t soff = l_off * ( seed_len + SEED_STEP - 1) + s_off;
    index_t i;
    sindex_t seed = 0;
    for( i = soff; i < seed_len + soff;  i ++){
        switch(read[i]){
            case('A'):
                seed  |= 0;
                break;
            case('C'):
                seed  |= 1U << ((i - soff)*2);
                break;
            case('G'):
                seed  |= 2U << ((i - soff)*2);
                break;
            case('T'):
                seed  |= 3U << ((i - soff)*2);
                break;
            default:
                seed  |= 0;
                break;
        }
    }
    return seed;
        
}
sindex_t get_start_loc(sindex_t loc, sindex_t lind, sindex_t sind, sindex_t seed_len){
    return loc -  (lind * (seed_len + SEED_STEP - 1) + sind); 
}
char int2char(sindex_t code){
    switch(code){
        case(0):
            return 'A';
        case(1):
            return 'C';
        case(2):
            return 'G';
        case(3):
            return 'T';
    }
}
void print_long_seeds(char* seq, sindex_t start_loc, sindex_t seed_len, sindex_t s_len){
    index_t i, k, j;
    sindex_t lseed_len = seed_len + SEED_STEP - 1;
    for(i = 0; i < s_len; i ++){
        printf("%c", seq[i]);
        if((i + 1) % lseed_len == 0){
            printf("\t");
        }
    }
    printf("\n");

    for(i = 0; i < s_len; i ++){
        index_t id_int = ((start_loc + i) / 32) * 2;
        index_t shift = (start_loc + i) %32;

        sindex_t lbit = (base_2bits[id_int] >> shift)&0x1;
        sindex_t hbit = (base_2bits[id_int + 1] >> shift)&0x1;
        sindex_t bs = (hbit << 1) | lbit;
        char base = int2char(bs);
        printf("%c", base);
        if((i + 1) % lseed_len == 0){
            printf("\t");
        }
    }
    printf("\n");


}

void print_ref2char(sindex_t a, sindex_t b){
    index_t i = 0;
    for( i = 0; i < 32; i ++){
        index_t tmp = (a >> i) & 0x1;
        index_t tmp1 = (b >> i) & 0x1;
        tmp = (tmp | (tmp1 << 1));
        switch(tmp){
            case(0):
                printf("%c", 'A');
                break;
            case(1):
                printf("%c", 'C');
                break;
            case(2):
                printf("%c", 'G');
                break;
            case(3):
                printf("%c", 'T');
                break;
        }
        //con >>= 2;
    }
}

void get_ref_seq(char *ref_seq, sindex_t start_loc, sindex_t s_len, sindex_t err){
    start_loc -= err;
    index_t i;
    for(i = 0; i < s_len + 2 * err; i ++){
        index_t id_int = ((start_loc + i) / 32) * 2;
        index_t shift = (start_loc + i) %32;

        sindex_t lbit = (base_2bits[id_int] >> shift)&0x1;
        sindex_t hbit = (base_2bits[id_int + 1] >> shift)&0x1;
        sindex_t bs = (hbit << 1) | lbit;
        char base = int2char(bs);
        ref_seq[i] = base;
    }

}
void print_ref(sindex_t pos, sindex_t s_len, sindex_t seed_len){
    sindex_t lseed_len = seed_len + SEED_STEP - 1;
    char *ref_path = "hg38/hg38.fa";
    kseq_t *kref = kseq_open(ref_path);
    int flag ;
    sindex_t ref_off = 0;
    sindex_t i;
    while(flag = kseq_read(kref)>=0){
       if(pos < ref_off + kref->seq.l){
           sindex_t off = pos - ref_off;
           for(i = 0; i < s_len; i ++){
               printf("%c", toupper(kref->seq.s[off + i]));
               if((i + 1)%lseed_len == 0)
                   printf("\t");
           }
           printf("\n");
           return;
       }
       ref_off += kref->seq.l;
    }
    kseq_destroy(kref);

}

typedef struct lseeds_info_h{
    sindex_t num;
    sindex_t id;
}lseeds_info_t;
void find_pos(char *seq, sindex_t s_len){
    index_t seed_len = 12;
    index_t ls_num = s_len / (seed_len + SEED_STEP - 1); 
    sindex_t err = 4;
    index_t i, k, j;
    sindex_t seeds[256];
    sindex_t seeds_off[256];
    sindex_t seeds_num[256];
    lseeds_info_t lseeds_info[256];
    for( i = 0; i < ls_num; i ++){
        lseeds_info[i].num = 0;
        lseeds_info[i].id = i;
        for(j = 0; j < SEED_STEP; j ++){
            sindex_t seed = get_seed(seq, i, j, seed_len); 
            seeds[i*SEED_STEP + j] = seed;
            seeds_num[i*SEED_STEP + j] = hash_offset[i*SEED_STEP + j + 1] - hash_offset[i*SEED_STEP + j];  
            lseeds_info[i].num += seeds_num[i*SEED_STEP + j];
        }
    }
    for(i = 0; i < ls_num; i ++){
        for(j = i + 1; j < ls_num; j ++){
            if(lseeds_info[i].num > lseeds_info[j].num){
                lseeds_info_t tmp  = lseeds_info[i];
                lseeds_info[i] = lseeds_info[j];
                lseeds_info[j] = tmp;
            }
        }
    }
//#define FIND 18135451  
    print_long_seeds(seq, FIND, seed_len, s_len); 
    print_ref(FIND, s_len, seed_len);
    for(i = 0; i < err + 1; i ++){
        sindex_t lind = lseeds_info[i].id;
        for(j = 0; j < SEED_STEP; j ++){
            sindex_t sind = lind * SEED_STEP + j;
            for(k = hash_offset[seeds[sind]]; k < hash_offset[seeds[sind] + 1]; k ++){
                sindex_t start_loc = get_start_loc(positions[k], lind, j, seed_len);
                sindex_t  f_loc = start_loc >  FIND ? start_loc - FIND : FIND - start_loc; 
                if(f_loc <= err){
                    printf("FIND %u %u %u %u %u\n", i, lind, j, start_loc,  positions[k]);
                    print_ref2char(pos_lseeds[k], pos_lseeds[k] >> (2 * (SEED_STEP - 1))); 
                    print_long_seeds(seq, start_loc, seed_len, s_len); 
                    char ref_seq[256];
                    get_ref_seq(ref_seq, start_loc, s_len, err);
                    myersEditDistance(ref_seq, seq, start_loc, s_len, err);
                }
            }
        }
    }
}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    index_t reads_name_length = 20;
    index_t reads_seq_length = 120;
    index_t reads_num = 10000;
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


    //open reads file
    //char *reads_path = "ERR013135_1.fastq";
    //char *reads_path = "err_test_1.fastq";
    //char *reads_path = "data/data4.fq";
    char *reads_path = "err_test_3.fq";
    int my_id; 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    //char reads_path[40];
    //sprintf(reads_path, "ERRPART/x%04d.fastq", my_id);
    //char result_name[40];
    //sprintf(result_name, "sam_%04d.sam.gz", my_id);
    //kseq_t *kreads = kseq_open(argv[1]);
    // sam_ptr = gzopen(argv[2], "w1");
    printf("reads name: %s\n", reads_path); 
    //printf("result name: %s\n", result_name); 

    kseq_t *kreads = kseq_open(reads_path);

    //printf("reads name: %s\n", argv[1]); 
    //printf("result name: %s\n", argv[2]); 

    int l = 0;
    seed_length = 12;

    load_ref_info();


    index_t cur_block = 0;
    max_match_size = 0;


    char *reads_name = reads_name_buffer[cur_block]; 
    char *reads_seq = reads_seq_buffer[cur_block]; 
    char *reads_qual = reads_qual_buffer[cur_block];

    sindex_t *reads_nl = reads_nl_buffer[cur_block]; 
    sindex_t *reads_sl = reads_sl_buffer[cur_block]; 


    int read_flag  = 0;
    mapped_reads = 0;


    kreads->name.s = reads_name;
    kreads->seq.s = reads_seq;
    kreads->qual.s = reads_qual;

    kreads->name.m = reads_name_length;
    kreads->seq.m = reads_seq_length;
    kreads->qual.m = reads_seq_length;
    while(read_flag >= 0){

        read_flag = kseq_read(kreads); 
        if(read_flag == -1)
            break;
        //printf("%s\n%s\n%s\n", kreads->name.s, kreads->seq.s, kreads->qual.s);

        find_pos(kreads->seq.s, kreads->seq.l);

        char reverse_seq[256]; 

        reverse_read(kreads->seq.s, reverse_seq, kreads->seq.l);
        find_pos(reverse_seq, kreads->seq.l);

        //reads_nl[cur_reads_size] = kreads->name.l; 
        //reads_sl[cur_reads_size] = kreads->seq.l; 
        //reads_nl_buffer[cur_block][cur_reads_size] = kreads->name.l; 
        //reads_sl_buffer[cur_block][cur_reads_size] = kreads->seq.l; 
    }

    printf("max match size %lu\n", max_match_size);
    printf("mapped reads %lu\n", mapped_reads);


    //    lwpf_report_summary(stdout, &conf);



    MPI_Finalize();
    return 0;
}
