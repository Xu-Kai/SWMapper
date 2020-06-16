//#include<zlib.h>
#include<stdio.h>
#include<stdint.h>
#include<fcntl.h>
//#include "htslib/kseq.h"
#include "kseq.h"
#include<assert.h>
#include <athread.h>
#include<getopt.h>

#ifdef MPE_PROFILE
#include "swlu.h"
#endif

#define index_t uint64_t
#define sindex_t uint32_t
#define half_word_t uint16_t


#define SEED_STEP 4
#define INVALID_SEED 0xffffffffffffffff
#define SINVALID_SEED 0xffffffff
#define MAX_HASH_INDEX (8*10000*10000UL)
#define MB (1024*1024)

#define HOST_BUILD

index_t total_malloc_size;


sindex_t position;
index_t seed_length;
index_t cur_seed_length;
index_t seed;
index_t chromo_position;
index_t chromo_length;
index_t valid_seed_size;

sindex_t *seed_count;
sindex_t *hash_index;
//sindex_t *hash_index_offset;
sindex_t *positions;
sindex_t *positions_hash;
sindex_t *base_2bits;
sindex_t *seeds_2bits;
sindex_t *base_3bits;
half_word_t *pos_lseeds;

sindex_t enc_2bits;
sindex_t enc_hbits;
sindex_t enc_lbits;
sindex_t enc_seeds_2bits;
sindex_t enc_3bits;
index_t chromo_cnt;
char ref_name[500][25];
sindex_t ref_name_len[500];
sindex_t ref_length[500];

char ref_path[200];
char prefix[200]; 
char out_prefix[200];


void* __wrap_malloc(size_t size){
    void *ptr = (void*)__real_malloc(size);
    if(ptr == 0)
        printf("malloc wrong!!!!!!!\n");
    total_malloc_size += size;
    return ptr;
}

//#define OUTDIR "hg38_cpe/"
void write_base_2bits_to_file(index_t size){
    FILE *fptr; 
    printf("write base_2bits\n");
    sprintf(out_prefix, "%s%s", prefix, ".base_2bits");
    //fptr = fopen(OUTDIR"base_2bits", "wb");
    fptr = fopen(out_prefix, "wb");
    fwrite(base_2bits, sizeof(sindex_t), ((size + 31)/32) * 2, fptr);
    fclose(fptr);


}
void write_positions_to_file(index_t size){
    FILE *fptr; 
    printf("write positions\n");
    sprintf(out_prefix, "%s%s", prefix, ".positions");
    //fptr = fopen(OUTDIR"hash_positions", "wb");
    fptr = fopen(out_prefix, "wb");
    fwrite(positions, sizeof(sindex_t), valid_seed_size, fptr);
    fclose(fptr);


}

void write_to_file(){
    FILE *fptr; 
    printf("write hash_count\n");
    sprintf(out_prefix, "%s%s", prefix, ".hash_count");
    //fptr = fopen(OUTDIR"hash_count", "wb");
    fptr = fopen(out_prefix, "wb");
    fwrite(seed_count, sizeof(sindex_t), (1 << (seed_length * 2)), fptr);
    fclose(fptr);

    printf("write hash_offset\n");
    sprintf(out_prefix, "%s%s", prefix, ".hash_offset");
    //fptr = fopen(OUTDIR"hash_offset", "wb");
    fptr = fopen(out_prefix, "wb");
    fwrite(hash_index, sizeof(sindex_t), (1 << (seed_length * 2)) + 1, fptr);
    fclose(fptr);

    //printf("write hash_positions\n");
    //fptr = fopen(OUTDIR"hash_positions", "wb");
    //fwrite(positions, sizeof(sindex_t), valid_seed_size, fptr);
    //fclose(fptr);

    printf("write pos_lseeds\n");
    sprintf(out_prefix, "%s%s", prefix, ".pos_lseeds");
    //fptr = fopen(OUTDIR"pos_lseeds", "wb");
    fptr = fopen(out_prefix, "wb");
    fwrite(pos_lseeds, sizeof(half_word_t), valid_seed_size, fptr);
    fclose(fptr);


    //  printf("write base_2bits\n");
    //  fptr = fopen(OUTDIR"base_2bits", "wb");
    //  fwrite(base_2bits, sizeof(sindex_t), (position + 15)/16, fptr);
    //  fclose(fptr);


    printf("write base_3bits\n");
    sprintf(out_prefix, "%s%s", prefix, ".base_3bits");
    //fptr = fopen(OUTDIR"base_3bits", "wb");
    fptr = fopen(out_prefix, "wb");
    fwrite(base_3bits, sizeof(sindex_t), (position + 31)/32, fptr);
    fclose(fptr);

}


KSEQ_INIT(int,read)
    kseq_t *kseq_open(const char *path){
        int fd = open(path, O_RDONLY);
        return kseq_init(fd);
    }
char get_base(char *kseq){
    if(chromo_position <  chromo_length){
        // if((position >= 4364363UL) && position < (4364363UL + 118)){
        //   printf("%c", toupper(kseq[chromo_position]));

        //}
        return  toupper(kseq[chromo_position ++]);
    }else{
        return 0;
    }
}
index_t get_seed(char *kseq){
    index_t mask = (1UL << (seed_length * 2)) - 1; 
    char base; 
    while(cur_seed_length < seed_length){
        base = get_base(kseq);
        if(base == 0)
            break;
        switch(base){
            case 'A':
                seed = (seed >> 2) | (0UL << ((seed_length - 1) * 2));
                cur_seed_length ++;

                enc_lbits |= 0UL << (position %32);
                enc_hbits |= 0UL << (position %32);

                enc_seeds_2bits |=   0UL << ((position%16)*2);

                break;
            case 'C':
                seed = (seed >> 2) | (1UL << ((seed_length - 1) * 2));
                cur_seed_length ++;


                enc_lbits |= 1UL << (position%32);
                enc_hbits |= 0UL << (position%32);

                enc_seeds_2bits |=   1UL << ((position%16)*2);
                break;
            case 'G':
                seed = (seed >> 2) | (2UL << ((seed_length - 1) * 2));
                cur_seed_length ++;

                enc_lbits |= 0UL << (position%32);
                enc_hbits |= 1UL << (position%32);

                enc_seeds_2bits |=   2UL << ((position%16)*2);
                break;
            case 'T':
                seed = (seed >> 2) | (3UL << ((seed_length - 1) * 2));
                cur_seed_length ++;

                enc_lbits |= 1UL << (position%32);
                enc_hbits |= 1UL << (position%32);

                enc_seeds_2bits |=   3UL << ((position %16)*2);
                break;
            default:

                enc_lbits |= 0UL << (position %32);
                enc_hbits |= 0UL << (position %32);

                enc_3bits |=   1UL << (position %32);

                enc_seeds_2bits |=   0UL << ((position %16)*2);

                cur_seed_length = 0;
                seed = 0;
                break;
        }
        //printf("base %llu %c %x\n", position, base, enc_2bits);
        if((position + 1) % 32 == 0){
            base_2bits[2*(position/32) ] = enc_lbits;
            base_2bits[2*(position/32) + 1] = enc_hbits;
            enc_lbits = 0;
            enc_hbits = 0;
        }
        if((position + 1) % 16 == 0){
            //base_2bits[position /16 - 1] = enc_2bits;
            //enc_2bits = 0;
            seeds_2bits[position /16] = enc_seeds_2bits;
            enc_seeds_2bits = 0;
        }
        if((position + 1)% 32 == 0){
            base_3bits[position/32] = enc_3bits;
            enc_3bits = 0;
        }

        position ++;
    }
    if(cur_seed_length == seed_length){
        cur_seed_length --;
    }

    if(base ==  0){
        return INVALID_SEED; 
    }

    return seed;
}
index_t get_valid_seed(char *kseq){

    index_t valid_seed = get_seed(kseq); 
    return valid_seed;

}

void chromo_hash_index(char *kseq){
    while(1){
        index_t valid_seed = get_valid_seed(kseq);
        if(valid_seed == INVALID_SEED){
            //printf("pos %llu %x\n", position, valid_seed);
            break;
            //continue;
        }
        if((position - seed_length)%SEED_STEP == 0){
            seed_count[valid_seed] ++;
            valid_seed_size ++;
        }

    }


}

void construct_hash_offset(){
    hash_index[0] = 0; 
    index_t i = 0;
    for(i = 0; i < 1 << (seed_length * 2); i ++){
        hash_index[i+1] = hash_index[i] + seed_count[i];
    }

}
half_word_t lseed_enc;
sindex_t get_seed_from_enc_bits(index_t pos){
    //sindex_t valid_seed_mask = 0x24924924; 
    sindex_t valid_seed_mask = (1UL << (seed_length * 2)) - 1; 
    sindex_t valid_flag_mask = (1UL << (seed_length)) - 1;

    //sindex_t valid_seed = ((base_2bits[pos/16] >> ((pos % 16) * 2)) | (base_2bits[pos/16 + 1] << ((16 - pos%16) * 2))) & valid_seed_mask;
    sindex_t valid_seed = ((seeds_2bits[pos/16] >> ((pos % 16) * 2)) | (seeds_2bits[pos/16 + 1] << ((16 - pos%16) * 2))) & valid_seed_mask;
    assert(SEED_STEP <= 5);
    assert(SEED_STEP >= 1);

    index_t s_pos = 0;  
    sindex_t res_len = SEED_STEP - 1; 
    sindex_t mask = (1U<<res_len) - 1;
    if(pos >= res_len)
        s_pos = pos - res_len;
    sindex_t shift = s_pos % 32;
    s_pos /= 32;

    sindex_t flseed = (base_2bits[s_pos * 2]  >> shift) |  (base_2bits[(s_pos + 1) * 2] << (32 - shift));
    sindex_t fhseed = (base_2bits[s_pos * 2 + 1]  >> shift) |  (base_2bits[(s_pos + 1) * 2 + 1]  << (32 - shift));
    flseed = (flseed & mask); 
    fhseed = (fhseed & mask);

    s_pos = pos + seed_length;
    shift = s_pos % 32;
    s_pos /= 32;

    sindex_t blseed = (base_2bits[s_pos * 2]  >> shift) |  (base_2bits[(s_pos + 1) * 2]  << (32 - shift));
    sindex_t bhseed = (base_2bits[s_pos * 2 + 1]  >> shift) |  (base_2bits[(s_pos + 1) * 2 + 1]  << (32 - shift));
    blseed = (blseed & mask); 
    bhseed = (bhseed & mask);

    flseed = flseed | (blseed << res_len);
    fhseed = fhseed | (bhseed << res_len);

    lseed_enc = flseed | (fhseed << (res_len * 2));

    //  printf("cons1 : %x %x\n", base_2bits[pos/16], base_2bits[pos/16 + 1]);
    //  printf("cons2 : %x %x\n", base_2bits[pos/16] >> ((pos % 16) * 2), base_2bits[pos/16 + 1] << ((16 - pos%16) * 2));
    //  printf("cons : %x \n", valid_seed);

    sindex_t valid_flag = ((base_3bits[pos/32] >> (pos % 32)) | (base_3bits[pos/32 + 1] << (32 - pos%32))) & valid_flag_mask;

    //printf("cons3 : %x %x %x %x\n", base_3bits[pos/32], valid_flag, valid_seed_mask, valid_flag_mask);
    if(valid_flag != 0){
        valid_seed = SINVALID_SEED;
    }
    //printf("pos ind %llu %x\n", pos, valid_seed);
    return valid_seed;


}
void store_positions(){
    index_t  i = 0;
    sindex_t valid_seed = 0;
    index_t seed_cnt = 0;

    hash_index[0] = 0; 
    for(i = 0; i < 1 << (seed_length * 2); i ++){
        hash_index[i+1] = hash_index[i] + seed_count[i];
    }


    for(i = 0; i < position - seed_length + 1; i += SEED_STEP){
        valid_seed = get_seed_from_enc_bits(i);
        if(valid_seed != SINVALID_SEED){
            positions[hash_index[valid_seed] ] = i;
            pos_lseeds[hash_index[valid_seed]] = lseed_enc;
            hash_index[valid_seed] ++;
            //printf("valid seed  : %x\n", valid_seed);
            seed_cnt ++;
        }
    }


    printf("valid seed size : %llu\n", seed_cnt);
}


void construct_hash_index(){

    seed_count = (sindex_t*)malloc((1 << (seed_length * 2))  * sizeof(sindex_t));
    memset(seed_count, 0, (1 << (seed_length * 2)) * sizeof(sindex_t));
    hash_index = (sindex_t*)malloc(((1 << (seed_length * 2)) +  1)  * sizeof(sindex_t)); 

    //positions = (sindex_t*)malloc(MAX_HASH_INDEX * sizeof(sindex_t)); 
    base_2bits = (sindex_t*)malloc(800*MB);
    seeds_2bits = (sindex_t*)malloc(800*MB);
    base_3bits = (sindex_t*)malloc(400*MB);

    //pos_lseeds = (half_word_t*) malloc(MAX_HASH_INDEX * sizeof(half_word_t));

    //enc_2bits = 0;
    enc_lbits = 0;
    enc_hbits = 0;
    enc_seeds_2bits = 0;
    enc_3bits = 0;
    memset(seeds_2bits, 0, 800*MB);
    memset(base_2bits, 0, 800*MB);
    //positions_hash = (sindex_t*)malloc(MAX_HASH_INDEX * sizeof(sindex_t)); 

    valid_seed_size = 0; 
    position = 0;
    //char *path = "hg38_test.fa";
    //char *path = "hg38/hg38.fa";
    //char *path = "chr1/chr1.fa";
    printf("ref path %s\n", ref_path);
    kseq_t *kref = kseq_open(ref_path);
    int l ;
    cur_seed_length = 0;
    seed = 0;
    chromo_cnt = 0;
    //FILE* ri_ptr = fopen(OUTDIR"ref_info", "w");
    sprintf(out_prefix, "%s%s", prefix, ".ref_info");
    FILE* ri_ptr = fopen(out_prefix, "w");

    while(l = kseq_read(kref) >= 0){
        chromo_position = 0;
        chromo_length = kref->seq.l;
        memcpy(ref_name[chromo_cnt], kref->name.s, kref->name.l); 
        ref_name_len[chromo_cnt] = kref->name.l;
        ref_length[chromo_cnt] = kref->seq.l;
        //memcpy(ref_name[chromo_cnt], kref->name.s, kref->name.l); 
        printf("name: %s\n", kref->name.s); 
        printf("name len: %llu\n", kref->name.l); 
        printf("seq len: %llu\n", kref->seq.l);
        fprintf(ri_ptr, "%lu\n", kref->seq.l);
        fprintf(ri_ptr, "%lu\n", kref->name.l);
        fprintf(ri_ptr, "%s\n", kref->name.s);

        //if(kref->comment.l)  printf("comment : %s\n", kref->comment.s);
        //printf("seq length %s\n", kref->seq.s);

        chromo_hash_index(kref->seq.s);

        //printf("seq length %llu %llu\n", chromo_position, kref->seq.l);
        chromo_cnt += 1;
        //if(chromo_cnt == 10)

        //if(kref->qual.l) printf("qual : %s\n", kref->qual.s);
    }
    printf("chrom num %u\n", chromo_cnt);
    //base_2bits[(position + 15)/16 - 1] = enc_2bits;
    base_2bits[2*((position + 31)/32 - 1)] = enc_lbits;
    base_2bits[2*((position + 31)/32 - 1) + 1] = enc_hbits;
    base_3bits[(position + 31)/32 - 1] = enc_3bits;
    seeds_2bits[(position + 15)/16 - 1] = enc_seeds_2bits;

    printf("return value:%d %x\n", l, enc_3bits);
    printf("valid seed size : %llu %u\n", valid_seed_size, position);
    printf("store positions\n");
    write_base_2bits_to_file(position);

    kseq_destroy(kref);
    positions = (sindex_t*)malloc((valid_seed_size + 8) * sizeof(sindex_t)); 
    pos_lseeds = (half_word_t*)malloc((valid_seed_size + 8) * sizeof(half_word_t)); 
    //double mem_size = (double)total_malloc_size / 1024 / 1024 / 1024;
    //printf("memory usage: %.2f GB\n", mem_size);
    store_positions();
    write_positions_to_file(valid_seed_size);
    free(base_2bits);
    fclose(ri_ptr);



}
char nbase[] = { 'A', 'C', 'G', 'T'};
char int2base(index_t ind){
    sindex_t n = base_3bits[ind/32] >> (ind&0x1f);
    if(n & 1){
        return 'N';
    }
    sindex_t l = base_2bits[(ind/32) * 2]  >> (ind&0x1f);
    l &=0x1;
    sindex_t h = base_2bits[(ind/32) * 2 + 1] >> (ind&0x1f);
    h &= 0x1;
    l = (h << 1) | l;

    return nbase[l];
}
char hl_int2base(index_t ind){
    sindex_t n = base_3bits[ind/32] >> (ind&0x1f);
    if(n & 1){
        return 'N';
    }
    sindex_t l = seeds_2bits[ind/16] >> ((ind&0xf) * 2);
    return nbase[l&0x3];
}
void verify_base_2bits(char *ref_seq, sindex_t ref_len, sindex_t ind){
    index_t i;
    for(i = 0; i < ref_len; i ++){
        char m = hl_int2base(ind);
        char n = int2base(ind ++);
        if(m != toupper(ref_seq[i])){
            printf("base 2bits wrong! %u %u %c %c\n",i, ref_len, toupper(ref_seq[i]), m);
            break;
        }
        if(n != toupper(ref_seq[i])){
            printf("seeds 2bits wrong! %u %u %c %c\n",i, ref_len, toupper(ref_seq[i]), n);
            break;
        }
    }

}
//extern SLAVE_FUNC(build_hash_index)();
typedef struct build_info_struct{
    sindex_t *seeds_buffer[2];
    sindex_t *base_2bits;
    sindex_t *base_3bits;
    sindex_t *seeds_2bits;

    char *ref_seq;
    char *res_seq;

    index_t res_len;
    index_t ref_seq_len;

    index_t cur_buffer;  
    index_t total_len;
    index_t seed_length; 
}build_info_t;
typedef struct index_info_h{
    sindex_t *seeds_2bits;
    sindex_t *base_2bits;
    sindex_t *base_3bits;
    sindex_t *seeds_buffer;
    uint16_t *pos_lseeds;
    index_t total_len;
    index_t ref_len;
    index_t seed_length;
    index_t store_flag;
} index_info_t;

extern void slave_build_hash_index(build_info_t *h);
extern void slave_store_hash_index(index_info_t *h);
sindex_t get_right_seeds(char *ref, sindex_t ind){
    index_t i = 0;
    sindex_t seed_right = 0;
    index_t j = 0;
    for(i = ind;  i < ind + seed_length; i ++){
        char c= toupper(ref[i]) ;
        //if(ind == 972486*4)
        //    printf("%c",  c);
        if(c == 'N'){
            seed_right = 0xffffffff; 
            break;
        }
        switch(c){
            case ('A'):
                seed_right |= 0U << (j*2); 
                break;
            case ('C'):
                seed_right |= 1U << (j*2); 
                break;
            case ('G'):
                seed_right |= 2U << (j*2); 
                break;
            case ('T'):
                seed_right |= 3U << (j*2); 
                break;
        }
        j ++;
    }
    return seed_right;
}

sindex_t get_right_pos_lseed(sindex_t pos){

    index_t s_pos = 0;  
    sindex_t res_len = SEED_STEP - 1; 
    sindex_t mask = (1U<<res_len) - 1;
    if(pos >= res_len)
        s_pos = pos - res_len;
    sindex_t shift = s_pos % 32;
    s_pos /= 32;

    sindex_t flseed = (base_2bits[s_pos * 2]  >> shift) |  (base_2bits[(s_pos + 1) * 2] << (32 - shift));
    sindex_t fhseed = (base_2bits[s_pos * 2 + 1]  >> shift) |  (base_2bits[(s_pos + 1) * 2 + 1]  << (32 - shift));
    flseed = (flseed & mask); 
    fhseed = (fhseed & mask);

    s_pos = pos + seed_length;
    shift = s_pos % 32;
    s_pos /= 32;

    sindex_t blseed = (base_2bits[s_pos * 2]  >> shift) |  (base_2bits[(s_pos + 1) * 2]  << (32 - shift));
    sindex_t bhseed = (base_2bits[s_pos * 2 + 1]  >> shift) |  (base_2bits[(s_pos + 1) * 2 + 1]  << (32 - shift));
    blseed = (blseed & mask); 
    bhseed = (bhseed & mask);

    flseed = flseed | (blseed << res_len);
    fhseed = fhseed | (bhseed << res_len);

    sindex_t lseed_enc = flseed | (fhseed << (res_len * 2));
    return lseed_enc;

}

void construct_hash_index_cpe(){


    seed_count = (sindex_t*)malloc((1 << (seed_length * 2))  * sizeof(sindex_t));
    memset(seed_count, 0, (1 << (seed_length * 2)) * sizeof(sindex_t));
    hash_index = (sindex_t*)malloc(((1 << (seed_length * 2)) +  1)  * sizeof(sindex_t)); 
    memset(hash_index, 0, (1 << (seed_length * 2)) * sizeof(sindex_t));

    //positions = (sindex_t*)malloc(MAX_HASH_INDEX * sizeof(sindex_t)); 
    seeds_2bits = (sindex_t*)malloc(780*MB);
    base_3bits = (sindex_t*)malloc(390*MB);

    half_word_t *pos_lseeds_tmp = (half_word_t*) malloc(100*MB * sizeof(half_word_t));
    base_2bits = (sindex_t*)malloc(780*MB);

    sindex_t *seeds_buffer = (sindex_t*)malloc(100 * MB * sizeof(sindex_t));
    sindex_t res_bbuffer[32];

    index_info_t index_info;

    build_info_t build_info;
    build_info.seeds_buffer[0] = seeds_buffer;
    //build_info.seeds_buffer[1] = seeds_buffer + 3 * MB;

    build_info.base_2bits = base_2bits;
    build_info.base_3bits = base_3bits;
    build_info.seeds_2bits = seeds_2bits;
    build_info.res_seq = (char*)res_bbuffer;
    build_info.cur_buffer = 0;
    build_info.total_len = 0;
    build_info.res_len = 0;
    build_info.seed_length = seed_length;


    //enc_2bits = 0;
    enc_lbits = 0;
    enc_hbits = 0;
    enc_seeds_2bits = 0;
    enc_3bits = 0;
    memset(seeds_2bits, 0, 780*MB);
    memset(base_2bits, 0, 780*MB);

    valid_seed_size = 0; 
    position = 0;
    //char *path = "hg38_test.fa";
    //char *path = "hg38/hg38.fa";
    //char *path = "chr1/chr1.fa";
    kseq_t *kref = kseq_open(ref_path);
    int l ;
    cur_seed_length = 0;
    seed = 0;
    chromo_cnt = 0;
    sprintf(out_prefix, "%s%s", prefix, ".ref_info");
    ///FILE* ri_ptr = fopen(OUTDIR"ref_info", "w");
    FILE* ri_ptr = fopen(out_prefix, "w");
//#define DEBUG
#ifdef DEBUG
    index_info.total_len = 0;
    index_t seq_total_len = 0;
#endif
    athread_init();
    while(l = kseq_read(kref) >= 0){
        chromo_position = 0;
        chromo_length = kref->seq.l;
        memcpy(ref_name[chromo_cnt], kref->name.s, kref->name.l); 

        ref_name_len[chromo_cnt] = kref->name.l;
        ref_length[chromo_cnt] = kref->seq.l;

        fprintf(ri_ptr, "%lu\n", kref->seq.l);
        fprintf(ri_ptr, "%lu\n", kref->name.l);
        fprintf(ri_ptr, "%s\n", kref->name.s);

        printf("name: %s\n", kref->name.s); 
        //printf("name len: %llu\n", kref->name.l); 
        printf("seq len: %llu\n", kref->seq.l);

        build_info.ref_seq = kref->seq.s; 
        build_info.ref_seq_len = kref->seq.l;

        athread_spawn(build_hash_index, &build_info);
        athread_join();


        memcpy(build_info.res_seq, kref->seq.s + kref->seq.l - 32, 32);

        build_info.total_len += kref->seq.l;
        build_info.res_len = build_info.total_len  & 0x1f;

        //printf("total len %llu\n", build_info.total_len);

        index_t seed_num = (kref->seq.l - seed_length + 1) / SEED_STEP;
        index_t i = 0;
        for(i = 0; i < seed_num; i ++){
            if(seeds_buffer[i] != SINVALID_SEED){
                seed_count[seeds_buffer[i]] += 1;
            }
        }

#ifdef DEBUG


        verify_base_2bits(kref->seq.s, kref->seq.l, seq_total_len);
        seq_total_len += kref->seq.l;

        index_info.seeds_2bits = seeds_2bits;
        index_info.base_2bits = base_2bits;
        index_info.base_3bits = base_3bits;

        index_info.seeds_buffer =  seeds_buffer;
        index_info.seed_length = seed_length;
        index_info.pos_lseeds = pos_lseeds_tmp;
        index_info.ref_len = ref_length[chromo_cnt];

        index_info.store_flag = 2;
        athread_spawn(store_hash_index, &index_info);
        athread_join();


        
        for(i = 0; i < seed_num; i ++){
            sindex_t right_seed = get_right_seeds(kref->seq.s, i * 4);
            index_t offset = index_info.total_len;
            if(seeds_buffer[i] != right_seed){
                printf("seed wrong %u %x %x\n", i, seeds_buffer[i], right_seed);
                break;
                if(right_seed != SINVALID_SEED){
                    sindex_t right_lseed = get_right_pos_lseed(offset + i*SEED_STEP);
                    sindex_t clseed = pos_lseeds_tmp[i];
                    if(clseed != right_seed){
                        printf("pos lseed wrong %u %x %x\n", i, right_lseed, clseed);
                        break;
                    }
                }

            }

        }

        index_info.total_len += ref_length[chromo_cnt];
#endif


        chromo_cnt += 1;
    }
    printf("chrom num %u\n", chromo_cnt);

    fclose(ri_ptr);

    position = build_info.total_len; 
    write_base_2bits_to_file(position);
    kseq_destroy(kref);

    valid_seed_size = 0;
    index_t i  = 0; 
    for( i = 0; i< (1 << (2*seed_length)); i ++){
        valid_seed_size += seed_count[i];
    }
    printf("valid seed size : %llu\n", valid_seed_size);
    printf("process positions\n");

    //index_info_t index_info;
    index_info.seeds_2bits = seeds_2bits;
    index_info.base_2bits = base_2bits;
    index_info.base_3bits = base_3bits;
    index_info.pos_lseeds = pos_lseeds_tmp;

    index_info.seeds_buffer =  seeds_buffer;
    index_info.seed_length = seed_length;
    index_info.total_len = 0;

    positions = (sindex_t*)malloc((valid_seed_size + 8) * sizeof(sindex_t)); 

    hash_index[0] = 0; 
    for(i = 0; i < (1 << (seed_length * 2)); i ++){
        hash_index[i+1] = hash_index[i] + seed_count[i];
    }

   index_info.store_flag = 1;
    for(i = 0; i < chromo_cnt; i ++){
        index_info.ref_len = ref_length[i];
        athread_spawn(store_hash_index, &index_info);
        athread_join();
        index_t seed_num = (ref_length[i] - seed_length + 1) / SEED_STEP;
        sindex_t offset = index_info.total_len;
        index_info.total_len += ref_length[i];
        index_t j = 0;
        for(j = 0; j < seed_num; j ++){
            sindex_t seed = seeds_buffer[j];
            if(seed != SINVALID_SEED){
                index_t ind = hash_index[seed];
                positions[ind] = offset + j * SEED_STEP;
                hash_index[seed] ++;
            }
        }
        printf(".");
    }
    printf("\n");
    write_positions_to_file(valid_seed_size);
    free(positions);

    printf("process pos_lseeds\n");

    index_info.store_flag = 2;
    index_info.total_len = 0;
    pos_lseeds = (half_word_t*)malloc((valid_seed_size + 8) * sizeof(half_word_t)); 
    hash_index[0] = 0; 
    for(i = 0; i < (1 << (seed_length * 2)); i ++){
        hash_index[i+1] = hash_index[i] + seed_count[i];
    }
    for(i = 0; i < chromo_cnt; i ++){
        index_info.ref_len = ref_length[i];
        athread_spawn(store_hash_index, &index_info);
        athread_join();
        index_info.total_len += ref_length[i];
        index_t seed_num = (ref_length[i] - seed_length + 1) / SEED_STEP;
        index_t j = 0;
        for(j = 0; j < seed_num; j ++){
            sindex_t seed = seeds_buffer[j];
            if(seed != SINVALID_SEED){
                index_t ind = hash_index[seed];
                pos_lseeds[ind] = pos_lseeds_tmp[j];
                hash_index[seed] ++;
            }
        }

        printf(".");
    }
    printf("\n");
    free(base_2bits);
    double mem_size = (double)total_malloc_size / 1024 / 1024 / 1024;
    printf("memory usage: %.2f GB\n", mem_size);
}
void printf_help(){
    printf("-i\tinput reference\n");
    printf("-p\tindex prefix\n");
    printf("-s\t[opt, default:12] seed length\n");
    exit(0);
}
//#define MPE_PROFILE
int main(int argc, char *argv[]){
    char *optstring  = "h::i:s::p:";
    int opt;
#ifdef MPE_PROFILE
    MPI_Init(NULL, NULL);
    GPTLinitialize();

    fem_swlu_prof_init();
    fem_swlu_prof_init();
    fem_swlu_prof_start();
    GPTLstart("total time");
#endif
    seed_length = 12;
    while((opt = getopt(argc,argv,optstring)) != -1){
        switch(opt){
            case ('i'):
                memcpy(ref_path, optarg, strlen(optarg));
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
        }
    }
#ifdef HOST_BUILD
    construct_hash_index();
#else
    construct_hash_index_cpe();
#endif
    write_to_file();

#ifdef MPE_PROFILE
    fem_swlu_prof_stop();
    fem_swlu_prof_print();
#endif

#ifdef MPE_PROFILE
    GPTLstop("total time");
    GPTLpr(0);

    MPI_Finalize();
#endif
}
