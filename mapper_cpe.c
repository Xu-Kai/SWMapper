#include<slave.h>
#include<simd.h>
#include<stdint.h>
#include<dma.h>
#include<assert.h>
#include "fm_define.h"
#include "/home/export/online1/swmore/opensource/cal/cal.h"
#define DMA_FAST

//#define DMA_COUNT_BYTES
#include "dma_macros.h"

//#define CPE_PROFILE
//#define INVALID_FLAG 0xffffffffffffffff
#define READ_LEN 256
#define SIMD_WIDTH 8
#define MAXERROR 10000

#define word_t uint32_t 
#define WORD_BITS 32
#define READ_INT_LEN (256/16)
#define PRINT_CORE 16

//#define NO_ALIGN
#define BALANCE
#define DMA_ASYN
#define LSEED_FILTER
//#define LONG_SEED
#define SORT_POS
#define SIMD_BM
#define SIMD_BM_ASM
#define OUT_CIGAR


/*
#define index_t uint64_t
#define sindex_t uint32_t
#define half_word_t uint16_t
#define SEED_STEP 4
#define INVALID_SFLAG 0xffffffff
#define MAX_CIGAR 256
#define CIGAR_LEN 16
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
typedef struct hash_count{
    sindex_t hash_count;
    sindex_t index;
}hash_count_t;

typedef struct hash_steps{
    sindex_t step;
    sindex_t index;
    sindex_t pos_size;
}hash_steps_t;
typedef struct seed_heap{
    sindex_t index;
    sindex_t pos_size;
}seed_heap_t;


typedef struct error_s{
    int error;
    int location;
}error_t;
typedef struct cigar_str{
    //char cigar[15*MAX_CIGAR];
    char *cigar;
    //sindex_t location[MAX_CIGAR];
    sindex_t *location;
    sindex_t cigar_len;
    sindex_t cigar_num;
    index_t main_flag;
    index_t last_flag;
    index_t f_flag;
}cigar_t;
typedef struct seed_ind_s{
    sindex_t loc;
    sindex_t seed_ind;
}seed_ind_t;

typedef struct read_info_s{
    char *read;
    sindex_t *read_2bits;
    sindex_t *read_seeds;
    uint8_t *res_seeds;
    sindex_t read_length;
}read_info_t;

#define shrw(x, b) ((x) >> (b))
#define shlw(x, b) ((x) << (b))
#define andw(x, y) ((x) & (y))
#define xorw(x, y) ((x) ^ (y))
#define orw(x, y) ((x) | (y))
#define notw(x) (~(x))
#define subw(x, y) ((x) - (y))
#define addw(x, y) ((x) + (y))
#define eqw(x, y) (~((x) ^ (y)))
#define ZEROW 0
#define ONEW 1

#ifdef CPE_PROFILE
#define CPE
#define LWPF_KERNELS K(choose_seed_pf) K(alignment_pf) K(faaw_pf) K(bmyers_pf) K(heap_pf) K(pnum_pf) K(heap_sort_pf) \
    K(align_pf) K(myers_pf) K(load_refc_pf) K(shift_pf) K(cigar_pf) \
    K(load_pf) K(simd_pf) K(loopf_pf) K(loopb_pf) 
#define LWPF_UNIT U(MAPPING) 
#include "lwpf2/lwpf2.h"
#endif

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


void print_int2char(sindex_t con){
    index_t i = 0;
    for( i = 0; i < 16; i ++){
        index_t tmp = con & 0x3;
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
        con >>= 2;
    }

}

void print_hl2char(sindex_t con){
    index_t i = 0;
    for( i = 0; i < 16; i ++){
        index_t tmp = (con >> i) & 0x1;
        index_t tmp1 = (con >> (16 + i)) << 1;
        tmp = (tmp | tmp1) & 0x3;
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
inline int sprint_length(char *cigar, sindex_t rlen){
    int llen = 0;
    sindex_t  rlen_tmp = rlen;
    while(rlen_tmp != 0){
        llen += 1;
        rlen_tmp /= 10;
    }
    int lres = llen;
    while(rlen != 0){
        char c = '0' + (char)rlen%10;
        cigar[--llen] = c;
        rlen /= 10;
    }
    return lres; 
}
void align_reads(ref_info_t *ref_info_s, sindex_t *reads_2bits, sindex_t *ref_2bits, 
        sindex_t *ref_length, sindex_t read_length, index_t *ref_shift, index_t *ref_pos, 
        cigar_t *cigar, error_t *err_loc_info, index_t valid_size){
#ifdef NO_ALIGN
return;
#endif
    int iii, ii, i;
    char *cigar_s;
    int cigar_len = 0;
    sindex_t d0_arry_64[READ_LEN];
    sindex_t hp_arry_64[READ_LEN];
    int route_size_whole[READ_LEN];
    char route_char_whole[READ_LEN];

    for(iii = 0; iii < valid_size; iii ++){
        if(err_loc_info[iii].error <= ref_info_s->errors){
            cigar_len = 0; 
            cigar_s = cigar->cigar + cigar->cigar_len; 
            int start_location = err_loc_info[iii].location - read_length + 1;
            int tmp_err = 0;
            sindex_t hmask = 0xffff0000;
            sindex_t lmask = 0x0000ffff;
            //printf("start loc %u %u\n", start_location, ref_pos[iii] - 4);
            for(ii = 0; ii * WORD_BITS < read_length; ii ++){

                int res = WORD_BITS;
                sindex_t mask = 0xffffffff;
                sindex_t neq_mask = 0x1;
                if((ii + 1) * WORD_BITS > read_length){
                    res = read_length - ii * WORD_BITS;
                    mask >>= (32 - res); 
                }

                sindex_t lref = (ref_2bits[(ii*2)*SIMD_WIDTH + iii]  >> start_location)| (ref_2bits[(ii + 1)*2*SIMD_WIDTH + iii] << (32 - start_location));
                sindex_t href = (ref_2bits[(ii*2 + 1)*SIMD_WIDTH + iii]  >> start_location)| (ref_2bits[((ii + 1)*2 + 1)*SIMD_WIDTH + iii] << (32 - start_location));


                sindex_t lread = reads_2bits[READ_INT_LEN + 2*ii];
                sindex_t hread = reads_2bits[READ_INT_LEN + 2*ii + 1];
                sindex_t xor_res = (hread ^ href) | (lread ^ lref);
                xor_res &=  mask;
                while(xor_res != 0){
                    if(xor_res & neq_mask != 0){
                        tmp_err ++; 
                    }
                    xor_res >>= 1;
                }
            }
            if(tmp_err == err_loc_info[iii].error){
                cigar_s[cigar_len++] = (char)tmp_err; 
                cigar_len ++;

                cigar_len += sprint_length(cigar_s + cigar_len, read_length);
                cigar_s[cigar_len++] = 'M'; 
                cigar_s[1] = (char)(cigar_len - 2);
                cigar->cigar_len += cigar_len;

                cigar->location[cigar->cigar_num++] = ref_pos[iii] + start_location + 1;
            }
            if(tmp_err != err_loc_info[iii].error){
                int dist_error = (int) ref_info_s->errors;
                word_t band_length = 2 * (ref_info_s->errors) + 1;
                word_t band_down = 2 * (ref_info_s->errors);
                word_t peq[4];
                int j;
                word_t res = 0;
                peq[0] = 0;
                peq[1] = 0;
                peq[2] = 0;
                peq[3] = 0;
                word_t tmp = (word_t)1;

                word_t ref_word = 0;
                word_t ref_bits = ref_shift[iii];
                word_t one = 1;
                for(j = 0; j< band_length; j ++){
                    word_t base_enc1 = (ref_2bits[iii] >> (j & 0x1f)) & 0x1; 
                    word_t base_enc2= (ref_2bits[SIMD_WIDTH + iii] >> (j & 0x1f)) & 0x1; 
                    word_t base_enc = base_enc1 | (base_enc2  << 1);

                    peq[base_enc] |= tmp;
                    tmp <<= 1;
                }
                int last_high = band_length - read_length + ref_length[iii] - band_down - 1;
                word_t x = 0,  vp = 0, vn = 0, d0 = 0, hn = 0, hp = 0;
                word_t w_0 = ZEROW, w_1 = ONEW;
                word_t i_bd = band_down;
                word_t err_mask = (word_t)1;
                word_t mask = (word_t)1<<(band_length - 1);
                word_t base_mask = 0x3;
                for(ii = 0; ii * WORD_BITS/2 < read_length; ii ++){
                    res = WORD_BITS/2; 
                    if ((ii + 1) * WORD_BITS / 2 > read_length)
                        res = read_length - ii * WORD_BITS/2;
                    for (i = 0; i < res; i ++){
                        word_t base_enc = (reads_2bits[ii] >> (i * 2)) & base_mask; 
                        //x = match | vn;
                        x = orw(peq[base_enc],  vn);
                        //printf("%x\n", x);
                        //d0 = ((vp + (x & vp)) ^ vp) | x;
                        d0 = andw(x, vp);
                        d0 = addw(d0, vp);
                        d0 = xorw(d0, vp);
                        d0 = orw(d0, x);
                        /* hn = vp & d0; */
                        hn = andw(vp, d0);
                        /* hp = vn | ~(vp | d0); */
                        hp = orw(vp, d0);
                        hp = notw(hp);
                        hp = orw(vn, hp);
                        //y = d0 >> 1;
                        x = shrw(d0, 1);
                        //vn = y & hp;
                        vn = andw(x, hp);
                        //vp = hn | ~(y | hp);
                        vp = orw(x, hp);
                        vp = notw(vp);
                        vp = orw(hn, vp);
                        d0_arry_64[ii * WORD_BITS/2 + i] = d0; 
                        hp_arry_64[ii * WORD_BITS/2 + i] = hp; 
                        peq[0] = peq[0]>>1;
                        peq[1] = peq[1]>>1;
                        peq[2] = peq[2]>>1;
                        peq[3] = peq[3]>>1;
                        ++i_bd;
                        ref_word = i_bd >> 5;

                        word_t base_enc1 = (ref_2bits[(ref_word*2)*SIMD_WIDTH + iii] >> (i_bd & 0x1f)) & 0x1; 
                        word_t base_enc2= (ref_2bits[(ref_word * 2 + 1)*SIMD_WIDTH + iii] >> (i_bd & 0x1f)) & 0x1; 
                        base_enc = base_enc1 | (base_enc2  << 1);

                        peq[base_enc] |= mask;
                    }
                }
                int match_site = err_loc_info[iii].location;
                int site=ref_length[iii]-last_high-1;
                int search_site =  match_site - site;
                int pre_size = 1;
                char pre_char = 'N';
                sindex_t mask_1 = 1;
                i = read_length - 1;
                int sum_err = 0;

                int err = err_loc_info[iii].error;
                //#define ref_read_comp(rfi, rdi, op) (((ref_2bits[(rfi/16)*SIMD_WIDTH + iii] >> (2*(rfi%16)))&0x3) op ((reads_2bits[rdi/16] >> (2*(rdi%16)))&0x3))
#define REF_BITS(rfi) ref_2bits[(rfi)*SIMD_WIDTH + iii]
#define REF_LBIT(rfi) ((REF_BITS( (rfi>>5) * 2) >> (rfi&0x1f)) & 0x1)
#define REF_HBIT(rfi) (((REF_BITS( (rfi>>5) * 2 + 1) >> (rfi&0x1f)) & 0x1) << 1)

#define ref_read_comp(rfi, rdi, op) ((REF_HBIT(rfi) | REF_LBIT (rfi)) op ((reads_2bits[rdi/16] >> (2*(rdi%16)))&0x3))
                j=1;
                if(((d0_arry_64[i]>>search_site)&mask_1)&&(ref_read_comp(match_site, i, ==))){
                    i--;
                    pre_size=1;
                    pre_char='M';
                    match_site--;
                }
                else if(!((d0_arry_64[i]>>search_site)&mask_1)&&(ref_read_comp(match_site, i, !=))){
                    i--;
                    pre_size=1;
                    pre_char='S';
                    match_site--;
                    sum_err++;
                }
                else if((hp_arry_64[i]>>search_site)&mask_1){
                    i--;
                    search_site++;
                    pre_size=1;

                    pre_char='S';
                    sum_err++;
                    start_location++;
                }
                else{
                    search_site--;
                    pre_size=1;
                    pre_char='D';

                    match_site--;
                    sum_err++;
                    start_location--;
                }

                while(i>=0){
                    if(sum_err==err) break;

                    if(((d0_arry_64[i]>>search_site)&mask_1)&&(ref_read_comp(match_site, i, ==))){
                        i--;
                        match_site--;

                        if(pre_char!='M'){
                            route_size_whole[j]=pre_size;
                            route_char_whole[j++]=pre_char;
                            pre_size=1;
                            pre_char='M';
                        }
                        else pre_size++;
                    }
                    else if(!((d0_arry_64[i]>>search_site)&mask_1)&&(ref_read_comp(match_site, i, !=))){
                        i--;
                        match_site--;
                        sum_err++;

                        if(pre_char!='S'){
                            route_size_whole[j]=pre_size;
                            route_char_whole[j++]=pre_char;
                            pre_size=1;
                            pre_char='S';
                        }
                        else pre_size++;
                    }
                    else if((hp_arry_64[i]>>search_site)&mask_1){
                        i--;
                        search_site++;
                        sum_err++;
                        if(pre_char!='I'){
                            route_size_whole[j]=pre_size;
                            route_char_whole[j++]=pre_char;
                            pre_size=1;
                            pre_char='I';
                        }
                        else pre_size++;
                        start_location++;
                    }
                    else{
                        search_site--;
                        match_site--;
                        sum_err++;
                        if(pre_char!='D'){
                            route_size_whole[j]=pre_size;
                            route_char_whole[j++]=pre_char;
                            pre_size=1;
                            pre_char='D';
                        }
                        else pre_size++;
                        start_location--;
                    }
                }

                route_size_whole[j]=pre_size;
                route_char_whole[j++]=pre_char;
                if(i>=0){
                    route_size_whole[j]=i+1;
                    route_char_whole[j++]='M';
                }

                if(cigar->cigar_num > 0)
                    if(cigar->location[cigar->cigar_num - 1] == ref_pos[iii] + start_location + 1)
                        continue;

                int size_sm=0;

                cigar_s[cigar_len++] = (char)(err_loc_info[iii].error); 
                cigar_len ++;

                cigar->location[cigar->cigar_num++] = ref_pos[iii] + start_location + 1;

                for(j=j-1; j > 0; j--){
                    if(route_char_whole[j]=='M'||route_char_whole[j]=='S'){

                        size_sm=0;
                        while(j > 0 && (route_char_whole[j]=='M'||route_char_whole[j]=='S')){
                            size_sm=size_sm+route_size_whole[j];
                            j--;
                        }
                        j++;
                        //sprintf(cigar+strlen(cigar),"%d%c",size_sm,'m');
                        cigar_len += sprint_length(cigar_s + cigar_len, size_sm);
                        cigar_s[cigar_len++]='M';
                    }
                    else{
                        //sprintf(cigar+strlen(cigar),"%d%c",route_size_whole[j],route_char_whole[j]);
                        cigar_len += sprint_length(cigar_s + cigar_len, route_size_whole[j]);
                        cigar_s[cigar_len++]= route_char_whole[j];
                    }
                }

                cigar_s[1] = (char)(cigar_len - 2);
                cigar->cigar_len += cigar_len;

            }

#undef ref_read_comp
#undef REF_BITS
#undef REF_HBIT
#undef REF_LBIT
        }
    }
}

void banded_myers_simd(ref_info_t *ref_info_s, sindex_t *reads_2bits, sindex_t *ref_2bits, 
        sindex_t *ref_length, sindex_t read_length, index_t *ref_shift, index_t *ref_pos, 
        error_t *err_loc_info, index_t valid_size){

    int dist_error = (int) ref_info_s->errors;
    word_t band_length = 2 * (ref_info_s->errors) + 1;
    word_t band_down = 2 * (ref_info_s->errors);
    word_t peq[4];

    int iii, ii, i;
    int j;
    word_t res = 0;
    char base[4];
    base[0] = 'A';
    base[1] = 'C';
    base[2] = 'G';
    base[3] = 'T';



    for(iii = 0; iii < valid_size; iii ++){
        peq[0] = 0;
        peq[1] = 0;
        peq[2] = 0;
        peq[3] = 0;
        word_t tmp = (word_t)1;
        word_t base_mask = 0x3;
        word_t ref_word = 0;
        word_t ref_bits = ref_shift[iii];


            for(j = 0; j< band_length; j ++){
                word_t base_enc1 = (ref_2bits[iii] >> (j & 0x1f)) & 0x1; 
                word_t base_enc2= (ref_2bits[SIMD_WIDTH + iii] >> (j & 0x1f)) & 0x1; 
                word_t base_enc = base_enc1 | (base_enc2  << 1);

                peq[base_enc] |= tmp;
                tmp <<= 1;
            }

            int last_high = band_length - read_length + ref_length[iii] - band_down - 1;
            //printf("%d\n", last_high);
            word_t x = 0,  vp = 0, vn = 0, d0 = 0, hn = 0, hp = 0;
            //int err = read_length;
            word_t w_0 = ZEROW, w_1 = ONEW;
            int err = 0;
            word_t i_bd = band_down;

            word_t err_mask = (word_t)1;
            word_t pass_flag = 1;
            word_t mask = (word_t)1<<(band_length - 1);
            for(ii = 0; ii * WORD_BITS/2 < read_length; ii ++){

                res = WORD_BITS/2; 
                if ((ii + 1) * WORD_BITS / 2 > read_length)
                    res = read_length - ii * WORD_BITS/2;
                for (i = 0; i < res; i ++){
                    word_t base_enc = (reads_2bits[ii] >> (i * 2)) & base_mask; 

                    //x = match | vn;
                    x = orw(peq[base_enc],  vn);
                    //if(_MYID == PRINT_CORE){
                    //    //if(ii == 0 && i == 0)
                    //    //printf("%x ", peq[base_enc]);
                    //    //printf("%x ", x);
                    //    
                    //}

                    //printf("%x\n", x);
                    //d0 = ((vp + (x & vp)) ^ vp) | x;
                    d0 = andw(x, vp);
                    d0 = addw(d0, vp);
                    d0 = xorw(d0, vp);
                    d0 = orw(d0, x);

                    /* hn = vp & d0; */
                    hn = andw(vp, d0);

                    /* hp = vn | ~(vp | d0); */
                    hp = orw(vp, d0);
                    hp = notw(hp);
                    hp = orw(vn, hp);

                    //y = d0 >> 1;
                    x = shrw(d0, 1);

                    //vn = y & hp;
                    vn = andw(x, hp);

                    //vp = hn | ~(y | hp);
                    vp = orw(x, hp);
                    vp = notw(vp);
                    vp = orw(hn, vp);
                    if(!(d0&err_mask)){
                        ++err;
                    }

                    peq[0] = peq[0]>>1;
                    peq[1] = peq[1]>>1;
                    peq[2] = peq[2]>>1;
                    peq[3] = peq[3]>>1;
                    ++i_bd;


                    ref_word = i_bd >> 5;
                    //base_enc = (ref_2bits[ref_word*SIMD_WIDTH + iii] >> ( (i_bd % 16) * 2)) & base_mask;

                    word_t base_enc1 = (ref_2bits[(ref_word*2)*SIMD_WIDTH + iii] >> (i_bd & 0x1f)) & 0x1; 
                    word_t base_enc2= (ref_2bits[(ref_word * 2 + 1)*SIMD_WIDTH + iii] >> (i_bd & 0x1f)) & 0x1; 
                    base_enc = base_enc1 | (base_enc2  << 1);

                    peq[base_enc] |= mask;
               }
            }

           int site=ref_length[iii]-last_high-1;

            int location = -1;
            int error = MAXERROR;
            if((err <= dist_error)&&(err < error)){
                error = err;
                location = site;
            }

            for(i = 0; i < last_high; i++){
                err=err+((vp>>i)&(word_t)1);
                err=err-((vn>>i)&(word_t)1);
                if((err<=dist_error)&&(err<error)){
                    error = err;
                    location = site+i+1;
                }
            }
            if(error > dist_error){
                pass_flag = 0;
            }
            err_loc_info[iii].error = error;
            err_loc_info[iii].location = location;


    }

}

void print_simd_n(uintv8 in, sindex_t i){
    uint32_t v[8];
    simd_store(in, v);
    printf("%x\n", v[i]);
}
//#define FIND 2907861898 
void banded_myers_simd_asmv2(ref_info_t *ref_info_s, sindex_t *reads_2bits, sindex_t *ref_2bits, 
        sindex_t *ref_length, sindex_t read_length, index_t *ref_shift, index_t *ref_pos, 
        error_t *err_loc_info, index_t valid_size){

    int dist_error = (int)ref_info_s->errors;
    sindex_t band_length = 2 * (ref_info_s->errors) + 1;
    sindex_t band_down = 2 * (ref_info_s->errors);
    //uintv8  peqv[4];
    //    uintv8 ref_2bitsv[READ_LEN/16];
    uintv8 onev = 1;

    int iii, ii, i;
    uintv8 base_maskv = 0x3;
    uintv8 peqv_0 = 0;
    uintv8 peqv_1 = 0;
    sindex_t peq_mask = (1U << band_length) - 1;
    uintv8 peqv_mask = peq_mask;
   //sindex_t band_down_mask = (1U << band_down) - 1;
    if(band_length < 32){
        uintv8 tmp1;
        uintv8 tmp2;
        simd_load(tmp1, ref_2bits);
        simd_load(tmp2, ref_2bits + SIMD_WIDTH);
        peqv_0 = tmp1 & peqv_mask;
        peqv_1 = tmp2 & peqv_mask;
    }

    uintv8 ref_2bitsv[READ_LEN / 16];
    for(i = 0; i < (read_length + band_down + 31)/32; i ++){

        uintv8 tmp1;
        uintv8 tmp2;

        uintv8 tmp3;
        uintv8 tmp4;

        sindex_t r_shift = band_length;
        sindex_t l_shift = 32 - band_length;

        simd_load(tmp1, ref_2bits + (i * 2    )*SIMD_WIDTH);
        simd_load(tmp2, ref_2bits + ((i + 1)* 2)*SIMD_WIDTH);
        //simd_print_int256(tmp1);
        //simd_print_int256(tmp2);
        asm ("vsrlw  %2, %1, %0\n":"=r"(tmp1):"r"(r_shift), "r"(tmp1)); 
        asm ("vsllw  %2, %1, %0\n":"=r"(tmp2):"r"(l_shift), "r"(tmp2)); 
        tmp3 =  simd_vbisw(tmp1,tmp2);

        //simd_print_int256(tmp3);

        simd_load(tmp1, ref_2bits + (i * 2 + 1)*SIMD_WIDTH);
        simd_load(tmp2, ref_2bits + ((i + 1) * 2 + 1)*SIMD_WIDTH);

        //simd_print_int256(tmp1);
        //simd_print_int256(tmp2);

        asm ("vsrlw  %2, %1, %0\n":"=r"(tmp1):"r"(r_shift), "r"(tmp1)); 
        asm ("vsllw  %2, %1, %0\n":"=r"(tmp2):"r"(l_shift), "r"(tmp2)); 
        tmp4 =  simd_vbisw(tmp1,tmp2);

        //simd_print_int256(tmp4);

        uintv8 hmaskv = 0xffff0000;
        uintv8 lmaskv = 0x0000ffff;

        tmp1 = simd_vandw(tmp3, lmaskv);
        tmp2 = simd_vandw(tmp4, lmaskv);

        //simd_print_int256(tmp1);
        //simd_print_int256(tmp2);

        asm ("vsllw  %1, 16, %0\n":"=r"(tmp2):"r"(tmp2)); 
        tmp1 = simd_vbisw(tmp1, tmp2);
        ref_2bitsv[i*2] = tmp1;

        //simd_print_int256(tmp1);

        tmp1 = simd_vandw(tmp3, hmaskv);
        tmp2 = simd_vandw(tmp4, hmaskv);

        //simd_print_int256(tmp1);
        //simd_print_int256(tmp2);

        asm ("vsrlw  %1, 16, %0\n":"=r"(tmp1):"r"(tmp1)); 
        tmp1 = simd_vbisw(tmp1, tmp2);
        ref_2bitsv[i*2 + 1] = tmp1;

        //simd_print_int256(tmp1);

    }

#ifdef CPE_PROFILE
    lwpf_start(loopf_pf);
#endif
    sindex_t base_mask = 0x3;
    sindex_t ref_word = 0;

    uintv8 xv = 0,  vpv = 0, vnv = 0, d0v = 0, hnv = 0, hpv = 0;
    intv8 errv = 0;

    uintv8 err_maskv = (word_t)1;
    sindex_t i_bd = band_down; 

#ifdef CPE_PROFILE
    lwpf_stop(loopf_pf);
#endif
    //index_t pre_ref = 0;
    //uintv8 base_refv;
    //simd_load(base_refv, ref_2bits + (i_bd/16)*SIMD_WIDTH);

#ifdef CPE_PROFILE
    lwpf_start(simd_pf);
#endif
    sindex_t shift_band = band_length - 1;

    struct {
        uintv8 peqv_l, peqv_h, vp, vn, errv, peqv_mask, refv;
        int i_bd, i2, read, ref_bptr, shift_band, res, pad0, pad1;
    } loop_vars;

    loop_vars.peqv_l = peqv_0;
    loop_vars.peqv_h = peqv_1;

    loop_vars.vp = vpv;
    loop_vars.vn = vnv;

    loop_vars.errv = errv;
    loop_vars.peqv_mask = peqv_mask;

    loop_vars.shift_band = shift_band;
    loop_vars.i_bd = i_bd;
    loop_vars.ref_bptr = ref_2bits;


    asm ("wcsr %0, 0x80\n\t" : : "r"(0xAC));
    asm ("wcsr %0, 0x81\n\t" : : "r"(0x80));
    asm ("wcsr %0, 0x82\n\t" : : "r"(0xBE));
    asm ("wcsr %0, 0x83\n\t" : : "r"(0xF1));
    //sindex_t rshift = 32 - band_length; 
    //uintv8 refh = 0;
    //uintv8 refl = 0; 
    for(ii = 0; ii * WORD_BITS/2 < read_length; ii ++){
        sindex_t res = WORD_BITS/2; 
        if ((ii + 1) * WORD_BITS / 2 > read_length)
            res = read_length - ii * WORD_BITS/2;
        loop_vars.i2 = 0;
        loop_vars.refv = ref_2bitsv[ii];
        loop_vars.read = reads_2bits[ii];
        //for(iii = 0; iii < res ; iii ++){

        loop_vars.res = res;
        uintv8 t0, t1, t2, t3, t4, t5, t6, t7, ptrs;

        asm(//load pointers to PTRS
                "vldd    %[ptrs], 224(%[bptr])\n\t"
                //peqv in t2/3
                "vldd    %[t3], 0(%[bptr])\n\t"
                "vldd    %[t4], 32(%[bptr])\n\t"
                "vldd    %[t6], 64(%[bptr])\n\t"
                "vldd    %[t7], 96(%[bptr])\n\t"
                //ref in at
                "vldd    $at,  192(%[bptr])\n\t"
                //sp for peqv_mask, bptr for errv
                "vinsw   $30, %[ptrs], 6, %[ptrs]\n\t"
                "vinsw   %[bptr], %[ptrs], 7, %[ptrs]\n\t"
                "vldd    $30, 160(%[bptr])\n\t"
                "vldd    %[bptr], 128(%[bptr])\n\t"
                "1001:"
                //read2bits to t0, 2i to t1
                "vextw   %[ptrs], 2, %[t0]\n\t"
                "#vextw   %[ptrs], 1, %[t1]\n\t"
                //benc to t0, extend to vec"
                "#ldi     %[t1], 2(%[t1])\n\t"
                "vshfw   %[t0], %[t0], $31, %[t0]\n\t"
                //peqv in t3/4
                /* "vbis    %[t3], %[t3], %[t4]\n\t" */
                /* "vbis    %[t2], %[t2], %[t3]\n\t" */
                //read_1bitsv in t1/2
                "vand    %[t0], 1, %[t1]\n\t"
                "vsrlw   %[t0], 1, %[t2]\n\t"
                //shift the read
                "srl     %[t0], 2, %[t0]\n\t"
                "vinsw   %[t0], %[ptrs], 2, %[ptrs]\n\t"
                "vsubw   $31, %[t1], %[t1]\n\t"
                "vand    %[t2], 1, %[t2]\n\t"
                "vsubw   $31, %[t2], %[t2]\n\t"
                //t1, t2 for beqvl/h
                "veqv    %[t1], %[t3], %[t1]\n\t"
                "veqv    %[t2], %[t4], %[t2]\n\t"
                //peqv_mask in sp
                "vlog3r1 %[t1], %[t2], $30, %[t0]\n\t"
                //t0->eqv t7->vn, t1->x0
                "vbis    %[t0], %[t7], %[t1]\n\t"
                //t0->d0, t6->vp, t1->x0
                "vand    %[t6], %[t1], %[t0]\n\t"
                "vaddw   %[t0], %[t6], %[t0]\n\t"
                "vlog3r2 %[t0], %[t6], %[t1], %[t0]\n\t"
                //t1->hp, t2->hn, t6->vp, t7->vn, t0->d0
                "vand    %[t6], %[t0], %[t2]\n\t"
                "vlog3r3 %[t7], %[t6], %[t0], %[t1]\n\t"
                //update vp, vn, t5->x0
                "vsrlw   %[t0], 1, %[t5]\n\t"
                "vand    %[t5], %[t1], %[t7]\n\t"
                "vlog3r3 %[t2], %[t5], %[t1], %[t6]\n\t"
                "vlog2x2 %[t0], 1, %[t0]\n\t"
                "vaddw   %[bptr], %[t0], %[bptr]\n\t"
                //ires->t0, bencv->t1
                "#vextw   %[ptrs], 0, %[t1]\n\t"
                //increase i_bd by 1
                "#ldi     %[t1], 1(%[t1])\n\t"
                "#vinsw   %[t1], %[ptrs], 0, %[ptrs]\n\t"
                "#srl     %[t1], 4, %[t1]\n\t"
                "#sll     %[t1], 8, %[t1]\n\t"
                "#vextw   %[ptrs], 3, %[t0]\n\t"
                "#addl    %[t0], %[t1], %[t1]\n\t"
                "#vldd    %[t1], 0(%[t1])\n\t"
                "#and     %[ptrs], 0xf, %[t0]\n\t"
                //shift peqvl/h
                "vsrlw   %[t3], 1, %[t3]\n\t"
                "vsrlw   %[t4], 1, %[t4]\n\t"
                //shift bencvl/h t1/t0
                "#vsrlw   %[t1], %[t0], %[t1]\n\t"
                "vsrlw   $at, 16, %[t0]\n\t"
                "#halt\n\t"
                //extract shift_band to t0\n\t"
                "vand    $at, 1, %[t1]\n\t"
                "vand    %[t0], 1, %[t0]\n\t"
                "vsrlw   $at, 1, $at\n\t"
                //"vextw   %[ptrs], 4, %[t3]\n\t"
                //"vsllw   %[t1], %[t3], %[t1]\n\t"
                //"vsllw   %[t0], %[t3], %[t0]\n\t"
                //extract shift_band to t2\n\t
                "vextw   %[ptrs], 4, %[t2]\n\t"
                "vsllw   %[t1], %[t2], %[t1]\n\t"
                "vsllw   %[t0], %[t2], %[t0]\n\t"
                "vbis    %[t3], %[t1], %[t3]\n\t"
                "vbis    %[t4], %[t0], %[t4]\n\t"

                //extract res, consider jump back
                "vextw   %[ptrs], 5, %[t0]\n\t"
                "ldi     %[t0], -1(%[t0])\n\t"
                "vinsw   %[t0], %[ptrs], 5, %[ptrs]\n\t"
                "bgt     %[t0], 1001b\n\t"
                //store/restore values
                "vextw   %[ptrs], 7, %[t0]\n\t"
                "vstd    %[t3], 0(%[t0])\n\t"
                "vstd    %[t4], 32(%[t0])\n\t"
                "vstd    %[t6], 64(%[t0])\n\t"
                "vstd    %[t7], 96(%[t0])\n\t"
                "vstd    %[bptr], 128(%[t0])\n\t"
                "vextw   %[ptrs], 6, $30\n\t"
                "vinsw   $30, $31, 0, $30\n\t"
                "vinsw   %[bptr], $31, 0, %[bptr]\n\t"
                : [t0]"=&r"(t0), [t1]"=&r"(t1), [t2]"=&r"(t2), [t3]"=&r"(t3),
            [t4]"=&r"(t4), [t5]"=&r"(t5), [t6]"=&r"(t6), [t7]"=&r"(t7),
            [ptrs]"=&r"(ptrs)
                : [bptr]"r"(&loop_vars)
                   : "memory"
                      );

    }

#ifdef CPE_PROFILE
    lwpf_stop(simd_pf);

    lwpf_start(loopb_pf);
#endif
    int err_arr[SIMD_WIDTH];
    sindex_t vp_arr[SIMD_WIDTH];
    sindex_t vn_arr[SIMD_WIDTH];

    simd_store(loop_vars.errv, err_arr);
    simd_store(loop_vars.vp, vp_arr);
    simd_store(loop_vars.vn, vn_arr);

    for(iii = 0; iii < valid_size; iii ++){
        int last_high = band_length - read_length + ref_length[iii] - band_down - 1;
        int site=ref_length[iii]-last_high-1;
        int err = err_arr[iii];
        sindex_t vp = vp_arr[iii]; 
        sindex_t vn = vn_arr[iii];

        //if(ref_pos[iii] ==  FIND - 4){
        //    printf("error %u\n", err);
        //}


        int location = -1;
        int error = MAXERROR;
        if((err <= dist_error)&&(err < error)){
            error = err;
            location = site;
        }

        for(i = 0; i < last_high; i++){
            err=err+((vp>>i)&(word_t)1);
            err=err-((vn>>i)&(word_t)1);
            if((err<=dist_error)&&(err<error)){
                error = err;
                location = site+i+1;
            }
        }
       err_loc_info[iii].error = error;
        err_loc_info[iii].location = location;

    }
#ifdef CPE_PROFILE
    lwpf_stop(loopb_pf);
#endif
}

void banded_myers_simd_v(ref_info_t *ref_info_s, sindex_t *reads_2bits, sindex_t *ref_2bits, 
        sindex_t *ref_length, sindex_t read_length, index_t *ref_shift, index_t *ref_pos, 
        error_t *err_loc_info, index_t valid_size){

    int dist_error = (int)ref_info_s->errors;
    sindex_t band_length = 2 * (ref_info_s->errors) + 1;
    sindex_t band_down = 2 * (ref_info_s->errors);
    //uintv8  peqv[4];
    //    uintv8 ref_2bitsv[READ_LEN/16];
    uintv8 onev = 1;


    int iii, ii, i;
    int j;

    uintv8 ref_2bitsvh[READ_LEN/32];
    uintv8 ref_2bitsvl[READ_LEN/32]; 
    for(i = 0; i < (read_length + band_down + 31)/32; i ++){

        uintv8 base_encvl;
        uintv8 base_encvh;

        simd_load(base_encvl, ref_2bits + (i * 2    )*SIMD_WIDTH);
        simd_load(base_encvh, ref_2bits + (i * 2 + 1)*SIMD_WIDTH);

        ref_2bitsvl[i] = base_encvl;
        ref_2bitsvh[i] = base_encvh;

    }

#ifdef CPE_PROFILE
    lwpf_start(loopf_pf);
#endif
    sindex_t base_mask = 0x3;

    sindex_t ref_word = 0;

    uintv8 base_maskv = 0x3;
    uintv8 peqv_0 = 0;
    uintv8 peqv_1 = 0;
    sindex_t peq_mask = (1U << band_length) - 1;
    uintv8 peqv_mask = peq_mask;

    if(band_length < 32){
        peqv_0 = ref_2bitsvl[0] & peqv_mask;
        peqv_1 = ref_2bitsvh[0] & peqv_mask;
    }


    for(i = 0; i < (read_length + band_down + 31)/32; i ++){
        sindex_t r_shift = band_length;
        sindex_t l_shift = 32 - band_length;
        uintv8 tmp1 = ref_2bitsvh[i];
        uintv8 tmp2 = ref_2bitsvh[i + 1];
        asm ("vsrlw  %2, %1, %0\n":"=r"(tmp1):"r"(r_shift), "r"(tmp1)); 
        asm ("vsllw  %2, %1, %0\n":"=r"(tmp2):"r"(l_shift), "r"(tmp2)); 
        ref_2bitsvh[i] = simd_vbisw(tmp1,tmp2);

        tmp1 = ref_2bitsvl[i];
        tmp2 = ref_2bitsvl[i + 1];
        asm ("vsrlw  %2, %1, %0\n":"=r"(tmp1):"r"(r_shift), "r"(tmp1)); 
        asm ("vsllw  %2, %1, %0\n":"=r"(tmp2):"r"(l_shift), "r"(tmp2)); 
        ref_2bitsvl[i] = simd_vbisw(tmp1,tmp2);
    }


    uintv8 xv = 0,  vpv = 0, vnv = 0, d0v = 0, hnv = 0, hpv = 0;
    //int err = read_length;
    uintv8 w_0v = ZEROW, w_1v = ONEW;
    intv8 errv = 0;

    uintv8 err_maskv = (word_t)1;
    sindex_t i_bd = band_down; 
    uintv8 read_1bitv[2];
    read_1bitv[0] = 0;
    read_1bitv[1] = peqv_mask;

#ifdef CPE_PROFILE
    lwpf_stop(loopf_pf);
#endif
    //index_t pre_ref = 0;
    //uintv8 base_refv;
    //simd_load(base_refv, ref_2bits + (i_bd/16)*SIMD_WIDTH);

#ifdef CPE_PROFILE
    lwpf_start(simd_pf);
#endif
    sindex_t shift_band = band_length - 1;
    //+++++++++++++++++
    //uintv8 ref_nhv =  ref_2bitsvh[1];  
    //uintv8 ref_nlv =  ref_2bitsvl[1];  

    //peqv_0 = ref_2bitsvl[0]; 
    //peqv_1 = ref_2bitsvh[0]; 
    //-----------------

    sindex_t rshift = 32 - band_length; 
    j = 0;
    uintv8 refh = 0;
    uintv8 refl = 0; 
    for(ii = 0; ii * WORD_BITS/2 < read_length; ii ++){
        sindex_t res = WORD_BITS/2; 
        if ((ii + 1) * WORD_BITS / 2 > read_length)
            res = read_length - ii * WORD_BITS/2;
        if(ii % 2 == 0){
            refh = ref_2bitsvh[ii/2];
            refl = ref_2bitsvl[ii/2];
        }
        //i = ii & 0xf;
        for (i = 0; i < res; i ++){
            word_t base_enc = (reads_2bits[ii] >> (i * 2))  & base_mask; 
            uintv8 beqv_l = simd_veqvw(read_1bitv[base_enc&1], peqv_0);
            base_enc >>= 1;
            uintv8 beqv_h = simd_veqvw(read_1bitv[base_enc&1], peqv_1);
            /********/
            uintv8 beqv = beqv_l & beqv_h;
            beqv = beqv & peqv_mask;
//#define USEASM
#ifdef USEASM
            asm("vbis   %[eqv], %[vnv], %[x0v]\n\t"
                    "vand   %[vpv], %[x0v], %[d0v]\n\t"
                    "vaddw  %[d0v], %[vpv], %[d0v]\n\t"
                    "vxor   %[d0v], %[vpv], %[d0v]\n\t"
                    "vbis   %[d0v], %[x0v], %[d0v]\n\t"
                    "vand   %[vpv], %[d0v], %[hnv]\n\t"
                    "vbis   %[vpv], %[d0v], %[hpv]\n\t"
                    "vornot %[vnv], %[hpv], %[hpv]\n\t"
                    "vsrlw  %[d0v],      1, %[x0v]\n\t"
                    "vand   %[x0v], %[hpv], %[vnn]\n\t"
                    "vbis   %[x0v], %[hpv], %[vpn]\n\t"
                    "vornot %[hnv], %[vpn], %[vpn]\n\t"
                    : [d0v]"=&r"(d0v), [hnv]"=&r"(hnv), [hpv]"=&r"(hpv),
                    [x0v]"=&r"(xv ), [vnn]"=r"(vnv), [vpn]"=r"(vpv)
                    : [eqv]"r"(beqv), [vnv]"r"(vnv), [vpv]"r"(vpv));
            //x = match | vn;
            //xv = orw(beqv,  vnv);
#else


            //x = match | vn;
            //xv = orw(beqv,  vnv);
            xv = simd_vbisw(beqv, vnv);

            //printf("%x\n", x);
            //d0 = ((vp + (x & vp)) ^ vp) | x;
            //d0v = andw(xv, vpv);
            //d0v = addw(d0v, vpv);
            //d0v = xorw(d0v, vpv);
            //d0v = orw(d0v, xv);
            d0v = simd_vandw(xv, vpv);
            d0v = simd_vaddw(d0v, vpv);
            /********/
            d0v = simd_vxorw(d0v, vpv);
            //replace
            d0v = simd_vbisw(d0v, xv);
            /*******/

            /* hn = vp & d0; */
            //hnv = andw(vpv, d0v);
            hnv = simd_vandw(vpv, d0v);

            /* hp = vn | ~(vp | d0); */
            //hpv = orw(vpv, d0v);
            //hpv = notw(hpv);
            //hpv = orw(vnv, hpv);
            /*************/
            hpv = simd_vbisw(vpv, d0v);
            //hpv = notw(hpv);
            //hpv = orw(vnv, hpv);
            //replace
            hpv = simd_vornotw(vnv, hpv);
            /**************/

            //y = d0 >> 1;
            //xv = shrw(d0v, 1);
            xv = simd_vsrlw(d0v,1); 

            //vn = y & hp;
            //vnv = andw(xv, hpv);
            vnv = simd_vandw(xv, hpv);

            //vp = hn | ~(y | hp);
            //vpv = orw(xv, hpv);
            //vpv = notw(vpv);
            //vpv = orw(hnv, vpv);
            /**************/
            vpv = simd_vbisw(xv, hpv);
            //replace
            vpv = simd_vornotw(hnv, vpv);
            /*************/
#endif

            /***********/
            intv8 tmp_resv = andw(d0v, err_maskv);
            tmp_resv = xorw(tmp_resv, err_maskv);
            /***********/
            errv = errv + tmp_resv;


            uintv8 bh =  refh;
            uintv8 bl =  refl;

            asm ("vsrlw  %1, 1, %0\n":"=r"(peqv_0):"r"(peqv_0)); 
            asm ("vsrlw  %1, 1, %0\n":"=r"(peqv_1):"r"(peqv_1)); 

            asm ("vsllw  %2, %1, %0\n":"=r"(bh):"r"(shift_band),"r"(bh)); 
            asm ("vsllw  %2, %1, %0\n":"=r"(bl):"r"(shift_band),"r"(bl)); 

            peqv_1 = simd_vbisw(peqv_1, bh);
            peqv_0 = simd_vbisw(peqv_0, bl);

            asm ("vsrlw  %1, 1, %0\n":"=r"(refh):"r"(refh)); 
            asm ("vsrlw  %1, 1, %0\n":"=r"(refl):"r"(refl)); 
            peqv_0 = simd_vandw(peqv_0, peqv_mask);
            peqv_1 = simd_vandw(peqv_1, peqv_mask);
        }


    }

#ifdef CPE_PROFILE
    lwpf_stop(simd_pf);

    lwpf_start(loopb_pf);
#endif
    int err_arr[SIMD_WIDTH];
    sindex_t vp_arr[SIMD_WIDTH];
    sindex_t vn_arr[SIMD_WIDTH];

    simd_store(errv, err_arr);
    simd_store(vpv, vp_arr);
    simd_store(vnv, vn_arr);

    for(iii = 0; iii < valid_size; iii ++){
        int last_high = band_length - read_length + ref_length[iii] - band_down - 1;
        int site=ref_length[iii]-last_high-1;
        int err = err_arr[iii];
        sindex_t vp = vp_arr[iii]; 
        sindex_t vn = vn_arr[iii];

        int location = -1;
        int error = MAXERROR;
        if((err <= dist_error)&&(err < error)){
            error = err;
            location = site;
        }

        for(i = 0; i < last_high; i++){
            err=err+((vp>>i)&(word_t)1);
            err=err-((vn>>i)&(word_t)1);
            if((err<=dist_error)&&(err<error)){
                error = err;
                location = site+i+1;
            }
        }

        err_loc_info[iii].error = error;
        err_loc_info[iii].location = location;

    }

#ifdef CPE_PROFILE
    lwpf_stop(loopb_pf);
#endif
}



sindex_t get_seed_from_reads(char *reads_s, sindex_t *reads_2bits, sindex_t *enc_2bits, sindex_t pre_seed, index_t seed_length, sindex_t seed_n, sindex_t read_length){


    index_t i = seed_length - 1;
    if(seed_n % SEED_STEP == 0){
        i = 0;  
    } 
    index_t reads_offset = (seed_n / SEED_STEP) * (seed_length + SEED_STEP - 1);
    index_t seed_offset = seed_n % SEED_STEP;
    sindex_t seed = pre_seed;
    index_t read_pos = reads_offset + seed_offset + i;
    for(;i < seed_length; i ++){
        char base = reads_s[reads_offset + seed_offset + i];
        read_pos = reads_offset + seed_offset + i;
        switch(base){
            case 'A':
                seed = (seed >> 2) | (0UL << ((seed_length - 1) * 2));
                *enc_2bits |= 0UL << ((read_pos % 16) * 2); 
                //reads_2bits[READ_INT_LEN + read_pos/16] |=  0UL << (read_pos%16);               
                //reads_2bits[READ_INT_LEN + read_pos/16] |=  0UL << (read_pos%16 + 16);               
                reads_2bits[READ_INT_LEN + (read_pos>>5)*2] |=  0UL << (read_pos&0x1f);               
                reads_2bits[READ_INT_LEN + (read_pos>>5)*2 + 1] |=  0UL << (read_pos&0x1f);               
                break;
            case 'C':
                seed = (seed >> 2) | (1UL << ((seed_length - 1) * 2));
                *enc_2bits |= 1UL << ((read_pos % 16) * 2); 
                //reads_2bits[READ_INT_LEN + read_pos/16] |=  1UL << (read_pos%16);               
                //reads_2bits[READ_INT_LEN + read_pos/16] |=  0UL << (read_pos%16 + 16);               
                reads_2bits[READ_INT_LEN + (read_pos>>5)*2] |=  1UL << (read_pos&0x1f);               
                reads_2bits[READ_INT_LEN + (read_pos>>5)*2 + 1] |=  0UL << (read_pos&0x1);               
                break;
            case 'G':
                seed = (seed >> 2) | (2UL << ((seed_length - 1) * 2));
                *enc_2bits |= 2UL << ((read_pos % 16) * 2); 
                //reads_2bits[READ_INT_LEN + read_pos/16] |=  0UL << (read_pos%16);               
                //reads_2bits[READ_INT_LEN + read_pos/16] |=  1UL << (read_pos%16 + 16);               
                reads_2bits[READ_INT_LEN + (read_pos>>5)*2] |=  0UL << (read_pos&0x1f);               
                reads_2bits[READ_INT_LEN + (read_pos>>5)*2 + 1] |=  1UL << (read_pos&0x1f);               
                break;
            case 'T':
                seed = (seed >> 2) | (3UL << ((seed_length - 1) * 2));
                *enc_2bits |= 3UL << ((read_pos % 16) * 2); 
                //reads_2bits[READ_INT_LEN + read_pos/16] |=  1UL << (read_pos%16);               
                //reads_2bits[READ_INT_LEN + read_pos/16] |=  1UL << (read_pos%16 + 16);               
                reads_2bits[READ_INT_LEN + (read_pos>>5)*2] |=  1UL << (read_pos&0x1f);               
                reads_2bits[READ_INT_LEN + (read_pos>>5)*2 + 1] |=  1UL << (read_pos&0x1f);               
                break;
            default:
                seed = (seed >> 2) | ((i%4) << ((seed_length - 1) * 2));
                *enc_2bits |= 1UL << ((read_pos % 16) * 2); 
                //reads_2bits[READ_INT_LEN + read_pos/16] |=  1UL << (read_pos%16);               
                //reads_2bits[READ_INT_LEN + read_pos/16] |=  0UL << (read_pos%16 + 16);               
                reads_2bits[READ_INT_LEN + (read_pos>>5)*2] |=  1UL << (read_pos&0x1f);               
                reads_2bits[READ_INT_LEN + (read_pos>>5)*2 + 1] |=  0UL << (read_pos&0x1f);               
                break;
        }
        if(read_pos % 16 == 15){
            reads_2bits[read_pos/16] = *enc_2bits;
            *enc_2bits = 0;
        }


        //if(_MYID == PRINT_CORE){
        //  printf("read_pos %llu\n", read_pos);
        //}
    }
    if(read_length - read_pos < seed_length){
        for(;read_pos < read_length; read_pos ++){
            char base = reads_s[read_pos];
            switch(base){
                case 'A':
                    *enc_2bits |= 0UL << ((read_pos % 16) * 2); 
                    //reads_2bits[READ_INT_LEN + read_pos/16] |=  0UL << (read_pos%16);               
                    //reads_2bits[READ_INT_LEN + read_pos/16] |=  0UL << (read_pos%16 + 16);               
                    reads_2bits[READ_INT_LEN + (read_pos>>5)*2] |=  0UL << (read_pos&0x1f);               
                    reads_2bits[READ_INT_LEN + (read_pos>>5)*2 + 1] |=  0UL << (read_pos&0x1f);               
                    break;
                case 'C':
                    *enc_2bits |= 1UL << ((read_pos % 16) * 2); 
                    //reads_2bits[READ_INT_LEN + read_pos/16] |=  1UL << (read_pos%16);               
                    //reads_2bits[READ_INT_LEN + read_pos/16] |=  0UL << (read_pos%16 + 16);               
                    reads_2bits[READ_INT_LEN + (read_pos>>5)*2] |=  1UL << (read_pos&0x1f);
                    reads_2bits[READ_INT_LEN + (read_pos>>5)*2 + 1] |=  0UL << (read_pos&0x1f);               
                    break;
                case 'G':
                    *enc_2bits |= 2UL << ((read_pos % 16) * 2); 
                    //reads_2bits[READ_INT_LEN + read_pos/16] |=  0UL << (read_pos%16);               
                    //reads_2bits[READ_INT_LEN + read_pos/16] |=  1UL << (read_pos%16 + 16);               
                    reads_2bits[READ_INT_LEN + (read_pos>>5)*2] |=  0UL << (read_pos&0x1f);               
                    reads_2bits[READ_INT_LEN + (read_pos>>5)*2 + 1] |=  1UL << (read_pos&0x1f);               
                    break;
                case 'T':
                    *enc_2bits |= 3UL << ((read_pos % 16) * 2); 
                    //reads_2bits[READ_INT_LEN + read_pos/16] |=  1UL << (read_pos%16);               
                    //reads_2bits[READ_INT_LEN + read_pos/16] |=  1UL << (read_pos%16 + 16);               
                    reads_2bits[READ_INT_LEN + (read_pos>>5)*2] |=  1UL << (read_pos&0x1f);               
                    reads_2bits[READ_INT_LEN + (read_pos>>5)*2 + 1] |=  1UL << (read_pos&0x1f);               
                    break;
                default:
                    *enc_2bits |= 1UL << ((read_pos % 16) * 2); 
                    //reads_2bits[READ_INT_LEN + read_pos/16] |=  1UL << (read_pos%16);               
                    //reads_2bits[READ_INT_LEN + read_pos/16] |=  0UL << (read_pos%16 + 16);               
                    reads_2bits[READ_INT_LEN + (read_pos>>5)*2] |=  1UL << (read_pos&0x1f);               
                    reads_2bits[READ_INT_LEN + (read_pos>>5)*2 + 1] |=  0UL << (read_pos&0x1f);               
                    break;
            }
            if(read_pos % 16 == 15){
                reads_2bits[read_pos/16] = *enc_2bits;
                *enc_2bits = 0;
            }
            //if(core_id == PRINT_CORE){
            //  printf("read_pos %llu\n", read_pos);
            //}
        }

        reads_2bits[(read_pos - 1)/16] = *enc_2bits;
    }
    return seed;
}
void choose_seed(char *reads_s, sindex_t read_length, sindex_t *hash_offset, 
        sindex_t *hash_count, hash_count_t *hash_seed_info, ref_info_t *ref_info_s, sindex_t *reads_2bits) {
    dma_declare();
    index_t j, k;
    index_t seeds_size = (read_length/(ref_info_s->seed_length + SEED_STEP - 1))*SEED_STEP;
    sindex_t seed  = 0;
    sindex_t enc_2bits = 0;
    for(j = 0; j < seeds_size; j ++){
        seed = get_seed_from_reads(reads_s, reads_2bits, &enc_2bits, seed, ref_info_s->seed_length, j, read_length);

        pe_get(ref_info_s->hash_offset + seed, hash_offset + j*2, 2*sizeof(sindex_t)); 
#ifndef DMA_ASYN
        dma_syn();
#endif
    }
#ifdef DMA_ASYN
    dma_syn();
#endif
    for(j = 0; j < seeds_size; j ++)
        hash_count[j] = hash_offset[j*2 + 1] - hash_offset[j*2];
    for(j = 0; j < 16; j ++){
        hash_seed_info[j].hash_count = 0;
        hash_seed_info[j].index = j;
    }
    for(j = 0; j <seeds_size; j ++)
        hash_seed_info[j/SEED_STEP].hash_count += hash_count[j]; 
    for(j = 0; j < seeds_size/SEED_STEP; j ++){
        for(k = j + 1; k < seeds_size/SEED_STEP; k ++){
            if(hash_seed_info[j].hash_count > hash_seed_info[k].hash_count){
                hash_count_t tmp = hash_seed_info[j];
                hash_seed_info[j] = hash_seed_info[k];
                hash_seed_info[k] = tmp;
            }
        }
    }

}

index_t adjust_ref_len(ref_info_t *ref_info_s, index_t *shift, sindex_t base_pos, sindex_t read_length){
    index_t dist_error = ref_info_s -> errors;
    index_t shift_tmp = ref_info_s -> errors; 
    index_t ref_len = read_length + 2 * shift_tmp;
    if(base_pos <= shift_tmp){
        shift_tmp = base_pos;
        ref_len -= (dist_error - base_pos);
    }
    int64_t front = (int64_t)(ref_info_s->ref_length) - (int64_t)(base_pos - dist_error);
    int64_t back = (int64_t)(base_pos + read_length + dist_error) - (int64_t)(ref_info_s->ref_length);
    if(front > back){
        if(back >= 0) ref_len -= (back + 1);

    }else{
        shift_tmp = dist_error - front - 1;
        ref_len -= front;
    }
    *shift = shift_tmp;
    return ref_len;

}
index_t pos_start(ref_info_t *ref_info_s, sindex_t base_pos, index_t *ref_shift, sindex_t *ref_length,sindex_t read_length){

    index_t tmp = adjust_ref_len(ref_info_s, ref_shift, base_pos, read_length); 
    *ref_length = (sindex_t)tmp;

    if(*ref_length < read_length - ref_info_s->errors){
        return  INVALID_SFLAG;
    }


    return base_pos - *ref_shift;
}

void put_cigar(ref_info_t *ref_info_s, cigar_t *cigar){
    index_t k, j;
    dma_declare();

#ifdef CPE_PROFILE
    lwpf_start(cigar_pf);
#endif
#ifdef OUT_CIGAR
    sindex_t cur_block = ref_info_s->cur_block;
    index_t cigar_ind;
    ref_info_s->cigar_size += cigar->cigar_num; 
    int unmap = 0;
    if(cigar->last_flag){
        if(ref_info_s->cigar_size == 0){
            unmap = 1;
        }
    }

    sindex_t cigar_num = cigar->cigar_num;
    if(unmap){
        asm volatile("faaw %0, 0(%1)\n\t":"=r"(cigar_ind):"r"((ref_info_s->unmap_ind)));
        pe_put(ref_info_s->cigar_num[cur_block] + cigar_ind,
                &(ref_info_s->cur_read) , 4);
    }
    if(cigar_num > 0){
        if(cigar->main_flag){
            cigar->cigar_num |= (1U << 31) ;
            cigar->main_flag = 0;
        }
        asm volatile("faaw %0, 0(%1)\n\t":"=r"(cigar_ind):"r"((ref_info_s->cigar_ind)));

        index_t cigar_offset =  cigar_ind * MAX_CIGAR * CIGAR_LEN; 
        sindex_t cur_block = ref_info_s->cur_block;

        //put results to HOST
        if(cigar->cigar_len > 0){
            //index_t coff =  cigar_offset + ref_info_s->cigar_len;
            if(cigar->cigar_len + 3 >  512 * 16){
                printf("wrong\n");
            }
            pe_put(ref_info_s->cigar_buffer[cur_block] + cigar_offset, 
                    cigar->cigar, ((cigar->cigar_len + 3)/4)*4);
        }

        sindex_t pos_off = cigar_ind * (MAX_CIGAR + 2);

        pe_put(ref_info_s->ref_pos_buffer[cur_block] + pos_off,  
                &(ref_info_s->cur_read), 
                sizeof(sindex_t));
        pe_put(ref_info_s->ref_pos_buffer[cur_block] + pos_off + 1,  
                &(cigar->cigar_num), 
                sizeof(sindex_t));
        pe_put(ref_info_s->ref_pos_buffer[cur_block] + pos_off + 2,  
                cigar->location, 
                cigar_num* sizeof(sindex_t));
        //cigar->location[-2] = ref_info_s->cur_read;
        //cigar->location[-1] = cigar->cigar_num;
        //pe_put(ref_info_s->ref_pos_buffer[cur_block] + pos_off,  
        //         cigar->location - 2, 
        //        (cigar->cigar_num  + 2)* sizeof(sindex_t));
    }
    dma_syn();

#endif
    cigar->cigar_num = 0;
    cigar->cigar_len = 0;

#ifdef CPE_PROFILE
    lwpf_stop(cigar_pf);
#endif
}

sindex_t shift_ref_pos_v(ref_info_t *ref_info_s, sindex_t *ref_seq, index_t *ref_pos, 
        sindex_t *ref_seqv, sindex_t *ref_length,index_t *valid_num,
        index_t *ref_shift, read_info_t *read_info){
    //index_t ref_enc_len = READ_LEN/16;
    index_t i = 0;
    index_t j = 0;
    index_t k = 0;
    //int seed_num = ref_info_s->errors + 1;
    int seed_len = ref_info_s->seed_length + SEED_STEP - 1;
    int seed_num = read_info->read_length / seed_len; 

    //sindex_t hmask = 0xffff0000;
    //sindex_t lmask = 0x0000ffff;
    sindex_t seed_mask = (1U << seed_len) - 1;
    index_t ref_ind = valid_num[0];
    for(i = valid_num[0]; i < valid_num[1]; i ++){
        index_t base_pos = ref_pos[i];
        index_t eq_flag = 0;
        index_t eq_flag1 = 0;
        sindex_t h_ref, l_ref;
        sindex_t read_enc_len = (ref_length[i] + 31) >> 5;
        //#undef LONG_SEED
#ifdef LONG_SEED        
        for(k = 0; k < seed_num; k ++){

            sindex_t long_seed = seed_len * k;
            sindex_t ref_off = (base_pos & 0x1f) + ref_shift[i] + long_seed;
            sindex_t shift_n = ref_off & 0x1f;
            sindex_t j = ref_off >> 5;

            l_ref =  (ref_seq[i*READ_INT_LEN + j * 2] >> shift_n) | ( ref_seq[i * READ_INT_LEN + (j + 1) * 2]  << (32 - shift_n));
            l_ref = l_ref & seed_mask;

            h_ref =  (ref_seq[i*READ_INT_LEN + j * 2 + 1] >> shift_n) | ( ref_seq[i * READ_INT_LEN + (j + 1) * 2 + 1]  << (32 - shift_n));
            h_ref = h_ref & seed_mask;



            sindex_t l_read = read_info->read_seeds[k*2];
            sindex_t h_read = read_info->read_seeds[k*2 + 1];

            sindex_t heq = h_ref ^ h_read;
            sindex_t leq = l_ref ^ l_read;
            sindex_t eq = (heq | leq);
            if(eq == 0){
                eq_flag = 1;
                break;
            }
        }

        if(eq_flag){       
#endif
            sindex_t shift_n = base_pos&0x1f;
            for(j = 0; j < read_enc_len; j ++){


                l_ref =  (ref_seq[i*READ_INT_LEN + j * 2] >> shift_n) | ( ref_seq[i * READ_INT_LEN + (j + 1) * 2]  << (32 - shift_n));
                h_ref =  (ref_seq[i*READ_INT_LEN + j * 2 + 1] >> shift_n) | ( ref_seq[i * READ_INT_LEN + (j + 1) * 2 + 1]  << (32 - shift_n));

                ref_seqv[j*2 * SIMD_WIDTH + ref_ind] = l_ref; 
                ref_seqv[(j*2 + 1) * SIMD_WIDTH + ref_ind] = h_ref; 

            }
            if(i != ref_ind){
                ref_pos[ref_ind] = ref_pos[i];
                ref_pos[i] = INVALID_SFLAG;
            }
            ref_ind ++;

#ifdef LONG_SEED
        }
#endif
    }
    valid_num[0] = ref_ind;
    valid_num[1] = ref_ind;


}

//#undef FIND
//#define FIND 2995776183 
void banded_myers(ref_info_t *ref_info_s, cigar_t *cigar, read_info_t *read_info, 
        sindex_t *positions, sindex_t pos_size){
    char *reads_s = read_info->read; 
    sindex_t *reads_2bits = read_info->read_2bits;
    sindex_t read_length = read_info->read_length;
    dma_declare();

    volatile int reply_ref = 0; 
    int count_ref = 0;
    //long ref_get_size = ((read_length + 2 * ref_info_s->errors + 15)/16 + 1) * sizeof(int);
    long ref_get_size = ((read_length + 2 * ref_info_s->errors + 31)/32 + 1) * 2 * sizeof(int);

    dma_desc pe_get_ref_desc = 0; 
    dma_set_op(&pe_get_ref_desc, DMA_GET);
    dma_set_reply(&pe_get_ref_desc, &reply_ref); 
    dma_set_size(&pe_get_ref_desc, ref_get_size); 

    index_t simd_width = SIMD_WIDTH;

    index_t i = 0; 
    index_t k, j;
    uintv8 seqv_bufferv[2*(READ_LEN/16)];
    index_t step_index = 0;
    sindex_t ref_seq[2][(READ_LEN/16) * simd_width];
    sindex_t *ref_seqv[2]; //[(READ_LEN/16) * simd_width] __attribute__((aligned(64)));
    ref_seqv[0] = (sindex_t*)&seqv_bufferv[0]; 
    ref_seqv[1] = (sindex_t*)&seqv_bufferv[(READ_LEN/16)]; 

    sindex_t ref_length[2][simd_width];

    index_t ref_pos[2][simd_width];
    index_t ref_shift[2][simd_width];

    error_t err_loc_info[simd_width];

    index_t ref_block = 0;
    index_t ref_enc_len = READ_LEN/16;
    //index_t read_enc_len =  (read_length + 15) / 16;

    index_t valid_size = 0; 
    index_t pre_valid_size = 0; 
    index_t valid_num[2][2];
    valid_num[0][0] = 0;
    valid_num[0][1] = 0;
    valid_num[1][0] = 0;
    valid_num[1][1] = 0;

    //sindex_t pre_pos = 0;
    //sindex_t base_pos = positions[0];

    //printf("pos size %u\n", pos_size);
    sindex_t count_pos = 0;

    sindex_t cur_block = ref_info_s->cur_block;

    index_t cigar_offset = ref_info_s->cur_read * ref_info_s->cigars_per_read
        * ref_info_s->cigar_avg_len; 

    index_t cur_cigar = ref_info_s->cur_read * ref_info_s->cigars_per_read;

#ifdef CPE_PROFILE
    lwpf_start(load_refc_pf);
#endif
    while(valid_size != simd_width && i != pos_size){ 
        //if(pre_pos != positions[i]){
        sindex_t base_pos = pos_start(ref_info_s, positions[i], &ref_shift[ref_block][valid_size], &ref_length[ref_block][valid_size], read_length);
        //if(positions[i] == FIND){
        //    printf("find pos5 %u %u\n", positions[i], i);
        //}
        //count_pos ++;
        //}
        //else
        //    base_pos = INVALID_SFLAG;
        //pre_pos = positions[i];
        i ++;
        if(base_pos == INVALID_SFLAG)
            continue;
        //index_t get_len = (base_pos + ref_length[ref_block][valid_size] + 15)/16 - base_pos / 16  + 1;
        //index_t get_len = ((base_pos + ref_length[ref_block][valid_size] + 31)/32 - base_pos / 32  + 1)*2;
        //int core_id;
        //GET_MYID(core_id);

        //if(core_id == PRINT_CORE){
        //    printf("count pos size core_id %lu %u %d %lu %lu\n", valid_size, get_len, core_id, base_pos, ref_info_s->base_2bits);

        //}
        //get_len = (((base_pos + ref_length[ref_block][valid_size] + 31)>>5) - (base_pos >> 5)  + 1) << 1;
        //pe_get(ref_info_s->base_2bits + base_pos/16, &(ref_seq[ref_block][valid_size * ref_enc_len]), get_len * sizeof(sindex_t)); 
        dma_rpl(pe_get_ref_desc, ref_info_s->base_2bits + (base_pos/32)*2, &(ref_seq[ref_block][valid_size * ref_enc_len]), reply_ref); 
        count_ref ++;

        ref_pos[ref_block][valid_size] = base_pos;
        //++++++++
#ifndef DMA_ASYN
        while (reply_ref != count_ref) { 
        };
        asm volatile("memb\n\t");                   
#endif

        valid_size ++;
    }
    valid_num[ref_block][0] = 0;
    valid_num[ref_block][1] = valid_size;
#ifdef DMA_ASYN
    while (reply_ref != count_ref) { 
    };
    asm volatile("memb\n\t");                   
#endif

    //    dma_syn();

#ifdef CPE_PROFILE
    lwpf_stop(load_refc_pf);
#endif
    while(i != pos_size){

        index_t pre_ref_block = ref_block;
        ref_block = (ref_block + 1) % 2;
        //pre_valid_size  = valid_size;
        valid_size = valid_num[ref_block][0];
#ifdef CPE_PROFILE
        lwpf_start(load_refc_pf);
#endif
        while(valid_size != simd_width && i != pos_size){ 

            sindex_t base_pos = pos_start(ref_info_s, positions[i], &ref_shift[ref_block][valid_size], &ref_length[ref_block][valid_size], read_length);
            i ++;
            if(base_pos == INVALID_SFLAG)
                continue;

            // index_t get_len = (base_pos + ref_length[ref_block][valid_size] + 15)/16 - base_pos / 16  + 1;
            // index_t get_len = ((base_pos + ref_length[ref_block][valid_size] + 31)/32 - base_pos / 32  + 1)*2;
            //pe_get(ref_info_s->base_2bits + base_pos/16, &ref_seq[ref_block][valid_size * ref_enc_len], get_len * sizeof(sindex_t)); 
            // dma_rpl(pe_get_ref_desc, ref_info_s->base_2bits + base_pos/16, &(ref_seq[ref_block][valid_size * ref_enc_len]), reply_ref); 
            dma_rpl(pe_get_ref_desc, ref_info_s->base_2bits + (base_pos/32)*2, &(ref_seq[ref_block][valid_size * ref_enc_len]), reply_ref); 
            count_ref ++;
            ref_pos[ref_block][valid_size] = base_pos;
            //++++++++
            valid_size ++;
#ifndef DMA_ASYN
            while (reply_ref != count_ref) { 
            };
            asm volatile("memb\n\t");                   
#endif 
        }
        valid_num[ref_block][1] = valid_size;

#ifdef CPE_PROFILE
        lwpf_stop(load_refc_pf);

        lwpf_start(shift_pf); 
#endif

        //shift_ref_pos(ref_seq[pre_ref_block], ref_pos[pre_ref_block], ref_seqv, ref_length[pre_ref_block], pre_valid_size);
        shift_ref_pos_v(ref_info_s, ref_seq[pre_ref_block], ref_pos[pre_ref_block], ref_seqv[pre_ref_block], 
                ref_length[pre_ref_block], valid_num[pre_ref_block], 
                ref_shift[pre_ref_block], read_info);

#ifdef CPE_PROFILE
        lwpf_stop(shift_pf);
#endif
        pre_valid_size  = valid_num[pre_ref_block][1];

        if(pre_valid_size == SIMD_WIDTH || i == pos_size){

#ifdef CPE_PROFILE
            lwpf_start(myers_pf);
#endif
#if defined SIMD_BM_ASM
            banded_myers_simd_asmv2(ref_info_s, reads_2bits, ref_seqv[pre_ref_block], ref_length[pre_ref_block], read_length, 
                    ref_shift[pre_ref_block], ref_pos[pre_ref_block], err_loc_info, pre_valid_size);
#elif defined SIMD_BM
            banded_myers_simd_v(ref_info_s, reads_2bits, ref_seqv[pre_ref_block], ref_length[pre_ref_block], read_length, 
                     ref_shift[pre_ref_block], ref_pos[pre_ref_block], err_loc_info, pre_valid_size);



#else
            banded_myers_simd(ref_info_s, reads_2bits, ref_seqv[pre_ref_block], ref_length[pre_ref_block], read_length, 
                    ref_shift[pre_ref_block], ref_pos[pre_ref_block], err_loc_info, pre_valid_size);
#endif

#ifdef CPE_PROFILE
            lwpf_stop(myers_pf);
#endif


            //banded_myers_simd(ref_info_s, reads_2bits, ref_seqv, ref_length[pre_ref_block], read_length, 
            //        ref_shift[pre_ref_block], ref_pos[pre_ref_block], err_loc_info, pre_valid_size);
#ifdef CPE_PROFILE
            lwpf_start(align_pf);
#endif
            align_reads(ref_info_s, reads_2bits, ref_seqv[pre_ref_block], ref_length[pre_ref_block], 
                    read_length, ref_shift[pre_ref_block], ref_pos[pre_ref_block], cigar, err_loc_info, pre_valid_size);
#ifdef CPE_PROFILE
            lwpf_stop(align_pf);
#endif
            valid_num[pre_ref_block][0] = 0;
            valid_num[pre_ref_block][1] = 0;
        }  
        if(cigar->cigar_num + 8 > MAX_CIGAR || cigar->cigar_len + 16 * CIGAR_LEN > MAX_CIGAR * CIGAR_LEN){
            //lwpf_start(cigar_pf);
            put_cigar(ref_info_s, cigar);
            //lwpf_stop(cigar_pf);
        }
#ifdef CPE_PROFILE
        lwpf_start(load_refc_pf);
#endif

#ifdef DMA_ASYN
        while (reply_ref != count_ref) { 
        };
        asm volatile("memb\n\t");                   
#endif 
        //dma_syn();
#ifdef CPE_PROFILE
        lwpf_stop(load_refc_pf);
#endif
}





index_t pre_ref_block = ref_block;
#ifdef CPE_PROFILE
lwpf_start(shift_pf); 
#endif
//pre_valid_size  = valid_size;
//shift_ref_pos(ref_seq[pre_ref_block], ref_pos[pre_ref_block], ref_seqv, ref_length[pre_ref_block], pre_valid_size);
shift_ref_pos_v(ref_info_s, ref_seq[pre_ref_block], ref_pos[pre_ref_block], ref_seqv[pre_ref_block], 
        ref_length[pre_ref_block], valid_num[pre_ref_block], 
        ref_shift[pre_ref_block], read_info);
#ifdef CPE_PROFILE
lwpf_stop(shift_pf);
#endif
pre_valid_size  = valid_num[pre_ref_block][1];
#ifdef CPE_PROFILE
lwpf_start(myers_pf);
#endif


#if defined SIMD_BM_ASM
banded_myers_simd_asmv2(ref_info_s, reads_2bits, ref_seqv[pre_ref_block], ref_length[pre_ref_block], read_length, 
        ref_shift[pre_ref_block], ref_pos[pre_ref_block], err_loc_info, pre_valid_size);
#elif defined SIMD_BM
banded_myers_simd_v(ref_info_s, reads_2bits, ref_seqv[pre_ref_block], ref_length[pre_ref_block], read_length, 
ref_shift[pre_ref_block], ref_pos[pre_ref_block], err_loc_info, pre_valid_size);

#else
banded_myers_simd(ref_info_s, reads_2bits, ref_seqv[pre_ref_block], ref_length[pre_ref_block], read_length, 
        ref_shift[pre_ref_block], ref_pos[pre_ref_block], err_loc_info, pre_valid_size);
#endif
//banded_myers_simd_asmv(ref_info_s, reads_2bits, ref_seqv, ref_length[pre_ref_block], read_length, 
//        ref_shift[pre_ref_block], ref_pos[pre_ref_block], err_loc_info, pre_valid_size);
#ifdef CPE_PROFILE
lwpf_stop(myers_pf);
#endif

#ifdef CPE_PROFILE
lwpf_start(align_pf);
#endif
align_reads(ref_info_s, reads_2bits, ref_seqv[pre_ref_block], ref_length[pre_ref_block], 
        read_length, ref_shift[pre_ref_block], ref_pos[pre_ref_block], cigar, err_loc_info, pre_valid_size);
#ifdef CPE_PROFILE
lwpf_stop(align_pf);
#endif
//lwpf_start(cigar_pf);
if(cigar->cigar_num + 8 > MAX_CIGAR || cigar->cigar_len + 16 * CIGAR_LEN > MAX_CIGAR * CIGAR_LEN){

    put_cigar(ref_info_s, cigar);
}
//lwpf_stop(cigar_pf);
//put results to HOST
}
void adjust_ref_pos_rmdup_v(ref_info_t *ref_info_s, sindex_t *positions, hash_count_t* hash_seed_info, 
        sindex_t *pos_offset, index_t pre_i, sindex_t pos_off_ind){
    index_t num = (ref_info_s->errors + 1) * SEED_STEP;
    index_t i, j;
    if(pos_off_ind == num){
        index_t pos_off = 0;
        for(i = 0; i < num; i ++){
            index_t base_offset = hash_seed_info[i/SEED_STEP].index * (ref_info_s->seed_length + SEED_STEP - 1)
                + i%SEED_STEP; 
            for(j = pos_off; j < pos_offset[i]; j ++){
                positions[j] = positions[j] > base_offset ? positions[j] - base_offset : 0; 
            }
            pos_off = pos_offset[i];
        }
    }else if(pos_off_ind != 0){
        i = pre_i;
        index_t base_offset = hash_seed_info[i/SEED_STEP].index * (ref_info_s->seed_length + SEED_STEP - 1)
            + i%SEED_STEP; 
        index_t pos_off = 0;
        if(pos_off_ind > 1){
            pos_off = pos_offset[pos_off_ind - 2];
        }
        for(j = pos_off; j < pos_offset[pos_off_ind - 1]; j ++){
            positions[j] = positions[j] > base_offset ? positions[j] - base_offset : 0; 
        }
    }
}

void adjust_ref_pos_rmdup(ref_info_t *ref_info_s, sindex_t *positions, hash_count_t* hash_seed_info, 
        sindex_t *pos_offset, index_t pre_i, sindex_t ajr_flag){
    index_t i, j;
    if(ajr_flag == 0){
        index_t pos_off = 0;
        index_t num = (ref_info_s->errors + 1) * SEED_STEP;
        for(i = 0; i < num; i ++){

            index_t base_offset = hash_seed_info[i/SEED_STEP].index * (ref_info_s->seed_length + SEED_STEP - 1)
                + i%SEED_STEP; 
            for(j = pos_off; j < pos_offset[i]; j ++){
               positions[j] = positions[j] > base_offset ? positions[j] - base_offset : 0; 
           }
            pos_off = pos_offset[i];
        }
    }else if(ajr_flag = 1){
        i = pre_i;
        index_t base_offset = hash_seed_info[i/SEED_STEP].index * (ref_info_s->seed_length + SEED_STEP - 1)
            + i%SEED_STEP; 
        for(j = 0; j < pos_offset[0]; j ++){

           positions[j] = positions[j] > base_offset ? positions[j] - base_offset : 0; 
       }
    }
}

void seeds_heap_insert(seed_ind_t *seeds_heap, sindex_t *seeds_heap_size,  sindex_t pos, sindex_t seed_ind){
    *seeds_heap_size += 1;
    seeds_heap[*seeds_heap_size].loc = pos;
    seeds_heap[*seeds_heap_size].seed_ind = seed_ind;
    sindex_t i = *seeds_heap_size;
    while(i > 1){
        sindex_t ti = i / 2;
        if(seeds_heap[ti].loc > seeds_heap[i].loc){
            seed_ind_t tmp = seeds_heap[ti];
            seeds_heap[ti] = seeds_heap[i];
            seeds_heap[i] = tmp;
        }else{
            break;
        }
        i =  ti;
    }
}
void seeds_heap_pop_and_insert(seed_ind_t *seeds_heap, sindex_t *seeds_heap_size,  sindex_t pos, sindex_t del_flag){
    if(del_flag == 1){
        //if(seeds_heap_size > 0)
        seeds_heap[1] = seeds_heap[*seeds_heap_size]; 
        *seeds_heap_size -= 1;
    }else{
        seeds_heap[1].loc = pos;
    }
    sindex_t i = 1;
    sindex_t ti = 0;
    while(i * 2 < *seeds_heap_size){
        ti = i * 2;
        if(i * 2 + 1 <= *seeds_heap_size){
            if(seeds_heap[i*2].loc > seeds_heap[i*2 + 1].loc){
                ti = i * 2 + 1; 
            }
        }
        if(seeds_heap[ti].loc < seeds_heap[i].loc){
            seed_ind_t tmp = seeds_heap[ti];
            seeds_heap[ti] = seeds_heap[i];
            seeds_heap[i] = tmp;
        }else{
            break;
        }
        i =  ti;
    }
}

void qsort_v(sindex_t *v, int left, int right){
    int i, last;
    sindex_t  tmp;
    if (left >= right)
        return;
    //swap(v, left, (left + right)/2);
    tmp = v[left];
    v[left] = v[(left + right)/2];
    v[(left + right)/2] = tmp;


    last = left;
    for (i = left+1; i <= right; i++){
        if (v[i] < v[left]){
            tmp = v[++last];
            v[last] = v[i];
            v[i] = tmp;

        }
    }

    tmp = v[last];
    v[last] = v[left];
    v[left] = tmp;


    qsort_v(v, left, last-1);
    qsort_v(v, last+1, right);

}

typedef struct pos_heap_c{
    sindex_t pos;
    sindex_t pos_ind;
}pos_heap_t;
sindex_t block_heap_sort(sindex_t *positions, pos_heap_t *pos_heap, sindex_t *pos_ind, sindex_t block_size){
    sindex_t i, j, k;
    sindex_t heap_size = 0;
    sindex_t invalid_value = 0xffffffff;
    for( i = 0; i < block_size/2; i ++){
        if(pos_ind[i*2] < pos_ind[2*i+ 1]){
            pos_heap[++heap_size].pos = positions[pos_ind[2*i]]; 
            pos_heap[heap_size].pos_ind = 2*i;
            pos_ind[2*i] ++;

            index_t cur_ind = heap_size;
            index_t p_ind  = cur_ind/2;
            while(cur_ind > 1){
                if(pos_heap[p_ind].pos > pos_heap[cur_ind].pos){
                    sindex_t tpos = pos_heap[cur_ind].pos;
                    sindex_t tind = pos_heap[cur_ind].pos_ind;
                    pos_heap[cur_ind].pos = pos_heap[p_ind].pos;
                    pos_heap[cur_ind].pos_ind = pos_heap[p_ind].pos_ind;
                    pos_heap[p_ind].pos = tpos;
                    pos_heap[p_ind].pos_ind = tind;
                }else{
                    break;
                }
                cur_ind /= 2;
                p_ind = cur_ind / 2;
            }
        }
    }
    sindex_t pre_small = invalid_value; 
    //+++++
    sindex_t cur_ind = 0;
    //--------
    while(heap_size > 1){
        sindex_t pos_id =  pos_heap[1].pos_ind;
        if(pre_small != pos_heap[1].pos){
            positions[2048 + cur_ind ++] = pos_heap[1].pos;
            pre_small = pos_heap[1].pos;
        }
        if(pos_ind[pos_id] < pos_ind[pos_id + 1]){
            pos_heap[1].pos = positions[pos_ind[pos_id]];
            pos_ind[pos_id] ++;
        } else if(pos_ind[pos_id] == pos_ind[pos_id]){
            pos_heap[1].pos = pos_heap[heap_size].pos;
            pos_heap[1].pos_ind = pos_heap[heap_size].pos_ind;
            heap_size --;
            if(heap_size == 1){
                break;
            }
        }

        sindex_t heap_ind = 1;
        while(heap_ind * 2 <= heap_size){
            sindex_t l_ind = heap_ind * 2;
            sindex_t r_ind = l_ind  + 1;
            sindex_t s_ind = l_ind;
            if(r_ind <= heap_size){
                if(pos_heap[l_ind].pos > pos_heap[r_ind].pos){
                    s_ind  = r_ind; 
                }
            }
            if(pos_heap[heap_ind].pos > pos_heap[s_ind].pos){
                sindex_t tpos = pos_heap[heap_ind].pos;
                sindex_t tind = pos_heap[heap_ind].pos_ind;
                pos_heap[heap_ind].pos = pos_heap[s_ind].pos;
                pos_heap[heap_ind].pos_ind = pos_heap[s_ind].pos_ind;
                pos_heap[s_ind].pos = tpos;
                pos_heap[s_ind].pos_ind = tind;

            }else{
                break;
            }
            heap_ind = s_ind;
        }
    }

    sindex_t pos_id =  pos_heap[1].pos_ind;
    for(i = pos_ind[pos_id] - 1; i < pos_ind[pos_id + 1]; i ++){
        sindex_t pos = positions[i]; 
        if(pos != pre_small){
            positions[2048 + cur_ind ++] = pos;
            pre_small = pos;
        }
    }

    //sindex_t block_ind = 1;
    pos_ind[0] = cur_ind;
    return 1;

    //----------------

}
void long_seeds_filter(ref_info_t *ref_info_s, read_info_t *read_info, sindex_t *positions, half_word_t *pos_lseed, sindex_t *pos_ind, 
        half_word_t *pos_lseed_ind, half_word_t *pos_lseed_off, sindex_t block_size){

    assert(SEED_STEP <= 5);
    assert(SEED_STEP >= 1);
    sindex_t res_len = SEED_STEP - 1;
    sindex_t mask = (1 << res_len) - 1;
    sindex_t i, j, k;
    sindex_t pos_size = 0;
    for( i = 0; i < block_size/2; i ++){
        sindex_t lind  = pos_lseed_ind[i];
        sindex_t shift = SEED_STEP - 1 - pos_lseed_ind[i] % SEED_STEP;
        sindex_t lseed_off = pos_lseed_off[i];

        sindex_t rl_res_seed = read_info->res_seeds[2*lind];
        sindex_t rh_res_seed = read_info->res_seeds[2*lind + 1];
        for(j = pos_ind[i*2]; j < pos_ind[2*i+ 1]; j ++){
            sindex_t sl_res_seed = pos_lseed[lseed_off] >> shift; 
            sindex_t sh_res_seed = pos_lseed[lseed_off] >> (res_len * 2 + shift); 


            sindex_t result = (sl_res_seed ^ rl_res_seed) | (sh_res_seed ^ rh_res_seed);
            result &= mask;
            if(result == 0 && pos_size != j){
                positions[pos_size ++] = positions[j];
            }
            lseed_off ++;
        }
        if(i > 0)
            pos_ind[i*2] = pos_ind[i*2 - 1];
        pos_ind[i*2 + 1] = pos_size;
    }

}
void alignment_rmdup(ref_info_t *ref_info_s, read_info_t *read_info, hash_count_t *hash_seed_info,
        sindex_t *hash_offset, sindex_t *hash_count){
    //qsort+++++++++
    //index_t positions_size = 4096;
    //---------------
    index_t positions_size = 2048;
    sindex_t positions[4096];
    sindex_t pos_ind[264];
    //#undef LSEED_FILTER
#ifdef LSEED_FILTER
    half_word_t pos_lseed_ind[128 + 8];
    half_word_t pos_lseed_off[128 + 8];
    half_word_t *pos_lseeds = (half_word_t*)&positions[positions_size]; 
#endif
    pos_heap_t pos_heap[264];
    sindex_t block_size = 0;

    char *reads_s = read_info->read; 
    sindex_t *reads_2bits = read_info->read_2bits; 
    sindex_t read_length = read_info->read_length;

    // every type of seeds(seed_length) steps
    // [step:0~SEE_STEP, pos_size: offset of positions for this seeds, index:index of seeds(seed_length) 
    //hash_steps_t hash_steps[64];

    seed_ind_t seeds_heap[64]; 
    sindex_t seeds_heap_size = 0;

    index_t pos_count[64];

    sindex_t pos_offset[64];
    sindex_t pos_off_ind = 0;

    //sindex_t pos_offset = 0;

    //uint32_t cigar_int[15*MAX_CIGAR/4];
    uint32_t cigar_int[CIGAR_LEN*MAX_CIGAR/4];
    sindex_t location[512 + 8];
    cigar_t cigar;
    cigar.cigar = (char*)cigar_int;
    cigar.location = location + 2;
    cigar.cigar_len = 0;
    cigar.cigar_num = 0;
    cigar.main_flag = 1;
    cigar.last_flag = 0;
    cigar.f_flag = 0;


    dma_declare();
    index_t i, j, k;
    i = 0;
    j = 0;

    index_t seeds_num = (ref_info_s->errors + 1) * SEED_STEP; 
    index_t seeds_per_block = 32;

    block_size = 1;
    pos_ind[0] = 0;
    pos_ind[1] = 0;

    sindex_t pos_size = 0;
#ifdef LSEED_FILTER
    sindex_t lseed_size = 0;
    sindex_t lseed_off = 0;
    sindex_t lseed_get = 0;
#endif
#ifdef CPE_PROFILE
    lwpf_start(pnum_pf);
#endif

    for( i = 0;  i < seeds_num; i ++){
        index_t seed_index = hash_seed_info[i/SEED_STEP].index * SEED_STEP + i % SEED_STEP;
        index_t seeds_count = hash_count[seed_index];
        index_t pos_get = 0;
        if(seeds_count > seeds_per_block){
            pos_get = seeds_per_block;
        }else{
            pos_get = seeds_count;
        }
        if(pos_get > 0){
            pe_get(ref_info_s->hash_positions + hash_offset[seed_index * 2], positions + pos_size, pos_get * sizeof(sindex_t));

#ifdef LSEED_FILTER
            lseed_off = hash_offset[seed_index * 2] - (hash_offset[seed_index * 2] & 1U);
            lseed_get = pos_get +  (hash_offset[seed_index * 2] & 1U);
            lseed_get = lseed_get + (lseed_get & 1U);
            pe_get(ref_info_s->pos_lseeds + lseed_off, pos_lseeds + lseed_size, lseed_get * sizeof(half_word_t));
            //lseed_off = (hash_offset[seed_index * 2] & 1);
#endif
        }
        pos_size += pos_get;


        pos_count[i] = pos_get;
        pos_offset[i] = pos_size;

        pos_off_ind ++;
        if(pos_get > 0){
#ifdef LSEED_FILTER
            pos_lseed_off[block_size/2] = lseed_size + (hash_offset[seed_index * 2] & 1);
            pos_lseed_ind[block_size/2] = seed_index;
            lseed_size += lseed_get;
#endif
            pos_ind[block_size++] = pos_size;
            pos_ind[block_size++] = pos_size;
        }
        //pos_get[0] = 0;
        //if(_MYID == PRINT_CORE)
        //  printf("seeds count pos_get %lu %lu %lu\n", seeds_count, pos_get, i);
#ifndef  DMA_ASYN 
        dma_syn();
#endif
    }
#ifdef DMA_ASYN
    dma_syn();
#endif
#ifdef CPE_PROFILE
    lwpf_stop(pnum_pf);
    lwpf_start(heap_pf);
#endif
    for( i = 0; i < seeds_num; i ++){
        index_t seed_index = hash_seed_info[i/SEED_STEP].index * SEED_STEP + i % SEED_STEP;
        index_t seeds_count = hash_count[seed_index];
        if(pos_count[i] != seeds_count){
            seeds_heap_insert(seeds_heap, &seeds_heap_size,  positions[pos_offset[i] - 1], i);
            //if(_MYID == PRINT_CORE)
            //    printf("pos: %u %lu\n", positions[pos_offset[i]-1], i);
        }
    }

#ifdef CPE_PROFILE
    lwpf_stop(heap_pf);
#endif
    index_t pre_pos_size = 0;
    index_t ajr_flag = 0;
    index_t pre_i = 0;

    while(1){
        //qsort +++++
        //block_size = 0;
        //----------
        if(seeds_heap_size > 0 && pos_size + seeds_per_block <= positions_size && block_size < 256){
            i = seeds_heap[1].seed_ind;
            index_t seed_index = hash_seed_info[i/SEED_STEP].index * SEED_STEP + i % SEED_STEP;
            index_t seeds_count = hash_count[seed_index];
            index_t pos_get = 0;

            if(seeds_count > pos_count[i] + seeds_per_block){
                pos_get = seeds_per_block;
            }else{
                pos_get = seeds_count - pos_count[i];
            }


#ifdef CPE_PROFILE
            lwpf_start(pnum_pf);
#endif
            if(pos_get > 0){

                sindex_t pos_off = hash_offset[seed_index * 2] + pos_count[i]; 
                pe_get(ref_info_s->hash_positions + pos_off, positions + pos_size, pos_get * sizeof(sindex_t));

#ifdef LSEED_FILTER
                lseed_off = pos_off - (pos_off&1U); 
                lseed_get = pos_get +  ( pos_off & 1U);
                lseed_get = lseed_get + (lseed_get & 1U);
                pe_get(ref_info_s->pos_lseeds + lseed_off , pos_lseeds + lseed_size, lseed_get * sizeof(half_word_t));
                lseed_off = pos_off & 1U; 
#endif

            }

#ifndef DMA_ASYN
            dma_syn();
#endif
            if(pos_size > 0){
                adjust_ref_pos_rmdup(ref_info_s, positions + pre_pos_size, hash_seed_info, pos_offset, pre_i, ajr_flag);
            }


            //adjust_ref_pos_rmdup(ref_info_s, positions + pre_pos_size, hash_seed_info, pos_offset, pre_i, pos_off_ind);
            //adjust_ref_pos_rmdup(ref_info_s, positions, hash_seed_info, pos_offset, pre_i, pos_off_ind);

            // for adjust_ref_pos_rmdup
            pre_pos_size = pos_size;
            pos_offset[0] = pos_get;
            //pos_offset[pos_off_ind ++] = pos_get;
            //pos_offset[pos_off_ind ++] = pos_size + pos_get;
            pre_i = i;
            ajr_flag = 1;

            pos_size += pos_get;
            pos_count[i] += pos_get;

            if(pos_get > 0){
#ifdef LSEED_FILTER
                pos_lseed_off[block_size/2] = lseed_size + lseed_off;
                pos_lseed_ind[block_size/2] = seed_index;
                lseed_size += lseed_get;
#endif

                pos_ind[block_size++] = pos_size;
                pos_ind[block_size++] = pos_size;

                //pos_ind[block_size++] = pos_size;
                //pos_ind[block_size++] = pos_size;
            }


#ifdef DMA_ASYN
            dma_syn();
#endif

#ifdef CPE_PROFILE
            lwpf_stop(pnum_pf);
#endif
            sindex_t del_flag = 0;
            if(pos_count[i] == seeds_count){
                del_flag = 1; 
            }else{
                del_flag = 0;
            }

#ifdef CPE_PROFILE
            lwpf_start(heap_pf);
#endif
            seeds_heap_pop_and_insert(seeds_heap, &seeds_heap_size,  positions[pos_size - 1], del_flag);
#ifdef CPE_PROFILE
            lwpf_stop(heap_pf);
#endif

        }
        if(pos_size + seeds_per_block > positions_size || seeds_heap_size == 0 || block_size >= 256){
            adjust_ref_pos_rmdup(ref_info_s, positions + pre_pos_size, hash_seed_info, pos_offset, pre_i, ajr_flag);

            //adjust_ref_pos_rmdup(ref_info_s, positions + pre_pos_size, hash_seed_info, pos_offset, pre_i, pos_off_ind);
            //adjust_ref_pos_rmdup(ref_info_s, positions, hash_seed_info, pos_offset, pre_i, pos_off_ind);

#ifdef CPE_PROFILE
            lwpf_start(heap_sort_pf);
#endif
            if(pos_size > 0){
                // qsort_v(positions, 0, pos_size - 1);
                //printf("pos size %u\n", pos_size);
                //HeapSort(positions, pos_size);

#ifdef LSEED_FILTER
                long_seeds_filter(ref_info_s, read_info, positions, pos_lseeds, pos_ind, pos_lseed_ind, pos_lseed_off, block_size - 1);
                lseed_size = 0;
#endif 
                pos_size = pos_ind[block_size - 2];
#ifdef SORT_POS

                block_size = block_heap_sort(positions, pos_heap, pos_ind, block_size - 1);
                pos_size = pos_ind[0];

                //pos_size = pos_ind[block_size - 1];
#endif
                //printf("pos size %u\n", pos_size);
            }

#ifdef CPE_PROFILE
            lwpf_stop(heap_sort_pf);
#endif
           if(pos_size > 0){

#ifdef CPE_PROFILE
                lwpf_start(bmyers_pf);
#endif
                //banded_myers(ref_info_s, &cigar, read_info, positions, pos_size);
#ifdef SORT_POS
                banded_myers(ref_info_s, &cigar, read_info, positions + 2048, pos_size);
#else
                banded_myers(ref_info_s, &cigar, read_info, positions, pos_size);
#endif

#ifdef CPE_PROFILE
                lwpf_stop(bmyers_pf);
#endif
            }

            pos_size = 0;
            pos_offset[0] = 0;
            pre_pos_size = 0;
            ajr_flag = 2;

            block_size = 1;
            pos_ind[0] = 0;
            pos_ind[1] = 0;

            pos_off_ind = 0;

            if(seeds_heap_size == 0)
                break;
        }
    }

    if(ref_info_s->cur_read & 1)
        cigar.last_flag = 1;
    cigar.f_flag = 1;
    put_cigar(ref_info_s, &cigar);

}


void reverse_read(char *reverse_seq, char *read, sindex_t read_len){
    index_t i;
    index_t j = read_len - 1;
    for(i = 0; i < read_len; i ++){
        switch(read[i]){
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
                reverse_seq[j-i] = read[i];
        }
    }
}
void set_read_long_seeds(ref_info_t *ref_info, read_info_t* read_info){

    //sindex_t *read_2bits, sindex_t *read_seeds, sindex_t read_len){
    sindex_t *read_2bits = read_info->read_2bits;
    sindex_t *read_seeds = read_info->read_seeds;

    uint8_t *res_seeds = read_info->res_seeds;

    sindex_t read_len = read_info -> read_length;
    int i = 0;
    //int seed_num = ref_info->errors + 1;
    int seed_len = ref_info->seed_length + SEED_STEP - 1;
    //read_len / seed_num; 
    int seed_num = read_len / seed_len;
    sindex_t seed_mask = (1U << seed_len) - 1;
    sindex_t res_mask = (1U << (SEED_STEP - 1)) - 1;
    for(i = 0 ; i < seed_num; i ++){
        index_t seed_off = i * seed_len;
        index_t j = seed_off >> 5;
        index_t shift_n = seed_off & 0x1f;
        index_t h_read = 0;   
        index_t l_read = 0;

        l_read =  (read_2bits[READ_INT_LEN + j*2]  >>  shift_n ) | ( read_2bits[READ_INT_LEN + (j + 1)*2]   << ( 32 - shift_n));
        l_read = l_read & seed_mask;

        h_read =  (read_2bits[READ_INT_LEN + j*2 + 1]  >>  shift_n ) | ( read_2bits[READ_INT_LEN + (j + 1)*2 + 1]   << ( 32 - shift_n));
        h_read = h_read & seed_mask;


        read_seeds[i*2] = l_read;
        read_seeds[i*2 + 1] = h_read;
        for(j = 0; j < SEED_STEP; j ++){
            sindex_t shift = ref_info->seed_length + j;
            sindex_t mask = (1U << j) - 1;

            sindex_t fseed = l_read & mask; 
            sindex_t bseed = (l_read >> shift) << j; 
            res_seeds[2*(SEED_STEP*i + j)] = ( fseed | bseed) & res_mask;

            fseed = h_read & mask; 
            bseed = (h_read >> shift) <<  j; 

            res_seeds[2*(SEED_STEP*i + j) + 1] = ( fseed | bseed) & res_mask;

        }
    }

}
void get_positions(ref_info_t *ref_info_s, char *reads_s, sindex_t *reads_l, sindex_t reads_per_core){
    assert(READ_LEN <= 256);
    // every type of seeds(seed_length) info 
    sindex_t hash_offset[128];
    sindex_t hash_count[64];
    sindex_t reads_2bits[32];
    //sindex_t reads_2bits[16];
    sindex_t read_seeds[64];
    uint8_t res_seeds[128];
    read_info_t read_info;
    read_info.read_2bits = reads_2bits;
    read_info.read_seeds = read_seeds;
    read_info.res_seeds = res_seeds;

    // every type of long seeds(seed_length + SEED_STEP - 1) info 
    hash_count_t hash_seed_info[16];
    index_t i, j, k;
    //reads_per_core = 3;
    //index_t read_sn_length = ref_info_s->reads_name_length + ref_info_s->reads_seq_length; 
    index_t read_sn_length = ref_info_s->reads_seq_length; 
    char reverse_seq[READ_LEN];

    for( i = 0; i < reads_per_core; i ++){
        read_info.read_length = reads_l[i];
        for(j = 0; j < READ_INT_LEN; j ++){
            reads_2bits[READ_INT_LEN + j] = 0;
            reads_2bits[j] = 0;
        }
        ref_info_s->cigar_size = 0;
        ref_info_s->cigar_len = 0;

#ifdef CPE_PROFILE
        lwpf_start(choose_seed_pf);
#endif
        choose_seed(reads_s + i*read_sn_length, reads_l[i], hash_offset, hash_count, 
                hash_seed_info, ref_info_s, reads_2bits); 

#ifdef CPE_PROFILE
        lwpf_stop(choose_seed_pf);
#endif

#ifdef CPE_PROFILE
        lwpf_start(alignment_pf);
#endif
        set_read_long_seeds(ref_info_s, &read_info);
        read_info.read_length = reads_l[i];
        read_info.read = reads_s + i *read_sn_length;

        alignment_rmdup(ref_info_s, &read_info, 
                hash_seed_info, hash_offset, hash_count);

#ifdef CPE_PROFILE
        lwpf_stop(alignment_pf);
#endif
        ref_info_s->cur_read += 1; 

        for(j = 0; j < READ_INT_LEN; j ++){
            reads_2bits[READ_INT_LEN + j] = 0;
            reads_2bits[j] = 0;
        }

        ref_info_s->cigar_len = 0;
        reverse_read(reverse_seq, reads_s + i*read_sn_length, reads_l[i]);

#ifdef CPE_PROFILE
        lwpf_start(choose_seed_pf);
#endif
        choose_seed(reverse_seq, reads_l[i], hash_offset, hash_count, 
                hash_seed_info, ref_info_s, reads_2bits); 

#ifdef CPE_PROFILE
        lwpf_stop(choose_seed_pf);

        lwpf_start(alignment_pf);
#endif
        set_read_long_seeds(ref_info_s, &read_info);
        read_info.read_length = reads_l[i];
        read_info.read = reverse_seq;

        alignment_rmdup(ref_info_s, &read_info, 
                hash_seed_info, hash_offset, hash_count);

#ifdef CPE_PROFILE
        lwpf_stop(alignment_pf);
#endif

        ref_info_s->cur_read +=  1; 

    }

}

void find_positions(ref_info_t *ref_info_h){
#ifdef CPE_PROFILE
    long bget = bytecount_shadow[DMA_GET];
    long bput = bytecount_shadow[DMA_PUT];
    long bgeta = bytecount_shadow[DMA_GET + 8];
    long bputa = bytecount_shadow[DMA_PUT + 8];
    lwpf_enter(MAPPING);
#endif
    dma_declare();

    ref_info_t  ref_info_s;
    pe_get(ref_info_h, &ref_info_s, sizeof(ref_info_t));
    dma_syn();

    sindex_t reads_l[10]; 
    char reads_s[256]; 
    sindex_t reads_per_core = 1;
    sindex_t cur_read_offset = 0;
    sindex_t cur_reads_num = 0;


    //index_t reads_length = ref_info_s.reads_name_length + ref_info_s.reads_seq_length;
    index_t reads_length = ref_info_s.reads_seq_length;
    index_t cur_loop = 0;


    index_t cur_block = ref_info_s.cur_block;
#ifdef BALANCE
    asm volatile("faaw %0, 0(%1)\n\t":"=r"(cur_loop) : "r"(&(ref_info_h->lock_read)));
    //ref_info_s.cur_read = 2 * (cur_loop * 64  + core_id) * reads_per_core;
#else
    cur_loop = _MYID;
#endif
    while(cur_loop * reads_per_core < ref_info_s.reads_num){

        cur_read_offset = cur_loop * reads_per_core;
        ref_info_s.cur_read = 2 * cur_loop * reads_per_core;

        if((cur_loop + 1) * reads_per_core < ref_info_s.reads_num){
            cur_reads_num = reads_per_core;
        }else{
            cur_reads_num = ref_info_s.reads_num - cur_read_offset;
        }

        pe_get(ref_info_s.reads[cur_block] + cur_read_offset * reads_length, reads_s, cur_reads_num * reads_length); 
        pe_get(ref_info_s.reads_seq_len[cur_block] + cur_read_offset, reads_l, cur_reads_num * sizeof(sindex_t)); 
        dma_syn();
        get_positions(&ref_info_s, reads_s, reads_l, cur_reads_num);
#ifdef CPE_PROFILE
        lwpf_start(faaw_pf);
#endif
#ifdef BALANCE
        asm volatile("faaw %0, 0(%1)\n\t":"=r"(cur_loop):"r"(&(ref_info_h->lock_read)));
#else
        cur_loop += 64;
#endif
#ifdef CPE_PROFILE
        lwpf_stop(faaw_pf);
#endif
    }

#ifdef CPE_PROFILE
    bget -= bytecount_shadow[DMA_GET];
    bput -= bytecount_shadow[DMA_PUT];
    bgeta -= bytecount_shadow[DMA_GET + 8];
    bputa -= bytecount_shadow[DMA_PUT + 8];
    //cal_locked_printf("id: %d, get: %lld, put: %lld, geta: %lld, puta: %lld\n",
    //        _MYID, -bget, -bput, -bgeta, -bputa);

    lwpf_exit(MAPPING);
#endif
}
