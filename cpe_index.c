/*************************************************************************
	> File Name: build_index.c
	> Author: Xu Kai
	> Created Time: Sun 08 Dec 2019 02:02:38 PM CST
 ************************************************************************/
	
	
#include<stdio.h>
#include<stdlib.h>
#include<slave.h>
#include<assert.h>
#include "dma_macros.h"
//#include "/home/export/online1/swmore/opensource/cal/cal.h"
#include "fm_define.h"
//#define sindex_t uint32_t 
//#define index_t uint64_t
#define BITS_SIZE 32
//#define INVALID_SFLAG 0xffffffff
#define CORE_SIZE 64
//#define SEED_STEP 4
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
}build_info_t;
#define encode_base(c, href, lref, nref, ref_seeds, pos) { \
    seed_len ++; \
    switch(toupper(c)){ \
        case('A'): \
            lref |= 0U << (pos & 0x1f); \
            href |= 0U << (pos & 0x1f); \
            ref_seeds |= 0U << ((pos & 0xf) * 2); \
            seed = (seed >> 2) | (0U << (2 * seed_length - 2)); \
            break; \
        case('C'): \
            lref |= 1U << (pos & 0x1f); \
            href |= 0U << (pos & 0x1f); \
            ref_seeds |= 1U << ((pos & 0xf) * 2); \
            seed = (seed >> 2) | (1U << (2 * seed_length - 2)); \
            break; \
        case('G'): \
            lref |= 0U << (pos & 0x1f); \
            href |= 1U << (pos & 0x1f); \
            ref_seeds |= 2U << ((pos & 0xf) * 2); \
            seed = (seed >> 2) | (2U << (2 * seed_length - 2)); \
            break; \
        case('T'): \
            lref |= 1U << (pos & 0x1f); \
            href |= 1U << (pos & 0x1f); \
            ref_seeds |= 3U << ((pos & 0xf) * 2); \
            seed = (seed >> 2) | (3U << (2 * seed_length - 2)); \
            break; \
        default: \
            lref |= 0U << (pos & 0x1f); \
            href |= 0U << (pos & 0x1f); \
            nref |= 1U << (pos & 0x1f); \
            ref_seeds |= 0U << ((pos & 0xf) * 2); \
            seed_len = 0; \
            seed = 0;  \
            break; \
    } \
  }

void build_hash_index(build_info_t* build_info_h){
    build_info_t build_info;
    dma_declare();

    pe_get(build_info_h, &build_info, sizeof(build_info_t));
    dma_syn();

    index_t core_id = _MYID;

    sindex_t bits_num = (build_info.ref_seq_len + BITS_SIZE - 1) / BITS_SIZE; 
    sindex_t seq_off = core_id * ((bits_num  +  CORE_SIZE - 1)/ CORE_SIZE); 

    //sindex_t part_size = 0; 

    sindex_t part_size = ((bits_num  +  CORE_SIZE - 1) / CORE_SIZE) * BITS_SIZE;

    if(seq_off*BITS_SIZE + part_size > build_info.ref_seq_len)
        part_size = build_info.ref_seq_len - seq_off * BITS_SIZE;

    if(seq_off * BITS_SIZE > build_info.ref_seq_len)
        part_size = 0; 

    if(part_size == 0)
        return;
    
    sindex_t res_buffer[8];
    char *res_seq = (char*) res_buffer;

    //sindex_t part_num = 4096*4;
    sindex_t get_size = 2048 * 4;

    sindex_t seq_buffer[2048 + 8];
    char *ref_seq = (char*)seq_buffer;

    sindex_t base_2bits[1024 + 8];
    sindex_t base_3bits[1024 + 8];
    sindex_t seeds_2bits[1024 + 8];
    sindex_t seeds[2048 + 8];

    sindex_t href_2bits = 0; 
    sindex_t lref_2bits = 0;
    sindex_t nref_3bits = 0;
    sindex_t ref_seeds = 0;

    index_t i, k, j;

    sindex_t seed = 0;
    sindex_t seed_len = 0;

    sindex_t seed_length = 12;
    //sindex_t lseed_len = seed_length + SEED_STEP - 1;
    //sindex_t seed_mask = (1U << (seed_length * 2))  - 1;

    sindex_t ref_off  = seq_off * BITS_SIZE;
    sindex_t base_ref_off  = seq_off * BITS_SIZE;

    if(core_id == 0){
        pe_get(build_info.res_seq, res_seq, BITS_SIZE);
    }else{
        //index_t tmp = (build_info.ref_seq + ref_off);
        //    tmp = tmp&0x3U;
        //if(tmp != 0)
        //    printf("offload wrong \n");
        //if(part_size > 0)
        pe_get(build_info.ref_seq + ref_off - BITS_SIZE, res_seq, BITS_SIZE);
    }
    dma_syn();
    index_t start_loc = BITS_SIZE - build_info.res_len;
    sindex_t spos = 0;
    sindex_t lpos = 0;
    
    for(i = start_loc; i < BITS_SIZE; i ++){
        encode_base(res_seq[i], href_2bits, lref_2bits, nref_3bits, ref_seeds, lpos); 
        if(((lpos + 1) & 0x1f) == 0){
                base_2bits[(lpos >> 5) * 2] = lref_2bits;
                base_2bits[(lpos >> 5) * 2 + 1] = href_2bits;
                base_3bits[(lpos >> 5)] = nref_3bits;

                lref_2bits = 0;
                href_2bits = 0;
                nref_3bits = 0;
        }

        if(((spos + 1) & 0xf) == 0){
            seeds_2bits[spos >> 4] = ref_seeds;
            ref_seeds = 0;
        }

        lpos ++;
        spos ++;
    }

    index_t gpos = build_info.total_len + base_ref_off - build_info.res_len; 
    sindex_t load_size = get_size;
    seed_len = 0;
    sindex_t seeds_num = 0; 
    sindex_t seeds_off = 0;
    //if(core_id == 1){
    //    printf("num seeds %u\n", (base_ref_off + part_size) / SEED_STEP);
    //}
    int first_pass = 1;
    while(ref_off < base_ref_off + part_size){

        if(get_size + ref_off > base_ref_off + part_size){
            get_size = part_size - (ref_off - base_ref_off);
            load_size = get_size + BITS_SIZE - (get_size & 0x1f);
        }
        if(get_size > 0){
            pe_get(build_info.ref_seq + ref_off, ref_seq, load_size);
            dma_syn();
        }
        //#define FIND 2110748
        //if(ref_off + get_size > FIND && ref_off < FIND){
        //   printf("get char core id %u %c %u %u %u\n", _MYID, build_info.ref_seq[FIND], ref_off, get_size, load_size);
        //}
        for(i = 0; i < get_size; i ++){
            //if(core_id == 0)
            //    printf("%c", toupper(ref_seq[i]));
            encode_base(ref_seq[i], href_2bits, lref_2bits, nref_3bits, ref_seeds, lpos);    
            //seed_pos ++;
            if(ref_off + i + 1 >= base_ref_off + seed_length && ((ref_off + i + 1 - seed_length) % SEED_STEP == 0)){
                //if(core_id == 0 && base_ref_off/SEED_STEP +  seeds_off + seeds_num == 972487){
                //    printf("seeds salve %x %u %u\n", seed, ref_off + i, i);
                //    int kk;
                //    for(kk = 11; kk >= 0; kk --){
                //        printf("%c", ref_seq[i - kk]);
                //    }
                //    printf("\n");
                //    for(kk = 11; kk >= 0; kk --){
                //        printf("%c", build_info.ref_seq[ref_off + i - kk]);
                //    }
                //    printf("\n");
                //}
                if(seed_len >= seed_length){
                    seeds[seeds_num] = seed;
                }
                else{
                    seeds[seeds_num] = INVALID_SFLAG;    
                }
                seeds_num ++;
            }
            if(seed_len >= seed_length)
                seed_len -= 1;

            if(((lpos + 1) & 0x1f) == 0){
                //if(core_id == 0){
                //    printf("\n");
                //    printf("%x\n", lref_2bits);
                //    printf("%x\n", href_2bits);
                //}
                base_2bits[(lpos >> 5) * 2] = lref_2bits;
                base_2bits[(lpos >> 5) * 2 + 1] = href_2bits;
                base_3bits[(lpos >> 5)] = nref_3bits;

                lref_2bits = 0;
                href_2bits = 0;
                nref_3bits = 0;
            }
            if(((spos + 1) & 0xf) == 0){
                seeds_2bits[spos >> 4] = ref_seeds;
                ref_seeds = 0;
            }
            lpos ++;
            spos ++;
        }
        if(lpos >= 32){
            //if(2*(lpos /32) > 1024){
             //   printf("wrong %u %u\n", lpos, lpos/32);
            //}
            //if((gpos >> 5) * 2 <= 0x9c0)
            //cal_locked_printf("id: %d, gpos : %llu %lx\n", _MYID, gpos, (gpos >> 5) * 2) ;                                                                                                   
            pe_put(build_info.base_2bits + (gpos>>5)*2, base_2bits, 2 * (lpos >> 5) * sizeof(sindex_t));
            pe_put(build_info.base_3bits + ((gpos >> 5)), base_3bits, (lpos >> 5) * sizeof(sindex_t));
        }
        if(spos >= 16)
            pe_put(build_info.seeds_2bits + ((gpos >> 4)) , seeds_2bits, (spos >> 4) * sizeof(sindex_t));
        if(seeds_num > 0)
            pe_put(build_info.seeds_buffer[0] + base_ref_off / SEED_STEP + seeds_off, seeds, (seeds_num) * sizeof(sindex_t));
        dma_syn();

        ref_off += get_size;
        //gpos += (lpos >> 5) * BITS_SIZE;
        if(first_pass){

            gpos += build_info.res_len;
            first_pass = 0;
        }
        gpos += get_size;
        seeds_off += seeds_num;
        //if(_MYID == 0)
        //cal_locked_printf("id: %d  %llu\n", _MYID, seeds_off) ;                                                                                                   
        seeds_num = 0;
        lpos &= 0x1f;
        spos &= 0xf;
    }


    if(base_ref_off + part_size  == build_info.ref_seq_len){
        //printf("core id iiii %u %x %x\n", _MYID, lref_2bits, href_2bits);
        base_2bits[0] = lref_2bits;
        base_2bits[1] = href_2bits;

        pe_put(build_info.base_2bits + ((gpos >> 5))*2, base_2bits, 2*sizeof(sindex_t));
        pe_put(build_info.base_3bits + ((gpos >> 5)), &nref_3bits, 1*sizeof(sindex_t));

        pe_put(build_info.seeds_2bits + ((gpos >> 4)) , &ref_seeds, 1*sizeof(sindex_t));
        dma_syn();
    }
    if(base_ref_off + part_size  < build_info.ref_seq_len){
        index_t seeds_size = (base_ref_off + part_size)  / SEED_STEP - base_ref_off/SEED_STEP;
        index_t res_num_seed = seeds_size - seeds_off;

        //if(_MYID == 0)
        //    cal_locked_printf("id: %d  %llu %u\n", _MYID, seeds_off, seeds_size);  
        if(res_num_seed == 0)
            return;
        index_t res_base = SEED_STEP -  (ref_off - seed_length)  % SEED_STEP ;
        if(res_num_seed > 1)
            res_base +=  (res_num_seed - 1) * SEED_STEP;
        //res_base = res_num_seed  * SEED_STEP;

        index_t load_size = res_base + (BITS_SIZE - (res_base&0x1f));
        pe_get(build_info.ref_seq + base_ref_off + part_size, ref_seq, load_size);
        dma_syn();
        for(i = 0; i < res_base; i ++){
            encode_base(ref_seq[i], href_2bits, lref_2bits, nref_3bits, ref_seeds, lpos);    
            if(ref_off + i + 1 >= base_ref_off + seed_length && ((ref_off + i + 1 - seed_length) %SEED_STEP == 0)){
                if(seed_len >= seed_length){
                    seeds[seeds_num] = seed;
                }else{
                    seeds[seeds_num] = INVALID_SFLAG;    
                }
                seeds_num ++;
            }
            if(seed_len >= seed_length)
                seed_len -= 1;
        }
        pe_put(build_info.seeds_buffer[0] + base_ref_off / SEED_STEP + seeds_off, seeds,  seeds_num * sizeof(sindex_t));
        dma_syn();
    }
        
}

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
//#define 
void store_hash_index(index_info_t *index_info_h){
    index_info_t index_info;
    dma_declare();
    pe_get(index_info_h, &index_info, sizeof(index_info));
    dma_syn();
    index_t seeds_num = (index_info.ref_len - index_info.seed_length + 1) / SEED_STEP; 

    index_t part_seeds = (seeds_num + 2 - 1) / 2;
    part_seeds =  (part_seeds + CORE_SIZE - 1) / CORE_SIZE ;


    //index_t part_size = part_seeds * SEED_STEP * 2 + index_info.seed_length - 1;

    
    //if(base_ref_off + part_size > index_info.ref_len){
    //    part_size = index_info.ref_len - base_ref_off;
    //}
    //if(base_ref_off > index_info.ref_len)
    //    part_size = 0;
    //if(part_size == 0)
    //    return;

    part_seeds *= 2;

    index_t base_ref_off = _MYID * part_seeds; 
    index_t ref_off = _MYID * part_seeds; 

    if(base_ref_off + part_seeds  > seeds_num){
        part_seeds = seeds_num - _MYID * part_seeds;
    }
    if(base_ref_off > seeds_num)
        return;

    index_t get_size = 2048;
    sindex_t seeds[2048 + 8];
    uint16_t pos_lseeds[2048 + 8];

    sindex_t seeds_2bits[512 + 8];
    sindex_t base_2bits[512 + 8];
    sindex_t base_3bits[256 + 8];



    index_t load_16size = 0;
    index_t load_32size = 0;
    index_t seed_length = index_info.seed_length;

    sindex_t valid_seed_mask = (1UL << (seed_length * 2)) - 1; 
    sindex_t valid_flag_mask = (1UL << (seed_length)) - 1;
    
    //index_t seeds_cnt = 0;
    index_t seeds_total = 0;

    //cal_locked_printf("wrong seeds num %u %u %u\n", _MYID, seeds_total, part_seeds);
    while(seeds_total < part_seeds){
        if(get_size + seeds_total > part_seeds){
            get_size = part_seeds -  seeds_total;
        }

        sindex_t res_len = SEED_STEP - 1; 
        sindex_t mask = (1U<<res_len) - 1;

        index_t gpos = index_info.total_len + (base_ref_off + seeds_total) * SEED_STEP;
        if(gpos >= res_len){
            gpos -= res_len;
        }

        // for base_2bits, base_3bits
        index_t res_32len = gpos & 0x1f;
        load_32size = (get_size * SEED_STEP + seed_length - 1 + res_len + res_32len + 32 - 1) / 32;
        index_t pos32 = (index_info.total_len + (base_ref_off + seeds_total) * SEED_STEP)  - ((gpos >> 5) << 5); 

        // for seeds_2bits
        index_t res_16len = gpos & 0xf;
        load_16size = (get_size * SEED_STEP + res_16len + 2 * (seed_length - 1) + res_len * 2 + 16 - 1) / 16;
        index_t pos16 = index_info.total_len + (base_ref_off + seeds_total)*SEED_STEP - ((gpos >> 4) << 4); 
        //index_t pos16 = index_info.total_len + ref_off - ((gpos >> 4) << 4); 
        if(_MYID == 0){
        //    printf("pos32 %u %u\n", pos32, gpos);
        }
        if(get_size > 0){
            pe_get(index_info.seeds_2bits + (gpos >> 4), seeds_2bits, (load_16size) * sizeof(sindex_t)) ;

            pe_get(index_info.base_2bits + (gpos >> 5) * 2, base_2bits, (load_32size * 2) * sizeof(sindex_t)) ;

            pe_get(index_info.base_3bits + (gpos >> 5), base_3bits, (load_32size) * sizeof(sindex_t));
            /*
            if(_MYID ==  0)
                printf("get pos32 %u %u %u %u\n", gpos, gpos >> 5, load_32size, pos32);
            */
            dma_syn();
        }
        index_t i;
        for(i = 0; i < get_size; i ++){
            
               
                sindex_t valid_flag = ((base_3bits[pos32>>5] >> (pos32 & 0x1f)) 
                        | (base_3bits[(pos32 >> 5)+ 1] << (32 -(pos32 & 0x1f)))) & valid_flag_mask;

                if(valid_flag != 0){
                    //pos_lseeds[i] = 
                    seeds[i] = INVALID_SFLAG;
                    
                    pos16 += SEED_STEP;
                    pos32 += SEED_STEP;
                    continue;
                }

                sindex_t valid_seed = ((seeds_2bits[pos16 >>4] >> ((pos16 & 0xf) * 2)) 
                        | (seeds_2bits[(pos16 >> 4) + 1] << ((16 - (pos16 & 0xf)) * 2))) & valid_seed_mask;

                seeds[i] = valid_seed;
                assert(SEED_STEP <= 5);
                assert(SEED_STEP >= 1);

                index_t s_pos = 0;  
                if(pos32 >= res_len)
                    s_pos = pos32 - res_len;
                sindex_t shift = s_pos & 0x1f;
                s_pos >>= 5;

                sindex_t flseed = (base_2bits[s_pos * 2]  >> shift) |  (base_2bits[(s_pos + 1) * 2] << (32 - shift));
                sindex_t fhseed = (base_2bits[s_pos * 2 + 1]  >> shift) |  (base_2bits[(s_pos + 1) * 2 + 1]  << (32 - shift));
                flseed = (flseed & mask); 
                fhseed = (fhseed & mask);

                s_pos = pos32 + seed_length;
                shift = s_pos & 0x1f;
                s_pos >>= 5;

                sindex_t blseed = (base_2bits[s_pos * 2]  >> shift) |  (base_2bits[(s_pos + 1) * 2]  << (32 - shift));
                sindex_t bhseed = (base_2bits[s_pos * 2 + 1]  >> shift) |  (base_2bits[(s_pos + 1) * 2 + 1]  << (32 - shift));
                blseed = (blseed & mask); 
                bhseed = (bhseed & mask);

                flseed = flseed | (blseed << res_len);
                fhseed = fhseed | (bhseed << res_len);

                sindex_t lseed_enc = flseed | (fhseed << (res_len * 2));
                pos_lseeds[i] = lseed_enc;
/*
                if(base_ref_off + seeds_total + i == 2500){
#define FIND 10000
                    sindex_t valid_flag = ((index_info.base_3bits[FIND >>5] >> (FIND & 0x1f)) 
                        | (index_info.base_3bits[(FIND >> 5)+ 1] << (32 -(FIND & 0x1f)))) & valid_flag_mask;

                    printf("%x\n", index_info.base_3bits[FIND >> 5]);
                    printf("%x\n", (index_info.base_3bits[(FIND >> 5)+ 1]));

                    printf("%x\n", index_info.base_3bits[FIND >> 5] >> (32 - (FIND & 0x1f)));
                    printf("%x\n", (index_info.base_3bits[(FIND >> 5)+ 1] << (32 - (FIND & 0x1f))));

                    printf("pos32 %x %u %u %u %u %u\n", valid_flag, pos32, gpos, seeds_total, base_ref_off);

                    printf("%x\n", valid_flag_mask);

                    printf("%x\n", base_3bits[57]); 
                    printf("%x\n", base_3bits[58]); 

                    printf("%x\n", base_3bits[pos32>>5] >> (pos32 & 0x1f));
                    printf("%x\n", base_3bits[(pos32>>5) + 1] << (32 -(pos32 & 0x1f)));
                    valid_flag = ((base_3bits[pos32>>5] >> (pos32 & 0x1f)) 
                        | (base_3bits[(pos32 >> 5)+ 1] << (32 -(pos32 & 0x1f)))) & valid_flag_mask;
                    printf("pos32 %x %u %u %u %u %u\n", valid_flag, pos32, gpos, seeds_total, base_ref_off);
                }
*/ 
                pos16 += SEED_STEP;
                pos32 += SEED_STEP;
        }

        //ref_off += get_size;
        //if(gpos > res_len)
        //    gpos += res_len;
        //gpos += get_size;
        //seeds_total += seeds_cnt;

        //if(index_info.store_flag == 1){
            pe_put(index_info.seeds_buffer + base_ref_off + seeds_total, seeds, get_size * sizeof(sindex_t)) ;
        //}
        //if(index_info.store_flag == 2){
            if(get_size > 0)
            pe_put(index_info.pos_lseeds + base_ref_off + seeds_total, pos_lseeds, (get_size + (get_size & 1))  * sizeof(uint16_t)) ;
        //}
        dma_syn();

        seeds_total += get_size;
    }

    //cal_locked_printf("wrong seeds num %u %u %u\n", _MYID, seeds_total, part_seeds);
    //if(seeds_total != part_seeds){
    //    //cal_locked_printf("id: %d, gpos : %llu %lx\n", _MYID, gpos, (gpos >> 5) * 2) ;                                                                                                   
    //    cal_locked_printf("wrong seeds num %u %u %u\n", _MYID, seeds_total, part_seeds);
    //}

}
