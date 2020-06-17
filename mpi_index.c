/*************************************************************************
  > File Name: mpi_index.c
  > Author: Xu Kai
  > Created Time: Wed 26 Feb 2020 05:20:46 PM CST
 ************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<stdio.h>
#include<stdint.h>
#include<fcntl.h>
//#include "htslib/kseq.h"
#include "kseq.h"
#include<assert.h>
#include <athread.h>
#include<getopt.h>
#include<mpi.h>

//#define MPE_PROFILE
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

index_t ref_total_length;

index_t ref_offset;
index_t ref_end;
index_t g_ind;
sindex_t base_2bits_s0, base_2bits_s1;

sindex_t *seeds_check;

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
KSEQ_INIT(int,read)
  kseq_t *kseq_open(const char *path){
    int fd = open(path, O_RDONLY);
    return kseq_init(fd);
  }



void get_base_2bits_s(index_t off, index_t len, char *seq){
  index_t t = 0; 
  if(off >= SEED_STEP - 1){
    t = SEED_STEP - 1; 
  }else{
    t = off;
  }
  index_t i;
  for(i = 0; i < t; i ++){
    char base = seq[off - i]; 
    switch(base){
      case 'A':
        base_2bits_s0 = 0UL << (31 - i); 
        base_2bits_s1 = 0UL << (31 - i); 
        break;
      case 'C':
        base_2bits_s0 = 0UL << (31 - i); 
        base_2bits_s1 = 1UL << (31 - i); 
        break;
      case 'G':

        base_2bits_s0 = 1UL << (31 - i); 
        base_2bits_s1 = 0UL << (31 - i); 
        break;
      case 'T':
        base_2bits_s0 = 1UL << (31 - i); 
        base_2bits_s1 = 1UL << (31 - i); 
        break;
      default:
        base_2bits_s0 = 0UL << (31 - i); 
        base_2bits_s1 = 0UL << (31 - i); 
        break;
    }

  }



}



void printf_help(){
  printf("-i\tinput reference\n");
  printf("-p\tindex prefix\n");
  printf("-s\t[opt, default:12] seed length\n");
  exit(0);
}
char get_base(char *kseq){

  if(chromo_position <  chromo_length){
    g_ind ++;
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
        //++++++++++++++++++++++++++++++
        //seed = (seed >> 2) | (3UL << ((seed_length - 1) * 2));
        //cur_seed_length ++;

        cur_seed_length = 0;
        seed = 0;
        //-----------------------------
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
    if((g_ind - seed_length)%SEED_STEP == 0){
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

  if(g_ind < ref_offset + res_len){
    shift = 29;
    flseed = base_2bits_s0 >> shift;
    fhseed = base_2bits_s1 >> shift;
  }

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
  index_t ref_seed_end;
  if(ref_end + seed_length - 1 < ref_total_length)
    ref_seed_end = ref_end + seed_length - 1; 
  else
    ref_seed_end = ref_total_length;

  
  index_t j = 0;
  for(g_ind = ref_offset; g_ind < ref_seed_end - seed_length + 1; g_ind += SEED_STEP){
    valid_seed = get_seed_from_enc_bits(j);
    j += SEED_STEP;
    if(valid_seed != SINVALID_SEED){
      seeds_check[valid_seed] ++;
      index_t ind = hash_index[valid_seed];
      //if(ind > valid_seed_size){
      //    printf("valid seed  : %llu %llu\n", ind, valid_seed_size);
      //}
      positions[ind] = g_ind;
      pos_lseeds[ind] = lseed_enc;
      hash_index[valid_seed] ++;
      seed_cnt ++;
    }
  }

 // printf("valid seed size : %llu\n", seed_cnt);
}


void get_ref_length(){

  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  if(my_id == 0){
    kseq_t *kref = kseq_open(ref_path);
    sprintf(out_prefix, "%s%s", prefix, ".ref_info");
    FILE* ri_ptr = fopen(out_prefix, "w");
    int l = 0;
    while(l >= 0){
#ifdef MPE_PROFILE
      GPTLstart("diskio");
#endif
      l = kseq_read(kref);

#ifdef MPE_PROFILE
      GPTLstop("diskio");
#endif
      if(l < 0)
        break;
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

      chromo_cnt += 1;
    }
    printf("chrom num %u\n", chromo_cnt);

    fclose(ri_ptr);

    kseq_destroy(kref);
  }
#ifdef MPE_PROFILE 
  GPTLstart("MPIcomm");
#endif
  MPI_Bcast(&chromo_cnt, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(ref_length,  chromo_cnt, MPI_UINT32_T, 0, MPI_COMM_WORLD);
#ifdef MPE_PROFILE
  GPTLstop("MPIcomm");
#endif
  int i = 0;
  ref_total_length = 0;
  index_t total_seeds_n = 0;
  for(;i < chromo_cnt; i ++){
    ref_total_length += ref_length[i];
  }
  total_seeds_n += (ref_total_length - seed_length + 1)/4;
  //if(my_id == 0)
  //  printf("total seeds n %llu\n", total_seeds_n);
}
void encode_ref_to_bits(){

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

  index_t ref_per_cg  = ref_total_length / mpi_size;
  if(ref_per_cg % 32 != 0){
    ref_per_cg = ref_per_cg  + (32 - ref_per_cg % 32);
  }

  ref_offset = my_id * ref_per_cg;
  ref_end = (my_id + 1) * ref_per_cg;

  if(my_id == mpi_size - 1)
    ref_end = ref_total_length;

  index_t chromo_start = 0;
  index_t chromo_end = 0;


  int i;
  index_t tmp = 0;
  for( i = 0; i < chromo_cnt; i ++){
    tmp += ref_length[i];
    if(tmp >= ref_offset){
      chromo_start = i;
      break;
    }
  }
  tmp = 0;
  for(i = 0; i < chromo_cnt; i ++){
    tmp += ref_length[i];
    if(tmp >= ref_end){
      chromo_end = i;
      break;
    }
  }


  index_t length_per_cg = ref_end  - ref_offset;
  index_t seed_kinds = 1 << (seed_length * 2);

  base_2bits = (sindex_t*)malloc(256 + length_per_cg  / 16 * sizeof(sindex_t));
  seeds_2bits = (sindex_t*)malloc(256 + length_per_cg  / 16 *  sizeof(sindex_t));
  base_3bits = (sindex_t*)malloc(256 + length_per_cg  / 32 *  sizeof(sindex_t));

  seed_count = (sindex_t*)malloc( (1 << (seed_length * 2) ) * sizeof(sindex_t));
  sindex_t *seed_offset = (sindex_t*)malloc((1 << (seed_length * 2)) * sizeof(sindex_t));
  seeds_check = (sindex_t*)malloc((1 << (seed_length * 2)) * sizeof(sindex_t));

  hash_index = (sindex_t*)malloc(((1 << (seed_length * 2)) +  1)  * sizeof(sindex_t)); 

  memset(seed_count, 0, seed_kinds * sizeof(sindex_t));
  memset(seed_offset, 0, seed_kinds * sizeof(sindex_t));
  memset(seeds_check, 0,  seed_kinds * sizeof(sindex_t));

  kseq_t *kref = kseq_open(ref_path);
  index_t l;

  enc_lbits = 0;
  enc_hbits = 0;
  enc_seeds_2bits = 0;
  enc_3bits = 0;

  valid_seed_size = 0; 
  position = 0;

  index_t ref_read_len = 0;
  chromo_cnt = 0;

  //printf("my id se %d %llu %llu\n", my_id, chromo_start, chromo_end);
  index_t com_length = 0;
  g_ind = ref_offset;
  while(l = kseq_read(kref) >= 0){

    chromo_length = kref->seq.l;

    //cur_seed_length = 0;
    //seed = 0;
    if(chromo_start <= chromo_cnt){

      chromo_position = 0; 
      if(chromo_start == chromo_cnt){
        chromo_position = ref_offset - ref_read_len; 
        get_base_2bits_s(chromo_position, kref->seq.l, kref->seq.s);
      }

      chromo_length = kref->seq.l; 


      if(chromo_length + ref_read_len > ref_end){
        chromo_length = ref_end - ref_read_len + chromo_position;
        if(chromo_start == chromo_cnt)
          chromo_length =  ref_end - ref_read_len; 
      }

      //if(chromo_end == chromo_cnt)
      //    chromo_length = ref_end - ref_read_len;

      //if(chromo_start == chromo_cnt)
      //    chromo_length =  ref_end - ref_read_len; 

      com_length += chromo_length - chromo_position;
      if(chromo_length + seed_length - 1 <= kref->seq.l){
        chromo_length = chromo_length + seed_length - 1;
      }else{
        chromo_length = kref->seq.l;
      }

      chromo_hash_index(kref->seq.s);
      //printf("mpi id: %d name: %s seq len: %llu\n", my_id, kref->name.s, kref->seq.l); 

    }

    if(chromo_cnt == chromo_end)
      break;

    ref_read_len += ref_length[chromo_cnt];
    chromo_cnt += 1;

  }

  base_2bits[2*((position + 31)/32 - 1)] = enc_lbits;
  base_2bits[2*((position + 31)/32 - 1) + 1] = enc_hbits;
  base_3bits[(position + 31)/32 - 1] = enc_3bits;
  seeds_2bits[(position + 15)/16 - 1] = enc_seeds_2bits;

  //printf("com length %d %llu %llu\n", my_id, com_length, length_per_cg);
  //printf("chrom num %u\n", chromo_cnt);
  index_t total_seed_count = 0;
  index_t ii;



  for(ii= 0; ii < (1<< (seed_length * 2)); ii ++ )
    total_seed_count += seed_count[ii];

  index_t sum_seeds = total_seed_count;
  MPI_Reduce(&total_seed_count, &sum_seeds, 1, MPI_UINT64_T,MPI_SUM,0,MPI_COMM_WORLD);


  //if(my_id == 0){
    printf("valid seed size : %llu %llu %llu\n", sum_seeds, valid_seed_size, total_seed_count);
  //}

  //base_2bits_s0 = 0;
  //base_2bits_s1 = 0;

  positions = (sindex_t*)malloc((valid_seed_size + 8) * sizeof(sindex_t)); 
  pos_lseeds = (half_word_t*)malloc((valid_seed_size + 8) * sizeof(half_word_t)); 

  store_positions();

  MPI_File file_base, file_nbase;


  char out_prefix[200];
  int rc;
  MPI_Status m_status;
  index_t w_count;
#ifdef MPE_PROFILE
  GPTLstart("diskio");
#endif
  sprintf(out_prefix, "%s%s", prefix, ".base_2bits");
  w_count = ((ref_end  - ref_offset + 31) / 32) * 2;

  rc = MPI_File_open(MPI_COMM_WORLD,  out_prefix, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &file_base); 
  MPI_File_write_at(file_base, ref_offset / 16 * sizeof(sindex_t), base_2bits,  w_count, MPI_UINT32_T,  &m_status); 
  MPI_File_close(&file_base);

  sprintf(out_prefix, "%s%s", prefix, ".base_3bits");
  w_count = ((ref_end  - ref_offset + 31) / 32);

  rc = MPI_File_open(MPI_COMM_WORLD,  out_prefix, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &file_nbase); 
  MPI_File_write_at(file_nbase, ref_offset / 32 * sizeof(sindex_t), base_3bits,  w_count, MPI_UINT32_T,  &m_status); 
  MPI_File_close(&file_nbase);

#ifdef MPE_PROFILE
  GPTLstop("diskio");
#endif

  for(ii = 0; ii < seed_kinds; ii ++){
    if(seed_count[ii] != seeds_check[ii])
      printf("my id %d wrong! %u %u\n", my_id, seed_count[ii], seeds_check[ii]);
  }

#ifdef MPE_PROFILE
  GPTLstart("MPIcomm");
#endif
  MPI_Allreduce(seed_count, seed_offset, seed_kinds, MPI_UINT32_T, MPI_SUM,
      MPI_COMM_WORLD);

#ifdef MPE_PROFILE
  GPTLstop("MPIcomm");
#endif


  sindex_t *seed_all_count = (sindex_t*)malloc( (seed_kinds + seed_kinds) * sizeof(sindex_t));

  sindex_t recive_size[512];
  sindex_t disp[512];

  sindex_t seeds_mpi[512];
  sindex_t seeds_mpi_off[512];


  index_t jj = 0;
  index_t kk = 0;

  index_t seeds_start = my_id * (seed_kinds / mpi_size);
  index_t seeds_end = (my_id + 1)* (seed_kinds / mpi_size);
  index_t seed_n = seed_kinds / mpi_size;
  if(my_id == mpi_size - 1){ 
    seeds_end = seed_kinds;
    seed_n = seeds_end - seeds_start;
  }


  sindex_t seeds_ns[512];
  sindex_t seeds_ns_off[512];

  seeds_mpi_off[0] = 0;
  seeds_ns_off[0] = 0;
  for(jj = 0; jj < mpi_size; jj ++){
    index_t off = jj * (seed_kinds/mpi_size);
    index_t l_end = seed_kinds/mpi_size;
    if(jj == mpi_size - 1)
      l_end = seed_kinds - off;
    seeds_ns[jj] = l_end; 
    seeds_mpi[jj] = 0;
    for(kk = 0; kk  < l_end; kk++){
      seeds_mpi[jj] += seed_count[off + kk];
    }
    seeds_mpi_off[jj + 1] = seeds_mpi_off[jj] + seeds_mpi[jj];
    seeds_ns_off[jj + 1] = seeds_ns_off[jj] + seeds_ns[jj];
  }
  sindex_t re_ns[512],  re_ns_off[512];
  re_ns_off[0] = 0;
  for(i = 0; i < mpi_size; i ++){
    re_ns[i] = seed_n; 
    re_ns_off[i + 1] = re_ns_off[i] + seed_n;
  }

#ifdef MPE_PROFILE
    GPTLstart("MPIcomm");
#endif
  MPI_Alltoallv(seed_count, seeds_ns, seeds_ns_off, MPI_UINT32_T, seed_all_count, re_ns, re_ns_off, MPI_UINT32_T, MPI_COMM_WORLD);

#ifdef MPE_PROFILE
  GPTLstop("MPIcomm");
#endif

  sindex_t *recive_pos;
  half_word_t *recive_lseeds;
  sindex_t *recive_adj_pos;
  half_word_t *recive_adj_lseeds;

      for(jj = 0; jj < mpi_size; jj ++){
        index_t off = jj * seed_n;
        for(kk = 0; kk  < seed_n; kk++){
          recive_size[jj] += seed_all_count[off + kk];
        }
        disp[jj + 1] = disp[jj] + recive_size[jj];
      }
      recive_pos = (sindex_t*)malloc(disp[mpi_size] * sizeof(sindex_t));
      recive_lseeds = (half_word_t*)malloc(disp[mpi_size] * sizeof(half_word_t));

      recive_adj_pos = (sindex_t*)malloc(disp[mpi_size] * sizeof(sindex_t));
      recive_adj_lseeds = (half_word_t*)malloc(disp[mpi_size] * sizeof(half_word_t));


#ifdef MPE_PROFILE
    GPTLstart("MPIcomm");
#endif
  MPI_Alltoallv(positions, seeds_mpi, seeds_mpi_off, MPI_UINT32_T, recive_pos, recive_size, disp, MPI_UINT32_T, MPI_COMM_WORLD);
  MPI_Alltoallv(pos_lseeds, seeds_mpi, seeds_mpi_off, MPI_UINT16_T, recive_lseeds, recive_size, disp, MPI_UINT16_T, MPI_COMM_WORLD);

#ifdef MPE_PROFILE
  GPTLstop("MPIcomm");
#endif

  index_t adj_ind = 0;
  for(ii = 0; ii < seed_n; ii ++){
    for(jj = 0; jj < mpi_size; jj ++){
      for(kk = 0; kk < seed_all_count[jj * seed_n + ii]; kk ++){
        recive_adj_pos[adj_ind] = recive_pos[disp[jj]];  
        recive_adj_lseeds[adj_ind] = recive_lseeds[disp[jj]];  
        adj_ind ++;
        disp[jj] ++;
      } 
    }
  }

  hash_index[0] = 0; 
  for(i = 0; i < 1 << (seed_length * 2); i ++){
    hash_index[i+1] = hash_index[i] + seed_offset[i];
  }

  //if(my_id == 0)
  //  printf("hash index total %llu\n", hash_index[seed_kinds]);




#ifdef MPE_PROFILE
  GPTLstart("diskio");
#endif
  MPI_File m_file;
  w_count = hash_index[seeds_end] - hash_index[seeds_start];

  sprintf(out_prefix, "%s%s", prefix, ".positions");
  rc = MPI_File_open(MPI_COMM_WORLD,  out_prefix, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &m_file); 
  MPI_File_write_at(m_file, hash_index[seeds_start] * sizeof(sindex_t), recive_adj_pos,  w_count, MPI_UINT32_T,  &m_status); 
  MPI_File_close(&m_file);

  sprintf(out_prefix, "%s%s", prefix, ".pos_lseeds");
  rc = MPI_File_open(MPI_COMM_WORLD,  out_prefix, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &m_file); 
  MPI_File_write_at(m_file, hash_index[seeds_start] * sizeof(half_word_t), recive_adj_lseeds, w_count, MPI_UINT16_T,  &m_status); 
  MPI_File_close(&m_file);


  if(my_id == 0){

    FILE *fptr; 
    sprintf(out_prefix, "%s%s", prefix, ".hash_count");
    fptr = fopen(out_prefix, "wb");
    fwrite(seed_offset, sizeof(sindex_t), (1 << (seed_length * 2)), fptr);
    fclose(fptr);

    sprintf(out_prefix, "%s%s", prefix, ".hash_offset");
    fptr = fopen(out_prefix, "wb");
    fwrite(hash_index, sizeof(sindex_t), (1 << (seed_length * 2)) + 1, fptr);
    fclose(fptr);

  }


#ifdef MPE_PROFILE
  GPTLstop("diskio");
#endif
}


int main(int argc, char *argv[]){
  char *optstring  = "h::i:s::p:";
  int opt;
  MPI_Init(NULL, NULL);
#ifdef MPE_PROFILE
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
  get_ref_length();
  encode_ref_to_bits();

#ifdef MPE_PROFILE
  fem_swlu_prof_stop();
  fem_swlu_prof_print();
#endif

#ifdef MPE_PROFILE
  GPTLstop("total time");
  GPTLpr(0);
#endif

  MPI_Finalize();
}
