# SWMapper

SWMapper is a scalable and efficient read mapper for the Sunway TaihuLight supercomputer. 

## Build

```bash
git clone git@github.com:Xu-Kai/SWMapper.git 

cd SWMapper

make
```

## Usage
### Build hash index
We provide two methods to buid hash index. One is seiral version. The other is distributed version. 


#### Serial version

Parameters:

```
-i input reference

-p prefix index name

-s  [opt, default:12] seed length. We recommend the default value. 

-h print help information
```

Example:

```bash
bsub -debug -I -b -m 1 -q q_sw_expr -share_size 7200 -host_stack 256 -n 1 -cgsp 64 -o print.out ./serial_index.out -ihg38.fa  -phg38
```


#### Distributed version

Parameters:

```
-i input reference

-p prefix name of index file

-s  [opt, default:12] seed length. We recommend the default value. 

-h print help information

```


Example:

```bash
bsub -debug -I -b -m 1 -q q_sw_expr -share_size 7200 -host_stack 256 -n 4 -cgsp 64 -o print.out ./mpi_index.out -ihg38.fa  -phg38
```

### Mapping

Parameters:

```
-l max read length (less than 250)

-o prefix name of output file 

-e error size

-p prefix name of index file

-s [opt, default:12] seed length

-i single data input

-f file which includes input data file names
```

Example:

```bash
bsub -debug -I -b -m 1 -q q_sw_share -share_size 7200 -host_stack 256 -n 8 -cgsp 64 -o print.out ./mpi_index.out -ihg38.fa  -phg38
```

## Bug Report

All bug reports, comments and suggestions are welcome. 

Feel free to open a new issue.
