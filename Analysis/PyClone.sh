set -eo pipefail

export LD_LIBRARY_PATH=/share/public/software/llvm/6.0.1/lib/:$LD_LIBRARY_PATH

sample=$1
mut=$2
cnv=$3
purity=$4

dir=$PWD
mkdir -p $dir
cd $dir
perl /share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/result/CCF/PyClone/Purple/merge.cnv_AF.pl $mut $cnv $sample.in
/share/work1/wangrr/local/simple/bin/PyClone setup_analysis --in_files $sample.in --working_dir $dir/ --tumour_contents $(tail -1 $purity | cut -f 1) --samples $sample --num_iters 1000 --prior major_copy_number --init_method connected
/share/work1/wangrr/local/simple/bin/PyClone run_analysis --seed 5 --config_file config.yaml
/share/work1/wangrr/local/simple/bin/PyClone build_table --config_file config.yaml --out_file $sample.pyclone_cluster.txt --table_type old_style --max_clusters 50
