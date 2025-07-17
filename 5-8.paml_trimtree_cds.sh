#!/bin/bash
#SBATCH -c 1
#SBATCH --job-name=PAML
#SBATCH --time=14-00:00:00
#SBATCH --partition=fos
#SBATCH --qos=fos

# ==== 🔧 配置 ====
PPN=1
fasta_dir="trimal_cds_rerun2_edit/"
output_dir="trimmed_trees"
paml_script=~/shell/prune_species_tree.py
whole_species_tree=/data/scratch/projects/punim2285/ziwei_data/TEST-8genes/Craig_Kumar_Hedges_final_timetree.nwk
paml_inputs=paml_inputs_branch-site_ma1_vs_ma2_batch #paml_inputs_m012 #paml_inputs #paml_inputs_real

# =================

# 加载模块
module load Biopython
module load PAML

# 确保所需的 Python 包已安装
#pip install --user ete3 biopython
#conda activate /data/scratch/projects/punim2285/ziwei_data/software/anaconda3/envs/rna/
#rename -v 's/fasta/fas/g' ${fasta_dir}/*

# 执行物种树修剪脚本
#python ${paml_script} \
#  --fasta_dir ${fasta_dir} \
#  --tree_file ${whole_species_tree} \
#  --output_dir ${output_dir}

#python ~/shell/replace_fasta_species_names.py \
#  --simplified_log ${output_dir}/simplified_matches.log \
#  --fuzzy_log ${output_dir}/fuzzy_matches.log \
#  --fasta_dir ${fasta_dir} \
#  --output_dir updated_fasta_for_tree

#设置前景支
#python ~/shell/rename_nwk.py

#准备paml
paml_inputs=paml_inputs_m012_Lemuroidea #omnivore_deepest #Colobinae#1_Lemuroidea#1 #_Colobinae #Cercopithecinae/ #carnivore #herbivore #_omnivore #_Colobinae #_Cercopithecinae #paml_inputs_m012_Lemuroidea #paml_inputs_m012
labelled_trees=labelled_trees_Lemuroidea #omnivore #_Colobinae_Lemuroidea #_Colobinae #Cercopithecinae/ #labeled_diet_trees/carnivore_labeled_trees/ #labelled_trees_Colobinae #_Cercopithecinae #labelled_trees #labelled_trees_Lemuroidea
#python ~/shell/prepare_paml_inputs_cds.py \
#  --fasta_dir updated_fasta_for_tree \
#  --tree_dir trimmed_trees \
#  --labelled_tree_dir ${labelled_trees} \
#  --output_dir ${paml_inputs}
#cp /data/scratch/projects/punim2285/ziwei_data/software/paml/dat/wag.dat ${paml_inputs}

# 运行 PAML
#cd ${paml_inputs}
#for ctl in *.ctl; do
#  codeml "$ctl"
#done
# 多批次运行PAML,ppn设置越小越快
#python ~/shell/submit_codeml_in_batches.py --ctl_dir ${paml_inputs} --cpu ${PPN}

#计算LTR并用卡方检验显著性
#python ~/shell/run_paml_lrt_m012.py --input_dir ${paml_inputs} --output_tsv lrt_summary_branch.tsv

#运行 branch site modela1.vs.a2
#paml_inputs=paml_inputs_branch-site.a1a2_herbivore_shallowest #_Colobinae_Lemuroidea #_Colobinae #_Cercopithecinae #Lemuroidea #Colobinae #paml_inputs_branch-site.a1a2_Cercopithecinae
#paml_inputs=paml_inputs_branch-site.a1a2_omnivore #_herbivore #_carnivore #_Lemuroidea #_Cercopithecinae #_Colobinae
#labelled_trees=labelled_trees_herbivore #omnivore_shallowest #_Colobinae_Lemuroidea #_Colobinae #_Cercopithecinae #_Colobinae #_Lemuroidea
#准备 paml 输入文件
#python ~/shell/prepare_paml_inputs_cds_branch-site.a1_vs.a2.py --fasta_dir updated_fasta_for_tree --labelled_tree_dir ${labelled_trees} --output_dir ${paml_inputs}
# 多批次运行PAML,ppn设置越小越快
#fos节点
#python ~/shell/submit_codeml_in_batches.py --ctl_dir ${paml_inputs} --cpu ${PPN}
#非fos节点
#python ~/shell/submit_codeml_in_batches_nofos.py --ctl_dir ${paml_inputs} --cpu ${PPN}
#计算LTR并用卡方检验显著性
#python ~/shell/run_paml_lrt_branch-site.py --input_dir ${paml_inputs} --output_tsv branch_site_LRT_results.tsv

#运行 site modela1.vs.a2
#paml_inputs=paml_inputs_site_m1a_m2a
#M1a vs M2a
#python ~/shell/prepare_paml_inputs_cds_site-model.ma1_vs.ma2.py \
#  --fasta_dir updated_fasta_for_tree \
#  --tree_dir trimmed_trees \
#  --output_dir ${paml_inputs}
#运行PAML,ppn设置越小越快
python ~/shell/submit_codeml_in_batches_nofos.py --ctl_dir ${paml_inputs} --cpu ${PPN}
#计算LTR并用卡方检验显著性
#python ~/shell/run_paml_lrt-site.py --input_dir ${paml_inputs}   --prefix1 modelM1a --prefix2 modelM2a --df 2 --output_tsv lrt_m1a_vs_m2a.tsv

#运行 site modelm7.vs.m8
#paml_inputs=paml_inputs_site_m7_m8
#M7 vs M8
#python ~/shell/prepare_paml_inputs_cds_site-model.m7_vs.m8.py \
#  --fasta_dir updated_fasta_for_tree \
#  --tree_dir trimmed_trees \
#  --output_dir ${paml_inputs}
#运行PAML,ppn设置越小越快
#python ~/shell/submit_codeml_in_batches_nofos.py --ctl_dir ${paml_inputs} --cpu ${PPN}
#计算LTR并用卡方检验显著性
#python ~/shell/run_paml_lrt_site.py --input_dir ${paml_inputs} --prefix1 modelM7 --prefix2 modelM8 --df 2 --output_tsv lrt_m7_vs_m8.tsv
