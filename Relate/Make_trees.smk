import os
##input
##input info extraction
workdir: os.getcwd()
include: os.getcwd()+'/config.sk'

chrall = [str(i) for i in range(1, 23)]

rule all:
    input:
        expand("pop_size/{pop_to_extract}/{pop_to_extract}_popsize.coal", pop_to_extract=POP_TO_EXTRACT)


#prepare input files

rule ConvertFromVcf:
    input:
        vcf+".vcf.gz"
    output:
        multiext("input_relate/haps_sample/allpop_chr{chr}",".haps.gz",".sample.gz")
    params:
        prefix= vcf,
        haps= "input_relate/haps_sample/allpop_chr{chr}.haps",
        sample= "input_relate/haps_sample/allpop_chr{chr}.sample",
        dir = Relate_Dir
    resources:
        mem_mb = 40000
    shell:
        """
		{params.dir}/bin/RelateFileFormats \
                	--mode ConvertFromVcf \
                	--haps {params.haps} \
                	--sample {params.sample} \
                	-i {params.prefix}
		gzip {params.haps}
		gzip {params.sample}
		"""
rule CreatePopLabel:
    input:
        sample="input_relate/haps_sample/allpop_chr1.sample.gz",
        pop_to_sort= popfile
    output:
        order="input_relate/order_relate.txt",
        poplabel="input_relate/poplabel.txt"
    resources:
        mem_mb = 4000
    shell:
        """
		zcat {input.sample} |cut -f1|grep -v -x "0\|ID_1"> {output.order}
		echo "sample population group sex" > {output.poplabel}
		while read LINE; do
    			awk -v var=$LINE '$1==var' {input.pop_to_sort} >> {output.poplabel}
		done < {output.order}
		"""

rule CreateRelateInput:
    input:
        poplabel="input_relate/poplabel.txt",
        ancestor=ancestor_file,
        mask= mask_file,
        haps="input_relate/haps_sample/allpop_chr{chr}.haps.gz",
        sample="input_relate/haps_sample/allpop_chr{chr}.sample.gz"
    params:
        prefix ="input_relate/RELATE_input_chr{chr}",
        dist="input_relate/RELATE_input_chr{chr}.dist",
        dir = Relate_Dir
    output:
        haps="input_relate/RELATE_input_chr{chr}.haps.gz",
        sample="input_relate/RELATE_input_chr{chr}.sample.gz",
        annot="input_relate/RELATE_input_chr{chr}.annot",
        dist="input_relate/RELATE_input_chr{chr}.dist.gz"
    resources:
        mem_mb = 100000
    shell:
        """
		{params.dir}/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps {input.haps} --sample {input.sample} --ancestor {input.ancestor} --mask {input.mask} --poplabels {input.poplabel} -o {params.prefix}
		"""
#make tree

rule make_tree:
    input:
        haps="input_relate/RELATE_input_chr{chr}.haps.gz",
        sample="input_relate/RELATE_input_chr{chr}.sample.gz",
        annot="input_relate/RELATE_input_chr{chr}.annot",
        dist="input_relate/RELATE_input_chr{chr}.dist.gz",
        map=map_file
    output:
        multiext("tree_chr{chr}",".mut.gz",".anc.gz")
    params:
        prefix="tree_chr{chr}",
        mut="tree_chr{chr}.mut",
        anc="tree_chr{chr}.anc",
        dir = Relate_Dir
    resources:
        threads = 10,
        mem_mb= 100000 #This is the total memory. You want that mem_mb/threads > 5000 because it is the minimum mem per cpu needed by relate
    shell:
        """
		{params.dir}/scripts/RelateParallel/RelateParallel.sh\
		            -m 1.25e-8 \
		            -N 30000 \
		            --haps {input.haps} \
		            --sample {input.sample} \ 
		            --annot {input.annot} \
		            --map {input.map} \
		            --dist {input.dist} \
		            --seed 1 \
		            -o {params.prefix} \ 
		            --threads {resources.threads} 
		
		gzip {params.mut}
		gzip {params.anc}
		gzip {input.annot}
		"""

rule move_tree:
    input:
        mut="tree_chr{chr}.mut.gz",
        anc="tree_chr{chr}.anc.gz"
    output:
        mut="trees_all/tree_chr{chr}.mut.gz",
        anc="trees_all/tree_chr{chr}.anc.gz"
    resources:
        mem_mb = 2000
    shell:
        """
		mv {input.mut} {output.mut}
		mv {input.anc} {output.anc}
		"""

#extract tree

rule extract_pop:
    input:
        mut="trees_all/tree_chr{chr}.mut.gz",
        anc="trees_all/tree_chr{chr}.anc.gz",
        poplabel="input_relate/poplabel.txt"
    output:
        poplabel="extracted_trees/{pop_to_extract}/{pop_to_extract}_chr{chr}.poplabels",
        zip=multiext("extracted_trees/{pop_to_extract}/{pop_to_extract}_chr{chr}",".mut.gz", ".anc.gz")
    params:
        prefix= "extracted_trees/{pop_to_extract}/{pop_to_extract}_chr{chr}",
        mut="extracted_trees/{pop_to_extract}/{pop_to_extract}_chr{chr}.mut",
        anc="extracted_trees/{pop_to_extract}/{pop_to_extract}_chr{chr}.anc",
        dir = Relate_Dir
    resources:
        mem_mb= 50000
    shell:
        """
		{params.dir}/bin/RelateExtract \
			--mode SubTreesForSubpopulation \
			--anc {input.anc} \
			--mut {input.mut} \
			--poplabels {input.poplabel} \
			--pop_of_interest {wildcards.pop_to_extract} \
			-o {params.prefix}

		gzip {params.mut}
		gzip {params.anc}
		"""

rule pop_size:
    input:
        expand("extracted_trees/{pop_to_extract}/{pop_to_extract}_chr{chr}.mut.gz",chr=chrall, pop_to_extract=POP_TO_EXTRACT),
        expand("extracted_trees/{pop_to_extract}/{pop_to_extract}_chr{chr}.anc.gz",chr=chrall, pop_to_extract=POP_TO_EXTRACT),
        poplabel="extracted_trees/{pop_to_extract}/{pop_to_extract}_chr22.poplabels"
    output:
        "pop_size/{pop_to_extract}/{pop_to_extract}_popsize.coal"
    params:
        input= "extracted_trees/{pop_to_extract}/{pop_to_extract}",
        output="pop_size/{pop_to_extract}/{pop_to_extract}_popsize",
        dir = Relate_Dir
    resources:
        threads = 10,
        mem_mb = 60000 #This is the total memory. You want that mem_mb/threads > 5000 because it is the minimum mem per cpu needed by relate
    shell:
        """
		{params.dir}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
			-i {params.input}\
			-m 1.25e-8 \
			--bins 2,7,0.1 \
			--poplabels {input.poplabel} \
			-o {params.output} \
			--threads {resources.threads} \
			--first_chr 1 \
			--last_chr 22
		"""
