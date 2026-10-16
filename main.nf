process FASTQC1 {
        tag "FASTQC1 on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/qc/", mode:'copy'
        container "staphb/fastqc:0.12.1"

        label "l_cpu"
        label "l_mem"
 
    input:
    tuple val(name), val(sample), path(read1), path(read2)

    output:
    tuple val(name), val(sample),
          path("${name}_R1_fastqc.html"),
          path("${name}_R2_fastqc.html"),
          path("${name}_R1_fastqc.zip"),
          path("${name}_R2_fastqc.zip")

    script:
    """
    echo FASTQC1 $name
    fastqc $read1 $read2
    """        
}

process MULTIQC1 {
        tag "MULTIQC1 on ${run_name} using $task.cpus CPUs and $task.memory memory"
        container "staphb/multiqc:1.30"

        input:
        tuple val(run_name), path(fastqc_files)

        publishDir { "${params.outDirectory}/${run_name}/qc/" }, mode: 'copy'

        output:
        path("multiqc_report_rawdata.html")

        script:
        """
        echo Running MultiQC for $run_name
        ls -lh
        multiqc . -n multiqc_report_rawdata.html
        """

}

process ALIGN {
	tag "ALIGN on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
      
        label "l_cpu"
        label "l_mem"
 
	
	input:
	tuple val(name), val(sample), path(fwd), path(rev)

	output:
    tuple val(name), val(sample), path("${name}.mdup.sorted.bam"), path("${name}.mdup.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo ALIGN $name
	source activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex}.fa $fwd $rev | samblaster |samtools view -Sb - | sambamba sort /dev/stdin -o ${name}.mdup.sorted.bam
        samtools index ${name}.mdup.sorted.bam ${name}.mdup.sorted.bai
	"""
}

process FLAGSTAT {

        tag "flagstat on $name using $task.cpus CPUs and $task.memory memory"
        label "m_mem"
        label "s_cpu"
        publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.flagstat")


        script:
        """
        source activate bwa
        samtools flagstat $bam > ${name}.flagstat
        """
}


process QUALIMAP {

        tag "QUALIMAP on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
        container "staphb/qualimap:2.3"
        label "l_cpu"
        label "xl_mem"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.qualimap")

        script:
        """
        echo QUALIMAP $name
        qualimap --java-mem-size=16G bamqc -bam $bam --feature-file ${params.varbed} -outformat PDF -outfile ${name}.qualimap.pdf -outdir ${name}.qualimap 
        """
}

process GATK {
       tag "GATK on $name"
       publishDir "${params.outDirectory}/${sample.run}/mutect/", mode:'copy'
        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.mutect.vcf"), path("${name}.mutect-filtered.vcf")

        script:
        """
        echo GATK $name
        source activate gatk4610
        gatk --java-options "-Xmx4g" Mutect2 -R ${params.ref}.fa -I $bam -L ${params.varbed} --max-mnp-distance 0  --dont-use-soft-clipped-bases true -A StrandBiasBySample -minimum-mapping-quality 0 -O ${name}.mutect.vcf
 
        gatk --java-options "-Xmx4g" FilterMutectCalls -R ${params.ref}.fa -V ${name}.mutect.vcf -O ${name}.mutect-filtered.vcf
        """
}


process VAFaNORMALIZACE1 {
        tag "VAFaNORMALIZACE on $name"
        publishDir "${params.outDirectory}/${sample.run}/mutect/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.mutect.vcf"), path("${name}.mutect-filtered.vcf")

        output:
        tuple val(name), val(sample), path("${name}.norm.mutect.vcf.gz"), path("${name}.norm.mutect.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo VAFaNORMALIZACE $name

        grep "#" ${name}.mutect-filtered.vcf > head
        grep -v "#" ${name}.mutect-filtered.vcf | grep "PASS"  > x
        cat head x > ${name}.mutect-filtered.PASS.vcf
 

        bcftools +fill-tags ${name}.mutect-filtered.PASS.vcf -Ob -o ${name}.pom2.mutect.bcf -- -t FORMAT/VAF
        bcftools convert -O v -o ${name}.vaf.mutect.vcf ${name}.pom2.mutect.bcf 
         
        bcftools norm -m-both -f ${params.ref}.fa -o ${name}.norm.mutect.vcf ${name}.vaf.mutect.vcf
        bgzip ${name}.norm.mutect.vcf
        tabix ${name}.norm.mutect.vcf.gz

        """
}

process ANOTACE_annovar1 {
       tag "ANOTACE on $name"
       //publishDir "${params.outDirectory}/${sample.run}/mutect/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.mutect.vcf.gz"), path("${name}.norm.mutect.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.norm.mutect.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.mutect.vcf.gz.hg38_multianno.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo ANOTACE $name

        ${params.annovar} -vcfinput ${name}.norm.mutect.vcf.gz ${params.annovardb}  -buildver hg38 -protocol refGeneWithVer,ensGene,1000g2015aug_all,1000g2015aug_eur,exac03nontcga,avsnp150,clinvar_20250721,dbnsfp41c,gnomad41_exome,gnomad41_genome,cosmic70,revel,GTEx_v8_eQTL \
        -operation gx,g,f,f,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -polish -xreffile ${params.gene_fullxref.txt} -arg '-splicing 50 -exonicsplicing',,,,,,,,,,,, --remove
        bgzip ${name}.norm.mutect.vcf.gz.hg38_multianno.vcf
        tabix ${name}.norm.mutect.vcf.gz.hg38_multianno.vcf.gz
        """
}

process VCF2TXT1 {
       tag "VCF2TXT on $name"
       publishDir "${params.outDirectory}/${sample.run}/mutect/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.mutect.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.mutect.vcf.gz.hg38_multianno.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.final.mutect.txt")

        script:
        """
        echo VCF2TXT $name
        source activate gatk4610
        gatk --java-options "-Xmx4g" VariantsToTable -R ${params.ref}.fa  --show-filtered  -V ${name}.norm.mutect.vcf.gz.hg38_multianno.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF SB -GF VAF -F Func.refGeneWithVer -F Gene.refGeneWithVer -F GeneDetail.refGeneWithVer -F ExonicFunc.refGeneWithVer -F AAChange.refGeneWithVer -F 1000g2015aug_all -F 1000g2015aug_eur  -F gnomad41_exome_AF -F gnomad41_exome_AF_nfe -F gnomad41_genome_AF -F gnomad41_genome_AF_nfe -F avsnp150 -F CLNSIG -F REVEL -F SIFT_pred -F MutationTaster_pred -F Gene_full_name.refGeneWithVer -F FATHMM_pred -F PROVEAN_pred -F Function_description.refGeneWithVer -F Disease_description.refGeneWithVer -F Tissue_specificityUniprot.refGeneWithVer -F Expression-egenetics.refGeneWithVer --output ${name}.final.mutect.txt
        """
}

process VARDICT {

         tag "VARDICT on $name"
         publishDir "${params.outDirectory}/${sample.run}/vardict/", mode:'copy'

         input:
         tuple val(name), val(sample), path(bam), path(bai)
         output:
         tuple val(name), val(sample), path("${name}.vardict.vcf")

         script:
         """
         echo VARDICT $name

        source activate vardict
        vardict -G ${params.refindex}.fa -f 0.01 -N ${name} -b ${bam} -c 1 -S 2 -E 3 -g 4 -U ${params.varbed} | Rscript ${params.teststrandbias} | perl ${params.var2vcf_valid} -N ${name} -f 0.01 -A > ${name}.vardict.vcf
       """
}


process VAFaNORMALIZACE2 {
        tag "VAFaNORMALIZACE on $name"
        publishDir "${params.outDirectory}/${sample.run}/vardict/", mode:'copy'

        input:
        tuple val(name), val(sample), path(vardict)

        output:
        tuple val(name), val(sample), path("${name}.norm.vardcit.vcf.gz"), path("${name}.norm.vardict.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo VAFaNORMALIZACE $name
        
        bgzip $vardict
        tabix  ${name}.vardict.vcf.gz
         
        bcftools norm -m-both -o ${name}.bcf1.vcf.gz ${name}.vardict.vcf.gz
        bcftools norm -f ${params.refindex} -o ${name}.norm.vardict.vcf ${name}.bcf1.vcf.gz
        bgzip ${name}.norm.vardict.vcf
        tabix ${name}.norm.vardict.vcf.gz
 
        """
}


process ANOTACE_annovar2 {
       tag "ANOTACE on $name"
       //publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vardict.vcf.gz"), path("${name}.norm.vardict.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.norm.vardict.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.vardict.vcf.gz.hg38_multianno.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo ANOTACE $name

        ${params.annovar} -vcfinput ${name}.norm.vardcit.vcf.gz ${params.annovardb}  -buildver hg38 -protocol refGeneWithVer,ensGene,1000g2015aug_all,1000g2015aug_eur,exac03nontcga,avsnp150,clinvar_20250721,dbnsfp41c,gnomad41_exome,gnomad41_genome,cosmic70,revel,GTEx_v8_eQTL \
        -operation gx,g,f,f,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -polish -xreffile ${params.gene_fullxref.txt} -arg '-splicing 50 -exonicsplicing',,,,,,,,,,,, --remove
        bgzip ${name}.norm.vardict.vcf.gz.hg38_multianno.vcf
        tabix ${name}.norm.vardict.vcf.gz.hg38_multianno.vcf.gz
        """
}

process VCF2TXT2 {
       tag "VCF2TXT on $name"
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vardict.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.vardict.vcf.gz.hg38_multianno.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.final.vardict.txt")

        script:
        """
        echo VCF2TXT $name
        source activate gatk4610
        gatk --java-options "-Xmx4g" VariantsToTable -R ${params.ref}.fa  --show-filtered  -V ${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF SB -F Func.refGeneWithVer -F Gene.refGeneWithVer -F GeneDetail.refGeneWithVer -F ExonicFunc.refGeneWithVer -F AAChange.refGeneWithVer -F 1000g2015aug_all -F 1000g2015aug_eur  -F gnomad41_exome_AF -F gnomad41_exome_AF_nfe -F gnomad41_genome_AF -F gnomad41_genome_AF_nfe -F avsnp150 -F CLNSIG -F REVEL -F SIFT_pred -F MutationTaster_pred -F Gene_full_name.refGeneWithVer -F FATHMM_pred -F PROVEAN_pred -F Function_description.refGeneWithVer -F Disease_description.refGeneWithVer -F Tissue_specificityUniprot.refGeneWithVer -F Expression-egenetics.refGeneWithVer --output ${name}.final.vardict.txt
        """
}


process COVERAGE1 {
          tag "COVERAGE1 on $name"
       publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
         container "staphb/samtools:1.20"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.coveragefin.txt")

        script:
        """
        echo COVERAGE1 $name
        samtools bedcov ${params.varbed} $bam -d 20 > ${name}.COV
        awk '{print \$5/(\$3-\$2)}'  ${name}.COV >  ${name}.COV-mean
        awk '{print (\$6/(\$3-\$2))*100"%"}' ${name}.COV > ${name}-procento-nad-20
        paste ${name}.COV-mean ${name}-procento-nad-20 > vysledek
        echo "chr" "start" "stop" "name" ${name}.COV-mean ${name}-procento-nad-20 > hlavicka
        sed -i 's/ /\t/'g hlavicka
        paste ${params.varbed} vysledek > coverage
        cat hlavicka coverage > ${name}.coveragefin.txt
        sed -i -e "s/\r//g" ${name}.coveragefin.txt
        """
}

process COMBINECOVERAGEMEAN {
    tag "COMBINECOVERAGEMEAN"

    input:
    tuple val(run_name), path(coverage_files)

    publishDir { "${params.outDirectory}/${run_name}/mapped/" }, mode: 'copy'

    output:
    path "coveragemeanALL"


    script:
    """
    i=0
    tmpfiles=""
    for file in ${coverage_files}; do
        awk '{print \$5}' "\$file" > column_\$i.txt
        tmpfiles="\$tmpfiles column_\$i.txt"
        i=\$((i+1))
    done
    echo "chr" "start" "stop" "name" > hlavicka
    sed -i 's/ /\t/'g hlavicka
    cat hlavicka ${params.varbed} > bedshlavickou
    paste  bedshlavickou \$tmpfiles > coveragemeanALL
    """
}

process COMBINECOVERAGEPROCENTA {
    tag "COMBINECOVERAGEPROCENTA"

    input:
    tuple val(run_name), path(coverage_files)

    publishDir { "${params.outDirectory}/${run_name}/mapped/" }, mode: 'copy'

    output:
    path "coverageprocentoALL"


    script:
    """
    i=0
    tmpfiles=""
    for file in ${coverage_files}; do
        awk '{print \$6}' "\$file" > column_\$i.txt
        tmpfiles="\$tmpfiles column_\$i.txt"
        i=\$((i+1))
    done
    echo "chr" "start" "stop" "name" > hlavicka
    sed -i 's/ /\t/'g hlavicka
    cat hlavicka ${params.varbed} > bedshlavickou
    paste  bedshlavickou \$tmpfiles > coverageprocentoALL
    """
}

process LOH {
        tag "LOH"
        publishDir { "${params.outDirectory}/${sample.run}/LOH/fily" }, mode: 'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.chr1.smooth"), path("${name}.chr2.smooth"), path("${name}.chr3.smooth"),path("${name}.chr4.smooth"),path("${name}.chr5.smooth"),path("${name}.chr6.smooth"),path("${name}.chr7.smooth"),path("${name}.chr8.smooth"),path("${name}.chr9.smooth"),path("${name}.chr10.smooth"),path("${name}.chr11.smooth"),path("${name}.chr12.smooth"),path("${name}.chr13.smooth"),path("${name}.chr14.smooth"),path("${name}.chr15.smooth"),path("${name}.chr16.smooth"),path("${name}.chr17.smooth"),path("${name}.chr18.smooth"),path("${name}.chr19.smooth"),path("${name}.chr20.smooth"),path("${name}.chr21.smooth"),path("${name}.chr22.smooth"), path("${name}.chrX.smooth")

        script:
        """
        source activate bcftoolsbgziptabix
        bcftools query   -f '%CHROM\t%POS\t[%GT]\t[%VAF]\n' ${name}.norm.vcf.gz > ${name}.baf.tsv

        awk '{print "chr"\$0}' ${name}.baf.tsv > ${name}.baf-chr.tsv

       grep -w chr1 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr1.smooth
       grep -w chr2 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr2.smooth
       grep -w chr3 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr3.smooth
       grep -w chr4 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr4.smooth
       grep -w chr5 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr5.smooth
       grep -w chr6 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr6.smooth
       grep -w chr7 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr7.smooth
       grep -w chr8 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr8.smooth
       grep -w chr9 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr9.smooth
       grep -w chr10 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr10.smooth
       grep -w chr11 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr11.smooth
       grep -w chr12 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr12.smooth
       grep -w chr13 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr13.smooth
       grep -w chr14 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr14.smooth
       grep -w chr15 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr15.smooth
       grep -w chr16 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr16.smooth
       grep -w chr17 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr17.smooth
       grep -w chr18 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr18.smooth
       grep -w chr19 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr19.smooth
       grep -w chr20 ${name}.baf-chr.tsv  |awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr20.smooth
       grep -w chr21 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr21.smooth
       grep -w chr22 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr22.smooth
       grep -w chrX ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chrX.smooth



       """
}

process GNUPLOT {
        tag "GNUPLOT"
        publishDir { "${params.outDirectory}/${sample.run}/LOH/" }, mode: 'copy'
        container "xreblova/gnuplot:1"

        input:
        tuple val(name), val(sample), path("${name}.chr1.smooth"), path("${name}.chr2.smooth"), path("${name}.chr3.smooth"),path("${name}.chr4.smooth"),path("${name}.chr5.smooth"),path("${name}.chr6.smooth"),path("${name}.chr7.smooth"),path("${name}.chr8.smooth"),path("${name}.chr9.smooth"),path("${name}.chr10.smooth"),path("${name}.chr11.smooth"),path("${name}.chr12.smooth"),path("${name}.chr13.smooth"),path("${name}.chr14.smooth"),path("${name}.chr15.smooth"),path("${name}.chr16.smooth"),path("${name}.chr17.smooth"),path("${name}.chr18.smooth"),path("${name}.chr19.smooth"),path("${name}.chr20.smooth"),path("${name}.chr21.smooth"),path("${name}.chr22.smooth"), path("${name}.chrX.smooth")

        output:
        tuple val(name), val(sample), path("${name}.loh.png")

        script:
        """
        echo 'set term png size 2000,900' > ${name}.gnuplot.script
        echo 'set output "${name}.loh.png"' >> ${name}.gnuplot.script
        echo 'set xtics font "Arial,9"' >> ${name}.gnuplot.script
        echo 'set multiplot layout 6,4' >> ${name}.gnuplot.script

        echo 'plot "${name}.chr1.smooth" using 1:2 w l title "CHR1"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr2.smooth" using 1:2 w l title "CHR2"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr3.smooth" using 1:2 w l title "CHR3"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr4.smooth" using 1:2 w l title "CHR4"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr5.smooth" using 1:2 w l title "CHR5"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr6.smooth" using 1:2 w l title "CHR6"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr7.smooth" using 1:2 w l title "CHR7"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr8.smooth" using 1:2 w l title "CHR8"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr9.smooth" using 1:2 w l title "CHR9"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr10.smooth" using 1:2 w l title "CHR10"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr11.smooth" using 1:2 w l title "CHR11"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr12.smooth" using 1:2 w l title "CHR12"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr13.smooth" using 1:2 w l title "CHR13"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr14.smooth" using 1:2 w l title "CHR14"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr15.smooth" using 1:2 w l title "CHR15"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr16.smooth" using 1:2 w l title "CHR16"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr17.smooth" using 1:2 w l title "CHR17"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr18.smooth" using 1:2 w l title "CHR18"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr19.smooth" using 1:2 w l title "CHR19"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr20.smooth" using 1:2 w l title "CHR20"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr21.smooth" using 1:2 w l title "CHR21"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr22.smooth" using 1:2 w l title "CHR22"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chrX.smooth" using 1:2 w l title "CHRX"' >> ${name}.gnuplot.script

        echo 'unset multiplot' >> ${name}.gnuplot.script

gnuplot ${name}.gnuplot.script
      """
}



workflow {
        rawfastq = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    .splitCsv(header: true)
    .map { row ->
        def baseDir = new File("${params.baseDir}")
        def runDir = baseDir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(row.run)
            }
        })[0] //get the real folderName that has prepended date

        def fileR1 = file("${runDir}/processed_fastq/${row.name}_R1.fastq.gz", checkIfExists: true)
        def fileR2 = file("${runDir}/processed_fastq/${row.name}_R2.fastq.gz", checkIfExists: true)

                def meta = [name: row.name, run: row.run]
        [
            meta.name,
            meta,
            fileR1,
            fileR2,
                ]
    }
     . view()

// kontrola1 = FASTQC1(rawfastq)

 //Group the FASTQC1 output files by run name
 //   grouped_fastqc_files = kontrola1
 //       .map { name, sample, r1html, r2html, r1zip, r2zip ->
 //           tuple(sample.run, [r1html, r2html, r1zip, r2zip])
 //       }
 //       .groupTuple()           // groups by run name: [run_name: [[files], [files], ...]]
 //       .map { run_name, lists_of_files ->
 //           def all_files = lists_of_files.flatten()
 //           tuple(run_name, all_files)
 //       }

// MULTIQC1(grouped_fastqc_files)


aligned = ALIGN(rawfastq)
// kontrolabamu1 = FLAGSTAT(aligned)
// kontrolabamu2 = QUALIMAP(aligned)

varcalling = GATK(aligned)
normalizovany = VAFaNORMALIZACE1(varcalling)
anotovany = ANOTACE_annovar1(normalizovany)
anotovanyfin = VCF2TXT1(anotovany)

// varcalling2 = VARDICT(aligned)
// normalizovany2 = VAFaNORMALIZACE2(varcalling2)
// anotovany2 = ANOTACE_annovar2(normalizovany2)
// anotovanyfin2 = VCF2TXT2(anotovany2)

coverage_results = COVERAGE1(aligned)
coverage_files_collected = coverage_results
    .map { name, sample, f -> tuple(sample.run, file(f)) }
    .groupTuple() // groups by sample.run automatically!
finalcoverage = COMBINECOVERAGEMEAN(coverage_files_collected)
finalprocenta = COMBINECOVERAGEPROCENTA(coverage_files_collected)

loh = LOH(normalizovany)
graf = GNUPLOT(loh)


}
