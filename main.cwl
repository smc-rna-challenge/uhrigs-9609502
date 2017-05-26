arguments: []
baseCommand: [bash, run_arriba.sh, '12', '50000', /opt/arriba]
class: CommandLineTool
cwlVersion: sbg:draft-2
dct:creator: {'@id': 'http://orcid.org/0000-0002-7681-6415', 'foaf:mbox': uhrigs@synapse.org,
  'foaf:name': uhrigs}
description: ''
doc: 'SMC-RNA challenge fusion detection submission

  '
hints:
- {class: 'sbg:CPURequirement', value: 12}
- {class: 'sbg:MemRequirement', value: 50000}
- {class: DockerRequirement, dockerImageId: '', dockerPull: 'quay.io/smc-rna-challenge/uhrigs-9609502-smc-rna-challenge-5:1.0.0'}
id: https://cgc-api.sbgenomics.com/v2/apps/uhrigs/smc-rna-challenge-5/arriba-0-9-sensitive-pc-nocgc/2/raw/
inputs:
- description: path to tar archive containing STAR index
  id: '#STAR_INDEX_TAR'
  inputBinding: {position: 1, 'sbg:cmdInclude': true, separate: true}
  label: STAR_INDEX_TAR
  sbg:fileTypes: TAR,TAR.GZ
  type: [File]
- description: path to reference annotation in GTF format (may be gzip-compressed)
  id: '#REFERENCE_GTF'
  inputBinding: {position: 2, 'sbg:cmdInclude': true, separate: true}
  label: REFERENCE_GTF
  sbg:fileTypes: GTF.GZ,GTF
  type: [File]
- description: path to reference genome in FastA format (may be gzip-compressed)
  id: '#REFERENCE_GENOME'
  inputBinding: {position: 3, 'sbg:cmdInclude': true, separate: true}
  label: REFERENCE_GENOME
  sbg:fileTypes: FA.GZ,FA
  type: [File]
- description: path to gzip-compressed FastQ file containing first mates
  id: '#TUMOR_FASTQ_1'
  inputBinding: {position: 4, 'sbg:cmdInclude': true, separate: true}
  label: TUMOR_FASTQ_1
  sbg:fileTypes: FQ.GZ,FASTQ.GZ
  type: [File]
- description: path to gzip-compressed FastQ file containing second mates
  id: '#TUMOR_FASTQ_2'
  inputBinding: {position: 5, 'sbg:cmdInclude': true, separate: true}
  label: TUMOR_FASTQ_2
  sbg:fileTypes: FQ.GZ,FASTQ.GZ
  type: [File]
label: arriba-0.9-sensitive-pc-nocgc
outputs:
- description: predicted fusions in BEDPE format
  id: '#OUTPUT'
  label: OUTPUT
  outputBinding: {glob: fusions.bedpe}
  sbg:fileTypes: BEDPE
  type: [File]
requirements:
- class: CreateFileRequirement
  fileDef:
  - {fileContent: "#!/bin/bash\n\nset -o pipefail\n\nif [ $# -ne 8 ]; then\n\techo\
      \ \"Usage: $(basename $0) threads memory tools_dir STAR_index.tar.gz annotation.gtf.gz\
      \ assembly.fa.gz read1.fastq.gz read2.fastq.gz\" 1>&2\n\texit 1\nfi\n\nset -x\n\
      \n# fetch arguments\nTHREADS=\"$1\"\nMEMORY=\"$2\"\nTOOLS_DIR=\"$3\"\nSTAR_INDEX=\"\
      $4\"\nREFERENCE_GTF=\"$5\"\nREFERENCE_GENOME=\"$6\"\nREAD1=\"$7\"\nREAD2=\"\
      $8\"\nFUSIONS_OUT=\"fusions\"\n\nfunction autounzip() {\n\tif [[ \"$1\" =~ \\\
      .gz$ ]]; then\n\t\t\"$TOOLS_DIR/pigz\" -d -c \"$1\"\n\telse\n\t\tcat \"$1\"\n\
      \tfi\n}\n\n# extract STAR index\nmkdir STAR_index || exit 1\nautounzip \"$STAR_INDEX\"\
      \ | tar -x -C STAR_index -f - & PID_STAR_INDEX=$!\n\n# protein-coding genes\
      \ from GTF\nautounzip \"$REFERENCE_GTF\" |\n\"$TOOLS_DIR/awk\" -v FS='\\t' '\n\
      \t$2==\"protein_coding\" {\n\t\tgene_id = $9;\n\t\tsub(/.*gene_id \"/, \"\"\
      , gene_id);\n\t\tsub(/\".*/, \"\", gene_id);\n\t\tif ($3 == \"gene\")\n\t\t\t\
      protein_coding[gene_id]++;\n\t\tif (gene_id in protein_coding)\n\t\t\tprint;\n\
      \t}\n' > protein_coding.gtf & PID_REFERENCE_GTF=$!\n\nwait $PID_STAR_INDEX ||\
      \ exit 1\nwait $PID_REFERENCE_GTF || exit 1\n\n# run alignment\n\"$TOOLS_DIR/STAR-2.5.3a\"\
      \ \\\n\t--runThreadN \"$THREADS\" \\\n\t--genomeDir \"$(find -name SAindex -printf\
      \ %h)\" --genomeLoad NoSharedMemory \\\n\t--readFilesIn \"$READ1\" \"$READ2\"\
      \ --readFilesCommand zcat \\\n\t--outStd BAM_Unsorted --outSAMtype BAM Unsorted\
      \ SortedByCoordinate \\\n\t--outFilterMultimapNmax 1 --outFilterMismatchNmax\
      \ 3 --outFilterMismatchNoverLmax 0.3 \\\n\t--alignIntronMax 500000 --alignMatesGapMax\
      \ 500000 \\\n\t--chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimScoreMin\
      \ 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation\
      \ 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimMainSegmentMultNmax\
      \ 10 \\\n\t--limitBAMsortRAM ${MEMORY}000000 |\n\"$TOOLS_DIR/extract_read-through_fusions-0.9\"\
      \ -g protein_coding.gtf -G \"gene_name=gene_name gene_id=gene_id transcript_id=transcript_id\
      \ gene_status=gene_biotype status_KNOWN=protein_coding gene_type=gene_biotype\
      \ type_protein_coding=protein_coding feature_exon=exon feature_UTR=UTR feature_gene=gene\"\
      \ > read_through.bam || exit 1\n\n# index normal alignments\nmv Aligned.sortedByCoord.out.bam\
      \ rna.bam || exit 1\n\"$TOOLS_DIR/samtools\" index rna.bam & PID_SAMTOOLS_INDEX=$!\n\
      \n# convert chimeric alignments from SAM to BAM\n\"$TOOLS_DIR/samtools\" view\
      \ -Sb Chimeric.out.sam > chimeric.bam & PID_SAMTOOLS_VIEW=$!\n\n# index reference\
      \ genome\nautounzip \"$REFERENCE_GENOME\" > assembly.fa || exit 1\n\"$TOOLS_DIR/samtools\"\
      \ faidx assembly.fa || exit 1\n\nwait $PID_SAMTOOLS_INDEX || exit 1\nwait $PID_SAMTOOLS_VIEW\
      \ || exit 1\n\n# call arriba\n\"$TOOLS_DIR/arriba-0.9\" \\\n\t-c chimeric.bam\
      \ \\\n\t-r read_through.bam \\\n\t-x rna.bam \\\n\t-o \"$FUSIONS_OUT.tsv\" \\\
      \n\t-O \"$FUSIONS_OUT.discarded.tsv\" \\\n\t-a assembly.fa \\\n\t-g protein_coding.gtf\
      \ \\\n\t-G \"gene_name=gene_name gene_id=gene_id transcript_id=transcript_id\
      \ gene_status=gene_biotype status_KNOWN=protein_coding gene_type=gene_biotype\
      \ type_protein_coding=protein_coding feature_exon=exon feature_UTR=UTR feature_gene=gene\"\
      \ \\\n\t-b \"$TOOLS_DIR/database/blacklist_GRCh37.75_2017-03-05.tsv.gz\" \\\n\
      \t-f \"spliced,many_spliced\" || exit 1\n\n# SMC-RNA Challenge-specific filters\n\
      \"$TOOLS_DIR/awk\" '\n\t($1 ~ /^RP(L|S)[0-9]/ || $2 ~ /^RP(L|S)[0-9]/) &&\n\t\
      $12 > 0 && $13 > 0 && $14 > 0 &&\n\t$7 == \"splice-site\" && $8 == \"splice-site\"\
      \ &&\n\t!($9 ~ /deletion.read.through/)\n' \"$FUSIONS_OUT.discarded.tsv\" >>\
      \ \"$FUSIONS_OUT.tsv\"\n\"$TOOLS_DIR/awk\" '\n\t{\n\t\thomolog1 = $1; homolog2\
      \ = $2;\n\t\tsub(/-.*/, \"\", homolog1); sub(/-.*/, \"\", homolog2);\n\t\tsub(/([A-Z][0-9]+)*[A-Z]?$/,\
      \ \"\", homolog1);\n\t\tsub(/([A-Z][0-9]+)*[A-Z]?$/, \"\", homolog2);\n\t}\n\
      \t!duplicate[$5,$6]++ && !duplicate[$6,$5]++ &&\n\t!duplicate[$1,$2]++ && !duplicate[$2,$1]++\
      \ &&\n\t!duplicate[homolog1,homolog2]++ && !duplicate[homolog2,homolog1]++ &&\n\
      \t$15 != \"low\" &&\n\t($7 == \"splice-site\" && $8 == \"splice-site\" && $3\
      \ ~ /\\+\\/\\+|-\\/-/ && $4 ~ /\\+\\/\\+|-\\/-/ || $7 == \"exon\" && $8 == \"\
      exon\" && $12+$13 == 0 && $3 ~ /\\+\\/\\+|-\\/-|.\\/\\./ && $4 ~ /\\+\\/\\+|-\\\
      /-|.\\/\\./) &&\n\t!($9 ~ /5.-5.|3.-3./) &&\n\thomolog1 != homolog2 &&\n\t!($1\
      \ ~ /^HIST[0-9]/ && $2 ~ /^HIST[0-9]/) &&\n\t!($1 ~ /^COL[0-9]/ && $2 ~ /^COL[0-9]/)\
      \ &&\n\t!($1 ~ /^NBPF[0-9]/ && $2 ~ /^NBPF[0-9]/) &&\n\t!($1 ~ /^AC[0-9]+\\\
      .[0-9]+$/ && $2 ~ /^AC[0-9]+\\.[0-9]+$/) &&\n\t!($1 ~ /PDE4DIP|MALAT1/ || $2\
      \ ~/PDE4DIP|MALAT1/) &&\n\t$1 != $2\n' \"$FUSIONS_OUT.tsv\" > \"$FUSIONS_OUT.filtered.tsv\"\
      \ || exit 1\n\n# convert to BEDPE\n\"$TOOLS_DIR/awk\" '\nBEGIN{ OFS=\"\\t\"\
      \ }\n{\n\tchromosome1 = $5; position1 = $5;\n\tsub(/:.*/, \"\", chromosome1);\
      \ sub(/.*:/, \"\", position1);\n\tchromosome2 = $6; position2 = $6;\n\tsub(/:.*/,\
      \ \"\", chromosome2); sub(/.*:/, \"\", position2);\n\tif ($10 == \"downstream\"\
      ) position1--;\n\tif ($11 == \"downstream\") position2--;\n\tstrand1 = $3; if\
      \ (strand1 ~ /\\.$/) { sub(/\\/.*/, \"\", strand1) } else { sub(/.*\\//, \"\"\
      , strand1) }\n\tstrand2 = $4; if (strand2 ~ /\\.$/) { sub(/\\/.*/, \"\", strand2)\
      \ } else { sub(/.*\\//, \"\", strand2) }\n\tfilters = $18;\n\tgsub(/a-z_()/,\
      \ \"\", filters);\n\tsplit(filters, filtered_reads, \",\");\n\treads = $12+$13+$14;\n\
      \tfor (i in filtered_reads) reads += int(filtered_reads[i]);\n\tprint chromosome1,\
      \ position1-1, position1, chromosome2, position2-1, position2, $1\">>\"$2, $15,\
      \ strand1, strand2, reads;\n}' \"$FUSIONS_OUT.filtered.tsv\" > \"$FUSIONS_OUT.bedpe\"\
      \ || exit 1", filename: run_arriba.sh}
sbg:appVersion: ['sbg:draft-2']
sbg:cmdPreview: bash run_arriba.sh 12 50000 /opt/arriba  /path/to/STAR_index_dir.ext  /path/to/annotation.ext  /path/to/assembly.ext  /path/to/fastq1.ext  /path/to/fastq2.ext
sbg:contributors: [uhrigs]
sbg:createdBy: uhrigs
sbg:createdOn: 1494578126
sbg:id: uhrigs/smc-rna-challenge-5/arriba-0-9-sensitive-pc-nocgc/2
sbg:image_url: null
sbg:job:
  allocatedResources: {cpu: 12, mem: 50000}
  inputs:
    REFERENCE_GENOME:
      class: File
      path: /path/to/assembly.ext
      secondaryFiles:
      - {path: .fai}
      size: 0
    REFERENCE_GTF:
      class: File
      path: /path/to/annotation.ext
      secondaryFiles: []
      size: 0
    STAR_INDEX_TAR:
      class: File
      path: /path/to/STAR_index_dir.ext
      secondaryFiles: []
      size: 0
    TUMOR_FASTQ_1:
      class: File
      path: /path/to/fastq1.ext
      secondaryFiles: []
      size: 0
    TUMOR_FASTQ_2:
      class: File
      path: /path/to/fastq2.ext
      secondaryFiles: []
      size: 0
sbg:latestRevision: 2
sbg:license: MIT License
sbg:modifiedBy: uhrigs
sbg:modifiedOn: 1494578842
sbg:project: uhrigs/smc-rna-challenge-5
sbg:projectName: smc-rna-challenge-5
sbg:revision: 2
sbg:revisionsInfo:
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1494578126, 'sbg:revision': 0, 'sbg:revisionNotes': Copy
    of uhrigs/smc-rna-challenge-4/arriba-0-8/2}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1494578817, 'sbg:revision': 1, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1494578842, 'sbg:revision': 2, 'sbg:revisionNotes': null}
sbg:sbgMaintained: false
sbg:toolAuthor: Sebastian Uhrig, DKFZ Heidelberg
sbg:toolkitVersion: '0.9'
sbg:validationErrors: []
stdin: ''
stdout: ''
successCodes: []
temporaryFailCodes: []
