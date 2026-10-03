#!/usr/bin/env bash
GTF=./TAIR10_re.gtf
ABUND=./C1.txt

echo "== GTF-level (re-assembled transcripts) =="
# 1) 转录本数 & 外显子数
awk '$3=="transcript"{t++} $3=="exon"{e++} END{printf("transcripts\t%d\nexons\t%d\n",t,e)}' "$GTF"

# 2) 每转录本外显子中位数（BSD awk 兼容：从第9列attributes中提取transcript_id）
awk -F'\t' '
  $3=="exon"{
    tid=$9
    if (tid ~ /transcript_id "/) {
      sub(/.*transcript_id "/,"",tid)
      sub(/".*/,"",tid)
      ex[tid]++
    }
  }
  END{
    for (t in ex) print ex[t]
  }
' "$GTF" \
| sort -n \
| awk '{a[NR]=$1}
       END{
         if(NR==0){print "median_exons_per_transcript\tNA"; exit}
         if(NR%2)  m=a[(NR+1)/2];
         else      m=(a[NR/2]+a[NR/2+1])/2;
         printf("median_exons_per_transcript\t%s\n", m)
       }'

# 3) 转录本跨度中位数（end-start+1）
awk '$3=="transcript"{print $5-$4+1}' "$GTF" \
| sort -n \
| awk '{a[NR]=$1}
       END{
         if(NR==0){print "median_transcript_span_bp\tNA"; exit}
         if(NR%2)  m=a[(NR+1)/2];
         else      m=(a[NR/2]+a[NR/2+1])/2;
         printf("median_transcript_span_bp\t%s\n", m)
       }'

echo -e "\n== Locus-level abundance from C1.txt =="
# 4) 已知基因 vs STRG 计数 + TPM 总和
awk -F'\t' 'NR>1{ if($1~/^gene:/) known++; else if($1~/^STRG\./) strg++; sum+=$9 }
END{printf("known_genes_rows\t%d\nstrg_loci_rows\t%d\nTPM_sum\t%.3f\n",known,strg,sum)}' "$ABUND"

# 5) TPM Top 5 已知基因
echo -e "Top5_known_by_TPM\tTPM\tGeneID\tGeneName"
awk -F'\t' 'NR>1 && $1~/^gene:/{printf("%.6f\t%s\t%s\n",$9,$1,$2)}' "$ABUND" \
| sort -k1,1nr | head -n 5 \
| awk -F'\t' '{printf("\t%.3f\t%s\t%s\n",$1,$2,$3)}'