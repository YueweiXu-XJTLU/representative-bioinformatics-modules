echo "===用map到已知基因上的高TPM的点位, 找与RAW一致的点位==="
awk -F'\t' 'NR>1 && $1~/^gene:/{
  printf("%.3f\t%s\t%s\tchr%s:%s-%s\tstrand=%s\n",$9,$1,$2,$3,$5,$6,$4)
}' C1.txt | sort -k1,1nr | head -n 10



echo "===用同一个gene出现多个转录本(更可能是异构体), 找与RAW不一致的点位==="
awk -F'\t' '$3=="transcript"{
  attr=$9; gid=attr; sub(/.*gene_id "/,"",gid); sub(/".*/,"",gid); cnt[gid]++
}
END{for(g in cnt) if(cnt[g]>=2) printf("%d\t%s\n",cnt[g],g)}' TAIR10_re.gtf \
| sort -k1,1nr | head -n 10



echo "===用RECONSTRUCTION中未map到抑制基因上, 但是表达量高的条目(STRG), 找RAW未出现的点位==="
awk -F'\t' 'NR>1 && $1~/^STRG\./{
  printf("%.3f\t%s\tchr%s:%s-%s\tstrand=%s\n",$9,$1,$3,$5,$6,$4)
}' C1.txt | sort -k1,1nr | head -n 10