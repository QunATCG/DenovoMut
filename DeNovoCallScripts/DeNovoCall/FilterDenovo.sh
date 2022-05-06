# input
INPUT=$1
OUTPUT=$2
# help
if [ -z $1 ] || [ -z $2 ] ;then
                echo -e "\033[32m{Useage: $0 INPUT_vcf OUTPUT_Dir}\033[0m"
                echo -e "\033[32m  Author: LiQun (liqun95@163.com) \033[0m"
                exit
fi
# filter
cat << 'EOF'
1)filter....
	FMT/GQ>40 && INFO/DP>10 && INFO/MQ>40 && INFO/QD>2 && INFO/FS<30 && INFO/VQSLOD>-10
EOF
TMPFILE=${RANDOM}.$$.vcf
TMPOUTVCF=${RANDOM}.$$.out.vcf
cat $INPUT | grep "#" >> $TMPFILE
cat $INPUT | grep "PASS" | grep -v "#" >> $TMPFILE
cat $TMPFILE | grep "#" >> $TMPOUTVCF
bcftools query -i 'FMT/GQ>40 && INFO/DP>10 && INFO/MQ>40 && INFO/QD>2 && INFO/FS<30 && INFO/VQSLOD>-10' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT[\t%GT:%AD:%GQ:%PL]\n'  $TMPFILE >> $TMPOUTVCF
rm $TMPFILE
# generate ped file
cat << 'EOF'
2)generate ped file....
EOF
TMPSAMPLE=${RANDOM}.$$.sample
TMPPED=${RANDOM}.$$.ped
echo -e "FID\tIID\tRID\tLOC" >> $TMPPED
bcftools query -l $TMPOUTVCF > $TMPSAMPLE
COUNT=0
while read line
do
        FID=`echo $line | awk -F "[_-]" '{ print $1 }'`
        RID=`echo $line | awk -F "[_-]" '{ print $2 }'`
        IID=`echo $line | awk -F "[_-]" '{ print $3 }'`
        COUNT=$((${COUNT}+1))
        echo -e "$FID\t$IID\t$RID\t$COUNT" >> $TMPPED    
done < $TMPSAMPLE
rm $TMPSAMPLE
# extract
cat << 'EOF'
3)extract De novo mutations by FID....
	filter mutations with DP > 10
EOF
python /home/liqun/denovo/extractDenovo.py --vcf $TMPOUTVCF --ped $TMPPED --out $OUTPUT
rm $TMPOUTVCF
#rm $TMPPED
cat << 'EOF'
3)finished
	please check output files
EOF
