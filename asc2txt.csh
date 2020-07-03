for g in subj*.asc

do

echo $g;

name=$(awk 'BEGIN{f=substr(ARGV[1],1,9);print f}' $g)

grep '^[0-9]' $g | awk '{print $1"\t"$2"\t"$3"\t"$4}' | awk -v OFS="\t" '$2=="."{$2=0} $3=="."{$3=0} {$1=$1}1' > $name.txt 

done
