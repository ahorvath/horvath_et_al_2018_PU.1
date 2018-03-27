if [[ $# -ne 1 ]]; then
        echo -e "USAGE: $0 <FILE>\n"
        echo -e "\tEXAMPLE: $0 <FILE>"
	echo -e "\tINPUT:"
	echo -e "\t\t 300 STAT6"
	echo -e "\t\t 150 H4ac"
	echo -e "\t\t 200 STAT6,H4ac"
        exit 1
fi

FILE=$1
awk 'BEGIN{row=1}{OFS="\t"; for (i=1;i<=$1;++i) {maxind=split($2, ids, ",");  for (j=1;j<=maxind;++j) {print "elem"row,ids[j];} ++row;}}' $FILE
