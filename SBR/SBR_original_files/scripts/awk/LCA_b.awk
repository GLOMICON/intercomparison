# Set field separator as tab
BEGIN{
	FS="\t"
	OFS="\t"
	print "amplicon","identity","seq_id","taxonomy"
}

# Extract taxo and id
{split($3,desc,"|")}

# Split taxo and compare last items until finding the lowest common ancestor
NR==1  {a=$1;b=$2;c=desc[1];n=split(desc[2],d,"_")}
$1==a && NR >1 {
	c = c","desc[1]
	split(desc[2],e,"_")
	for (i=1;i<=n;i++) {
		if (d[i]!=e[i]){
			d[i]="*"
		}
	}
}
$1!=a && NR>1 {
	output=sep=""
	for (i=1;i<=n;i++) {
		output=output sep d[i]
		sep="_" # here to change output taxo separator
	}
	print a,b,c,output
	a=$1;b=$2;c=desc[1]
	n=split(desc[2],d,"_")
}
END{
	output=sep=""
	for (i=1;i<=n;i++) {
		output=output sep d[i]
		sep="_" # here to change output taxo separator
	}
	print a,b,c,output
}
