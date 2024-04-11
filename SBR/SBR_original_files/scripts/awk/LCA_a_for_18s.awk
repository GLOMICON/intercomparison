# Set field separator as tab or ;
BEGIN{
	FS="\t"
	print "amplicon","identity","sequence","references","taxonomy"
}

NR==1  {a=$1;b=$2;c=$3;d=$4;n=split($5,e,",")}
$1==a && NR >1 {
	d = d","$4
	split($5,f,",")
	for (i=1;i<=n;i++) {
		if (e[i]!=f[i]){
			e[i]="*"
		}
	split("",tmp)
	}
}
$1!=a && NR>1 {
	output=sep=""
	for (i=1;i<=n;i++) {
		output=output sep e[i]
		sep="|"
	}
	print a,b,c,d,output
	a=$1;b=$2;c=$3;d=$4
	n=split($5,e,",")
}
END{
	output=sep=""
	for (i=1;i<=n;i++) {
		output=output sep e[i]
		sep="|"
	}
	print a,b,c,d,output
}
